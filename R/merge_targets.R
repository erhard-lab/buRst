# Merge two buRst .targets directories by cells.
# Supported: identical genes/features, disjoint cells.

.targets_read_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

.targets_write_lines <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(as.character(x), con, useBytes = TRUE)
  invisible(path)
}

.targets_read_mtx <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  Matrix::readMM(con)
}

.targets_write_mtx <- function(x, path) {
  tmp <- tempfile(fileext = ".mtx")
  on.exit(unlink(tmp), add = TRUE)
  Matrix::writeMM(x, tmp)

  input <- file(tmp, "rb")
  output <- gzfile(path, "wb")
  on.exit(close(input), add = TRUE)
  on.exit(close(output), add = TRUE)

  repeat {
    z <- readBin(input, "raw", n = 1024L * 1024L)
    if (!length(z)) break
    writeBin(z, output)
  }
  invisible(path)
}

.targets_find_h5 <- function(dir, h5_filename = NULL) {
  if (!is.null(h5_filename)) {
    p <- file.path(dir, h5_filename)
    if (!file.exists(p)) stop("Missing HDF5: ", p, call. = FALSE)
    return(p)
  }

  p <- list.files(dir, "\\.h5$", full.names = TRUE, ignore.case = TRUE)
  if (length(p) != 1L) {
    stop(
      "Expected exactly one .h5 file in ", dir,
      "; found ", length(p), ". Supply h5_filename.",
      call. = FALSE
    )
  }
  p
}

.targets_write_dataset <- function(
    h5,
    path,
    value,
    compression_level = 4L,
    chunk_target = 1e6
) {
  value <- as.vector(value)
  if (!length(value)) {
    # In this format a fully empty profile is represented by no profile values.
    # Store a one-element empty marker only when rhdf5 cannot create length zero.
    value <- integer(0)
    rhdf5::h5write(value, h5, path)
    return(invisible(path))
  }

  bytes <- if (is.integer(value)) 4L else 8L
  chunk <- max(1L, min(length(value), as.integer(chunk_target / bytes)))

  rhdf5::h5createDataset(
    file = h5,
    dataset = path,
    dims = length(value),
    chunk = chunk,
    level = as.integer(compression_level),
    storage.mode = if (is.integer(value)) "integer" else "double"
  )
  rhdf5::h5write(value, h5, path)
  invisible(path)
}

.targets_merge_meta <- function(a, b, n1, n2, gene) {
  if (is.null(a) && is.null(b)) return(NULL)
  if (is.null(a) || is.null(b)) {
    stop("/genes_meta/", gene, " exists in only one input.", call. = FALSE)
  }
  if (identical(a, b)) return(a)

  # Per-cell metadata: concatenate directly.
  if (length(a) == n1 && length(b) == n2) {
    return(c(a, b))
  }

  # Ragged-array offsets: one cumulative offset per cell boundary, hence
  # n_cells + 1 values. The second offset vector must be shifted by the
  # terminal offset of the first input; omit b[1] to avoid duplicating the
  # shared boundary.
  if (length(a) == n1 + 1L && length(b) == n2 + 1L) {
    a_num <- as.numeric(a)
    b_num <- as.numeric(b)

    valid_offsets <- function(x) {
      length(x) >= 1L &&
        all(is.finite(x)) &&
        x[1L] == 0 &&
        all(diff(x) >= 0)
    }

    if (valid_offsets(a_num) && valid_offsets(b_num)) {
      merged <- c(
        a_num,
        a_num[length(a_num)] + b_num[-1L]
      )

      if (is.integer(a) && is.integer(b)) {
        merged <- as.integer(merged)
      }

      return(merged)
    }
  }

  da <- dim(a)
  db <- dim(b)
  if (!is.null(da) && !is.null(db) && length(da) == length(db)) {
    different <- which(da != db)
    if (length(different) == 1L) {
      margin <- different[[1L]]
      if (da[margin] == n1 && db[margin] == n2) {
        if (margin == 1L) return(rbind(a, b))
        if (margin == 2L) return(cbind(a, b))
      }
    }
  }

  stop(
    "Cannot safely merge /genes_meta/", gene,
    ": unsupported metadata layout. ",
    "Input lengths are ", length(a), " and ", length(b),
    "; expected either n_cells (", n1, ", ", n2,
    ") or n_cells + 1 offsets (", n1 + 1L, ", ", n2 + 1L, ").",
    call. = FALSE
  )
}

#' Merge two buRst .targets directories
#'
#' This performs a cell-wise merge. Both inputs must have identical genes and
#' feature order, while their cell names must not overlap.
#'
#' @param targets1 First .targets directory.
#' @param targets2 Second .targets directory.
#' @param output_targets New output .targets directory.
#' @param h5_filename Optional HDF5 filename. If NULL, each input must contain
#'   exactly one .h5 file.
#' @param overwrite Replace an existing output directory.
#' @param compression_level HDF5 compression level.
#' @param chunk_target Approximate HDF5 chunk size in bytes.
#' @param verbose_every Progress interval in genes; zero disables messages.
#' @return Invisibly, a manifest describing the merged directory.
#' @seealso `SimulateReads()`, `write_10x_reference()`
#'
#' @examples
#' \dontrun{
#' # Inputs must share their gene order and have distinct cell names.
#' manifest <- merge_targets_dirs(
#'   targets1 = "reference.targets", targets2 = "perturbed.targets",
#'   output_targets = "combined.targets",
#'   h5_filename = "4sU.Binom.discrete.h5", overwrite = FALSE
#' )
#' manifest
#' }
#'
#' @export
merge_targets_dirs <- function(
    targets1,
    targets2,
    output_targets,
    h5_filename = NULL,
    overwrite = FALSE,
    compression_level = 4L,
    chunk_target = 1e6,
    verbose_every = 100L
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix.")
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Need rhdf5.")

  targets1 <- normalizePath(targets1, mustWork = TRUE)
  targets2 <- normalizePath(targets2, mustWork = TRUE)
  output_targets <- sub("[/\\\\]+$", "", output_targets)

  required <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
  for (d in c(targets1, targets2)) {
    missing <- required[!file.exists(file.path(d, required))]
    if (length(missing)) {
      stop("Missing in ", d, ": ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  h1 <- .targets_find_h5(targets1, h5_filename)
  h2 <- .targets_find_h5(targets2, h5_filename)
  output_h5_name <- if (is.null(h5_filename)) basename(h1) else h5_filename

  features1 <- .targets_read_lines(file.path(targets1, "features.tsv.gz"))
  features2 <- .targets_read_lines(file.path(targets2, "features.tsv.gz"))
  if (!identical(features1, features2)) {
    stop("features.tsv.gz files are not identical.", call. = FALSE)
  }

  barcodes1 <- .targets_read_lines(file.path(targets1, "barcodes.tsv.gz"))
  barcodes2 <- .targets_read_lines(file.path(targets2, "barcodes.tsv.gz"))
  overlap <- intersect(barcodes1, barcodes2)
  if (length(overlap)) {
    stop(
      "Overlapping cell names: ",
      paste(utils::head(overlap, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  m1 <- .targets_read_mtx(file.path(targets1, "matrix.mtx.gz"))
  m2 <- .targets_read_mtx(file.path(targets2, "matrix.mtx.gz"))
  if (nrow(m1) != length(features1) || nrow(m2) != length(features2)) {
    stop("Matrix rows do not match features.", call. = FALSE)
  }
  if (ncol(m1) != length(barcodes1) || ncol(m2) != length(barcodes2)) {
    stop("Matrix columns do not match barcodes.", call. = FALSE)
  }

  genes1 <- as.character(rhdf5::h5read(h1, "/meta/gene_names"))
  genes2 <- as.character(rhdf5::h5read(h2, "/meta/gene_names"))
  cells1 <- as.character(rhdf5::h5read(h1, "/meta/cell_names"))
  cells2 <- as.character(rhdf5::h5read(h2, "/meta/cell_names"))

  if (!identical(genes1, genes2)) stop("HDF5 gene names differ.", call. = FALSE)
  if (!identical(cells1, barcodes1)) stop("Input 1 HDF5 cells differ from barcodes.")
  if (!identical(cells2, barcodes2)) stop("Input 2 HDF5 cells differ from barcodes.")
  if (length(genes1) != nrow(m1) || nrow(m1) != nrow(m2)) {
    stop("Gene dimensions are inconsistent.", call. = FALSE)
  }

  tree1 <- rhdf5::h5dump(h1, load = FALSE)
  tree2 <- rhdf5::h5dump(h2, load = FALSE)
  if (!setequal(names(tree1$genes), genes1)) stop("Input 1 HDF5 is incomplete.")
  if (!setequal(names(tree2$genes), genes2)) stop("Input 2 HDF5 is incomplete.")

  tmp <- paste0(output_targets, ".tmp_", Sys.getpid())
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE, force = TRUE)
  if (dir.exists(output_targets) && !isTRUE(overwrite)) {
    stop("Output exists: ", output_targets, call. = FALSE)
  }
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  hout <- file.path(tmp, output_h5_name)
  if (!rhdf5::h5createFile(hout)) stop("Could not create ", hout)
  rhdf5::h5createGroup(hout, "genes")
  rhdf5::h5createGroup(hout, "genes_meta")
  rhdf5::h5createGroup(hout, "meta")

  merged_cells <- c(cells1, cells2)
  rhdf5::h5write(merged_cells, hout, "/meta/cell_names")
  rhdf5::h5write(genes1, hout, "/meta/gene_names")

  for (i in seq_along(genes1)) {
    gene <- genes1[[i]]
    if (verbose_every > 0L &&
        (i == 1L || i %% verbose_every == 0L || i == length(genes1))) {
      message("merge_targets_dirs: ", i, " / ", length(genes1), " genes")
    }

    p <- paste0("/genes/", gene)
    q1 <- rhdf5::h5read(h1, p)
    q2 <- rhdf5::h5read(h2, p)
    .targets_write_dataset(
      hout, p, c(q1, q2), compression_level, chunk_target
    )

    mp <- paste0("/genes_meta/", gene)

    meta_names1 <- names(tree1$genes_meta)
    meta_names2 <- names(tree2$genes_meta)

    has_meta1 <- !is.null(meta_names1) && gene %in% meta_names1
    has_meta2 <- !is.null(meta_names2) && gene %in% meta_names2

    if (xor(has_meta1, has_meta2)) {
      stop(
        "/genes_meta/", gene,
        " exists in only one input.",
        call. = FALSE
      )
    }

    if (has_meta1 && has_meta2) {
      meta1 <- rhdf5::h5read(h1, mp)
      meta2 <- rhdf5::h5read(h2, mp)

      expected_children <- c("counts_reads", "lens_reads_plus1")

      if (!identical(sort(names(meta1)), sort(expected_children)) ||
          !identical(sort(names(meta2)), sort(expected_children))) {
        stop(
          "Unexpected metadata fields for gene ", gene,
          ". Input 1: ", paste(names(meta1), collapse = ", "),
          "; input 2: ", paste(names(meta2), collapse = ", "),
          call. = FALSE
        )
      }

      for (field in expected_children) {
        if (length(meta1[[field]]) != length(cells1) ||
            length(meta2[[field]]) != length(cells2)) {
          stop(
            "Invalid /genes_meta/", gene, "/", field,
            " length. Expected ", length(cells1), " and ", length(cells2),
            ", got ", length(meta1[[field]]), " and ", length(meta2[[field]]),
            call. = FALSE
          )
        }
      }

      rhdf5::h5createGroup(hout, mp)

      .targets_write_dataset(
        hout,
        paste0(mp, "/counts_reads"),
        c(meta1$counts_reads, meta2$counts_reads),
        compression_level,
        chunk_target
      )

      .targets_write_dataset(
        hout,
        paste0(mp, "/lens_reads_plus1"),
        c(meta1$lens_reads_plus1, meta2$lens_reads_plus1),
        compression_level,
        chunk_target
      )
    }
  }

  merged_matrix <- cbind(
    Matrix::Matrix(m1, sparse = TRUE),
    Matrix::Matrix(m2, sparse = TRUE)
  )
  .targets_write_mtx(merged_matrix, file.path(tmp, "matrix.mtx.gz"))
  .targets_write_lines(merged_cells, file.path(tmp, "barcodes.tsv.gz"))
  .targets_write_lines(features1, file.path(tmp, "features.tsv.gz"))

  b1 <- file.path(targets1, "burstMCMC.txt")
  b2 <- file.path(targets2, "burstMCMC.txt")
  if (file.exists(b1) || file.exists(b2)) {
    z <- character()
    if (file.exists(b1)) z <- c(z, readLines(b1, warn = FALSE))
    if (file.exists(b2)) z <- c(z, readLines(b2, warn = FALSE))
    writeLines(z, file.path(tmp, "burstMCMC.txt"))
  }

  final_tree <- rhdf5::h5dump(hout, load = FALSE)
  stopifnot(
    identical(as.character(rhdf5::h5read(hout, "/meta/gene_names")), genes1),
    identical(as.character(rhdf5::h5read(hout, "/meta/cell_names")), merged_cells),
    setequal(names(final_tree$genes), genes1),
    nrow(merged_matrix) == length(genes1),
    ncol(merged_matrix) == length(merged_cells)
  )

  if (dir.exists(output_targets)) unlink(output_targets, recursive = TRUE, force = TRUE)
  if (!file.rename(tmp, output_targets)) {
    stop("Merge finished, but could not rename ", tmp, " to ", output_targets)
  }

  manifest <- list(
    targets_dir = normalizePath(output_targets, mustWork = FALSE),
    h5 = normalizePath(file.path(output_targets, output_h5_name), mustWork = FALSE),
    matrix = normalizePath(file.path(output_targets, "matrix.mtx.gz"), mustWork = FALSE),
    barcodes = normalizePath(file.path(output_targets, "barcodes.tsv.gz"), mustWork = FALSE),
    features = normalizePath(file.path(output_targets, "features.tsv.gz"), mustWork = FALSE),
    n_genes = length(genes1),
    n_cells = length(merged_cells),
    n_cells_input1 = length(cells1),
    n_cells_input2 = length(cells2)
  )

  message("Merged targets: ", manifest$n_genes, " genes, ", manifest$n_cells, " cells.")
  invisible(manifest)
}
