#' HDF5 I/O for discrete likelihood profiles
#'
#' Stores per-cell discrete likelihood vectors (quantized) in an HDF5 file
#' for fast downstream inference.
#'
#' @keywords internal
#' @noRd
NULL

dequantize_ell_u16_max0 <- function(q_raw, R_clip = 60) {
  q_raw * (R_clip / 65535)
}

quantize_ell_u16_max0 <- function(ell, R_clip = 60) {
  ell0 <- ell - max(ell)
  d <- pmin(pmax(-ell0, 0), R_clip)
  q <- as.integer(round(65535 * d / R_clip))
  q <- pmin(pmax(q, 0L), 65535L)
  q
}

write_all_gene_profiles_from_disc <- function(h5_path, disc, counts,
                                              gene_names = rownames(counts),
                                              cell_names = colnames(counts),
                                              overwrite_file = FALSE,
                                              chunk_target = 1e6,
                                              compression_level = 4) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Need rhdf5.")

  G <- nrow(counts)
  C <- ncol(counts)

  if (length(disc) != G * C) {
    stop(sprintf("disc length (%d) must equal nrow(counts)*ncol(counts) = %d",
                 length(disc), G * C))
  }
  if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(G))
  if (is.null(cell_names)) cell_names <- paste0("cell_", seq_len(C))

  if (overwrite_file && file.exists(h5_path)) file.remove(h5_path)
  if (!file.exists(h5_path)) rhdf5::h5createFile(h5_path)

  try(rhdf5::h5createGroup(h5_path, "genes"), silent = TRUE)
  try(rhdf5::h5createGroup(h5_path, "genes_meta"), silent = TRUE)
  try(rhdf5::h5createGroup(h5_path, "meta"), silent = TRUE)

  try(rhdf5::h5write(cell_names, h5_path, "/meta/cell_names"), silent = TRUE)

  h5_has_path <- function(file, path) {
    ls <- rhdf5::h5ls(file, recursive = TRUE)
    full <- paste0(ls$group, "/", ls$name)
    path %in% full
  }

  for (g_idx in seq_len(G)) {
    gene <- gene_names[g_idx]
    dset <- paste0("/genes/", gene)

    if (h5_has_path(h5_path, dset)) rhdf5::h5delete(h5_path, dset)

    idx <- seq(from = g_idx, to = G * C, by = G)

    counts_g <- as.integer(counts[g_idx, ])
    lens <- ifelse(counts_g == 0L, 0L, counts_g + 1L)

    ok <- vapply(seq_along(idx), function(ci) {
      q <- disc[[idx[ci]]]
      if (lens[ci] == 0L) length(q) == 0L else length(q) == lens[ci]
    }, logical(1))

    if (!all(ok)) {
      bad <- which(!ok)[1]
      stop(sprintf(
        "Length mismatch for gene '%s', cell '%s': counts=%d => expected %d, got %d",
        gene, cell_names[bad], counts_g[bad], lens[bad], length(disc[[idx[bad]]])
      ))
    }

    flat <- unlist(disc[idx], use.names = FALSE)
    if (length(flat) == 0L) flat <- as.integer(0)
    flat <- as.integer(flat)

    ch <- min(length(flat), as.integer(chunk_target))
    if (ch < 1L) ch <- 1L

    rhdf5::h5createDataset(
      file = h5_path,
      dataset = dset,
      dims = length(flat),
      storage.mode = "integer",
      chunk = c(ch),
      level = compression_level
    )
    rhdf5::h5write(flat, h5_path, dset)

    meta_grp <- paste0("/genes_meta/", gene)
    try(rhdf5::h5createGroup(h5_path, meta_grp), silent = TRUE)

    cpath <- paste0(meta_grp, "/counts_reads")
    lpath <- paste0(meta_grp, "/lens_reads_plus1")
    if (h5_has_path(h5_path, cpath)) rhdf5::h5delete(h5_path, cpath)
    if (h5_has_path(h5_path, lpath)) rhdf5::h5delete(h5_path, lpath)

    rhdf5::h5write(as.integer(counts_g), h5_path, cpath)
    rhdf5::h5write(as.integer(lens),     h5_path, lpath)

    if (g_idx %% 200 == 0) message(sprintf("Wrote %d / %d genes", g_idx, G))
  }

  invisible(TRUE)
}

read_gene_profiles <- function(h5_path, gene, counts, R_clip = 60) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Need rhdf5.")

  dataset_path <- paste0("/genes/", gene)
  raw_data <- as.vector(rhdf5::h5read(h5_path, dataset_path))
  diffs <- raw_data * (R_clip / 65535.0)

  lens <- ifelse(counts == 0, 0, counts + 1)
  if (sum(lens) != length(diffs)) stop("Data length does not match counts!")

  ends <- cumsum(lens)
  starts <- ends - lens + 1

  re <- lapply(seq_along(counts), function(i) {
    if (lens[i] == 0) numeric(0) else diffs[starts[i]:ends[i]]
  })
  names(re) <- names(counts)
  re
}

# read_gene_profiles_subset <- function(h5_path, gene, counts_all, subset_names, R_clip = 60) {
#   gp_all <- read_gene_profiles(h5_path, gene = gene, counts = counts_all, R_clip = R_clip)
#   gp_all[subset_names]
# }

load_10x_reference <- function(directory) {
  message(
    "LOADING 10X REFERENCE: ",
    directory,
    " | PID: ",
    Sys.getpid(),
    " | time: ",
    format(Sys.time(), "%H:%M:%S")
  )
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  barcode_path <- file.path(directory, "barcodes.tsv.gz")
  feature_path <- file.path(directory, "features.tsv.gz")
  matrix_path <- file.path(directory, "matrix.mtx.gz")

  if (!file.exists(feature_path)) {
    feature_path <- file.path(directory, "genes.tsv.gz")
  }

  barcodes <- readLines(gzfile(barcode_path))

  features <- read.delim(
    gzfile(feature_path),
    header = FALSE,
    stringsAsFactors = FALSE
  )

  counts <- Matrix::readMM(gzfile(matrix_path))

  if (ncol(counts) != length(barcodes)) {
    stop("Matrix columns do not match barcodes.")
  }

  if (nrow(counts) != nrow(features)) {
    stop("Matrix rows do not match features.")
  }

  rownames(counts) <- features[[1L]]
  colnames(counts) <- barcodes

  list(
    counts = counts,
    features = features,
    barcodes = barcodes
  )
}
read_gene_profiles_subset <- function(
    h5_path,
    gene,
    subset_names,
    reference,
    R_clip = 60
) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required.")
  }

  gene_ids <- sub("\\.[0-9]+$", "", reference$features[[1L]])
  gene_clean <- sub("\\.[0-9]+$", "", gene)

  gene_row <- match(gene_clean, gene_ids)

  if (is.na(gene_row)) {
    stop("Gene not found in features.tsv: ", gene)
  }

  all_barcodes <- reference$barcodes

  # Important: no4sU also contains "4sU"
  has_no4su_field <- vapply(
    strsplit(all_barcodes, ".", fixed = TRUE),
    function(fields) any(tolower(fields) == "no4su"),
    logical(1)
  )

  h5_barcodes <- all_barcodes[!has_no4su_field]

  counts_all <- as.numeric(
    reference$counts[gene_row, , drop = TRUE]
  )
  names(counts_all) <- all_barcodes

  counts_h5 <- counts_all[h5_barcodes]

  fit_indices <- match(subset_names, h5_barcodes)

  if (anyNA(fit_indices)) {
    missing <- unique(subset_names[is.na(fit_indices)])

    stop(
      "Fit cells not found among HDF5 4sU cells: ",
      paste(head(missing, 10L), collapse = ", ")
    )
  }

  dataset_gene <- reference$features[gene_row, 1L]
  dataset_gene <- sub("\\.[0-9]+$", "", dataset_gene)
  dataset_path <- paste0("/genes/", dataset_gene)

  raw_data <- as.vector(
    rhdf5::h5read(h5_path, dataset_path)
  )

  lens <- ifelse(
    counts_h5 == 0L,
    0L,
    counts_h5 + 1L
  )

  if (sum(lens) != length(raw_data)) {
    stop(
      "Data length does not match for gene ", gene, ".\n",
      "Expected from 4sU counts: ", sum(lens), "\n",
      "Observed in HDF5: ", length(raw_data)
    )
  }

  ends <- cumsum(lens)
  starts <- ends - lens + 1L

  profiles <- lapply(fit_indices, function(i) {
    if (lens[i] == 0L) {
      numeric(0)
    } else {
      raw_data[starts[i]:ends[i]] * (R_clip / 65535)
    }
  })

  names(profiles) <- subset_names
  profiles
}
