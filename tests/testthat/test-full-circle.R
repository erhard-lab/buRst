test_that("steady and non-steady simulated conditions can be compared end-to-end", {
  
  skip_on_cran()
  
  # ------------------------------------------------------------------
  # Temporary directory
  # ------------------------------------------------------------------
  
  td <- tempfile("burst-full-circle-")
  dir.create(td)
  
  on.exit(
    unlink(td, recursive = TRUE, force = TRUE),
    add = TRUE
  )
  
  
  # ------------------------------------------------------------------
  # Define a tiny simulation truth
  # ------------------------------------------------------------------
  
  genes <- paste0("gene_", 1:4)
  
  truth <- data.frame(
    gene = genes,
    
    # steady-state reference parameters
    km_0 = c(0.5, 1.0, 1.5, 2.0),
    B_0  = c(3, 5, 7, 9),
    d_0  = c(0.15, 0.20, 0.25, 0.30),
    
    # known changes in log2 scale
    true_lfc_km = c(
      1,     # gene_1: frequency up
      0,     # gene_2: frequency unchanged
      -1,    # gene_3: frequency down
      0.5    # gene_4: frequency up
    ),
    
    true_lfc_B = c(
      0,      # gene_1: size unchanged
      1,      # gene_2: size up
      0,      # gene_3: size unchanged
      -0.5    # gene_4: size down
    ),
    
    stringsAsFactors = FALSE
  )
  
  truth$km_2 <- truth$km_0 * 2^truth$true_lfc_km
  truth$B_2  <- truth$B_0 * 2^truth$true_lfc_B
  
  # degradation remains unchanged
  truth$d_2 <- truth$d_0
  
  
  # ------------------------------------------------------------------
  # Simulate steady-state reference
  # ------------------------------------------------------------------
  
  prefix_0 <- file.path(
    td,
    "sim_ref"
  )
  
  data_0 <- simulateBursting(
    kms = truth$km_0,
    ds = truth$d_0,
    Bs = truth$B_0,
    
    gene.info = truth$gene,
    
    t.end = 2,
    
    condition = "sim_0",
    
    no.cells = 100,
    
    parallel = FALSE,
    ntr.sim = TRUE,
    burstMCMC = TRUE,
    fit.beta = FALSE,
    
    seed = 1234,
    
    output_prefix = prefix_0
  )
  
  
  # ------------------------------------------------------------------
  # Simulate perturbed non-steady-state condition
  # ------------------------------------------------------------------
  
  prefix_2 <- file.path(
    td,
    "sim_perturbed"
  )
  
  data_2 <- simulateBursting(
    kms = truth$km_2,
    ds = truth$d_2,
    Bs = truth$B_2,
    
    gene.info = truth$gene,
    
    t.end = 2,
    tp.label = 2,
    
    steady = FALSE,
    
    # second condition starts from the reference counts
    init = as.matrix(
      data_0$data$count
    ),
    
    condition = "sim_2",
    
    no.cells = 100,
    
    parallel = FALSE,
    ntr.sim = TRUE,
    burstMCMC = TRUE,
    fit.beta = FALSE,
    
    seed = 5678,
    
    output_prefix = prefix_2
  )
  
  
  # ------------------------------------------------------------------
  # Basic simulation checks
  # ------------------------------------------------------------------
  
  expect_equal(
    nrow(data_0$data$count),
    length(genes)
  )
  
  expect_equal(
    nrow(data_2$data$count),
    length(genes)
  )
  
  expect_true(
    all(is.finite(data_0$data$count))
  )
  
  expect_true(
    all(is.finite(data_2$data$count))
  )
  
  
  # ------------------------------------------------------------------
  # Attach known truth to gene.info
  # ------------------------------------------------------------------
  
  add_truth <- function(data, truth, timepoint) {
    
    if (timepoint == "t0") {
      
      truth_sub <- data.frame(
        Gene = truth$gene,
        true_B.t0 = truth$B_0,
        true_km.t0 = truth$km_0,
        true_d.t0 = truth$d_0,
        stringsAsFactors = FALSE
      )
      
    } else if (timepoint == "t2") {
      
      truth_sub <- data.frame(
        Gene = truth$gene,
        true_B.t2 = truth$B_2,
        true_km.t2 = truth$km_2,
        true_d.t2 = truth$d_2,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      stop(
        "timepoint must be 't0' or 't2'"
      )
    }
    
    idx <- match(
      data$gene.info$Gene,
      truth_sub$Gene
    )
    
    expect_false(
      anyNA(idx)
    )
    
    truth_cols <- setdiff(
      names(truth_sub),
      "Gene"
    )
    
    data$gene.info[truth_cols] <-
      truth_sub[idx, truth_cols, drop = FALSE]
    
    data
  }
  
  
  data_0 <- add_truth(
    data_0,
    truth,
    "t0"
  )
  
  data_2 <- add_truth(
    data_2,
    truth,
    "t2"
  )
  
  
  # ------------------------------------------------------------------
  # Merge the two simulated conditions
  # ------------------------------------------------------------------
  
  data <- grandR::merge.grandR(
    data_0,
    data_2
  )
  
  expect_equal(
    length(grandR::Genes(data)),
    length(genes)
  )
  
  
  # ------------------------------------------------------------------
  # Merge discrete likelihood files
  # ------------------------------------------------------------------
  
  targets_0 <- paste0(
    prefix_0,
    ".targets"
  )
  
  targets_2 <- paste0(
    prefix_2,
    ".targets"
  )
  
  merged_targets <- file.path(
    td,
    "sim_merged.targets"
  )
  
  merge_targets_dirs(
    targets1 = targets_0,
    targets2 = targets_2,
    output_targets = merged_targets,
    h5_filename = "4sU.Binom.discrete.h5",
    overwrite = TRUE,
    verbose_every = 100
  )
  
  h5_path <- file.path(
    merged_targets,
    "4sU.Binom.discrete.h5"
  )
  
  expect_true(
    file.exists(h5_path)
  )
  
  
  # ------------------------------------------------------------------
  # Re-fit both conditions and compare them
  # ------------------------------------------------------------------
  
  out <- run_kinetics_batch(
    data = data,
    genes = grandR::Genes(data),
    
    cond1 = "sim_0",
    cond2 = "sim_2",
    
    t.label = 2,
    t.diff = 2,
    
    fit1_steady = TRUE,
    fit2_steady = FALSE,
    
    fit1_ref_cond = NULL,
    fit2_ref_cond = "sim_0",
    
    h5_path = h5_path,
    
    # keep small for an integration test
    ndraw = 500,
    seed = 123,
    
    do_compare = TRUE,
    
    truth_suffix1 = "t0",
    truth_suffix2 = "t2",
    
    require_h5 = TRUE,
    keep_only_unique_symbols = TRUE,
    
    parallel = FALSE,
    
    timeout_sec = 120,
    verbose_every = 50
  )
  
  
  # ------------------------------------------------------------------
  # Check top-level output
  # ------------------------------------------------------------------
  
  expect_type(
    out,
    "list"
  )
  
  expect_true(
    all(
      c(
        "gene_map",
        "genes_ok",
        "fits",
        "comparisons"
      ) %in% names(out)
    )
  )
  
  
  # ------------------------------------------------------------------
  # Check fitted parameters
  # ------------------------------------------------------------------
  
  expect_s3_class(
    out$fits,
    "data.frame"
  )
  
  expect_true(
    all(
      c(
        "gene",
        "fit_cond",
        "km_hat",
        "B_hat",
        "d_hat",
        "id_flag"
      ) %in% names(out$fits)
    )
  )
  
  expect_true(
    all(
      c(
        "sim_0",
        "sim_2"
      ) %in% unique(out$fits$fit_cond)
    )
  )
  
  expect_gt(
    sum(is.finite(out$fits$km_hat)),
    0
  )
  
  expect_gt(
    sum(is.finite(out$fits$B_hat)),
    0
  )
  
  expect_gt(
    sum(is.finite(out$fits$d_hat)),
    0
  )
  
  
  # ------------------------------------------------------------------
  # Check condition comparison
  # ------------------------------------------------------------------
  
  cmp <- out[["comparisons"]]
  
  expect_s3_class(
    cmp,
    "data.frame"
  )
  
  expect_equal(
    nrow(cmp),
    length(genes)
  )
  
  expect_true(
    all(
      c(
        "gene",
        "cond1",
        "cond2",
        "ok",
        "error_cmp",
        "call_change",
        "top_mechanism",
        "median_lfc_km",
        "median_lfc_B"
      ) %in% names(cmp)
    )
  )
  
  expect_true(
    all(cmp$ok)
  )
  
  expect_true(
    all(is.na(cmp$error_cmp))
  )
  
  expect_true(
    all(cmp$cond1 == "sim_0")
  )
  
  expect_true(
    all(cmp$cond2 == "sim_2")
  )
  
  expect_setequal(
    cmp$gene,
    truth$gene
  )
  
  
  # ------------------------------------------------------------------
  # Join inferred changes to known truth
  # ------------------------------------------------------------------
  
  cmp_truth <- merge(
    cmp,
    truth[
      c(
        "gene",
        "true_lfc_km",
        "true_lfc_B"
      )
    ],
    by = "gene",
    sort = FALSE
  )
  
  
  # ------------------------------------------------------------------
  # Scientific sanity checks:
  # inferred directions should agree with strong simulated changes
  # ------------------------------------------------------------------
  
  # gene_1: frequency increased
  expect_gt(
    cmp_truth$median_lfc_km[
      cmp_truth$gene == "gene_1"
    ],
    0
  )
  
  # gene_2: burst size increased
  expect_gt(
    cmp_truth$median_lfc_B[
      cmp_truth$gene == "gene_2"
    ],
    0
  )
  
  # gene_3: frequency decreased
  expect_lt(
    cmp_truth$median_lfc_km[
      cmp_truth$gene == "gene_3"
    ],
    0
  )
  
  # gene_4: frequency increased
  expect_gt(
    cmp_truth$median_lfc_km[
      cmp_truth$gene == "gene_4"
    ],
    0
  )
  
  # gene_4: burst size decreased
  expect_lt(
    cmp_truth$median_lfc_B[
      cmp_truth$gene == "gene_4"
    ],
    0
  )
  
  
  # ------------------------------------------------------------------
  # Loose quantitative recovery check
  # ------------------------------------------------------------------
  
  mean_abs_error_km <- mean(
    abs(
      cmp_truth$median_lfc_km -
        cmp_truth$true_lfc_km
    )
  )
  
  mean_abs_error_B <- mean(
    abs(
      cmp_truth$median_lfc_B -
        cmp_truth$true_lfc_B
    )
  )
  
  expect_lt(
    mean_abs_error_km,
    1
  )
  
  expect_lt(
    mean_abs_error_B,
    1
  )
})