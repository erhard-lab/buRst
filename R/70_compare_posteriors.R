# km/B posterior comparison
#
# Utilities for comparing burst-frequency (km) and burst-size (B)
# posterior samples between two fitted conditions. External functions are
# namespace-qualified so this file can live directly in a package R/ directory.

# Internal helpers ----

#' Construct paired km/B log2 fold-change draws
#'
#' Keeps `km` and `B` paired within each condition, while sampling the two
#' conditions independently. This preserves within-condition posterior
#' covariance without artificially pairing draw `i` from condition 1 with
#' draw `i` from condition 2.
#'
#' @param samples1 Named list containing posterior vectors `km` and `B` for
#'   condition 1.
#' @param samples2 Named list containing posterior vectors `km` and `B` for
#'   condition 2.
#' @param n Number of comparison draws. If `NULL`, uses the smaller number
#'   of valid posterior draws across the two conditions.
#' @param seed Optional random seed used for independent condition sampling.
#' @param eps Small positive threshold used when filtering non-positive draws.
#'
#' @return A list containing the joint log2FC posterior, cleaned draws,
#'   sampled draws, and the sampled row indices.
#' @keywords internal
#' @noRd
.make_lfc_kmB_independent <- function(
    samples1,
    samples2,
    n = NULL,
    seed = NULL,
    eps = 1e-12
) {

  required <- c("km", "B")

  if (!all(required %in% names(samples1))) {
    stop(
      "samples1 must contain km and B.",
      call. = FALSE
    )
  }

  if (!all(required %in% names(samples2))) {
    stop(
      "samples2 must contain km and B.",
      call. = FALSE
    )
  }

  # Keep km and B together within each condition.
  nraw1 <- min(
    length(samples1$km),
    length(samples1$B)
  )

  nraw2 <- min(
    length(samples2$km),
    length(samples2$B)
  )

  if (nraw1 < 10 || nraw2 < 10) {
    stop(
      "Too few posterior draws before cleaning.",
      call. = FALSE
    )
  }

  draws1 <- data.frame(
    km = as.numeric(
      samples1$km[seq_len(nraw1)]
    ),
    B = as.numeric(
      samples1$B[seq_len(nraw1)]
    )
  )

  draws2 <- data.frame(
    km = as.numeric(
      samples2$km[seq_len(nraw2)]
    ),
    B = as.numeric(
      samples2$B[seq_len(nraw2)]
    )
  )

  # Remove complete rows rather than filtering km and B
  # separately. This preserves within-condition covariance.
  valid1 <- with(
    draws1,
    is.finite(km) &
      is.finite(B) &
      km > eps &
      B > eps
  )

  valid2 <- with(
    draws2,
    is.finite(km) &
      is.finite(B) &
      km > eps &
      B > eps
  )

  draws1 <- draws1[
    valid1,
    ,
    drop = FALSE
  ]

  draws2 <- draws2[
    valid2,
    ,
    drop = FALSE
  ]

  n1 <- nrow(draws1)
  n2 <- nrow(draws2)

  if (n1 < 10 || n2 < 10) {
    stop(
      paste0(
        "Too few valid posterior draws after cleaning. ",
        "Condition 1: ", n1,
        "; condition 2: ", n2, "."
      ),
      call. = FALSE
    )
  }

  if (is.null(n)) {
    n <- min(n1, n2)
  }

  if (
    length(n) != 1 ||
    !is.numeric(n) ||
    !is.finite(n) ||
    n < 10
  ) {
    stop(
      "n must be one finite number of at least 10.",
      call. = FALSE
    )
  }

  n <- as.integer(n)

  if (!is.null(seed)) {
    set.seed(
      as.integer(seed)
    )
  }

  # Crucial correction:
  # sample condition 1 and condition 2 independently.
  idx1 <- sample.int(
    n = n1,
    size = n,
    replace = n > n1
  )

  idx2 <- sample.int(
    n = n2,
    size = n,
    replace = n > n2
  )

  paired1 <- draws1[
    idx1,
    ,
    drop = FALSE
  ]

  paired2 <- draws2[
    idx2,
    ,
    drop = FALSE
  ]

  posterior <- data.frame(
    lfc_km = log2(
      paired2$km /
        paired1$km
    ),
    lfc_B = log2(
      paired2$B /
        paired1$B
    )
  )

  posterior <- posterior[
    is.finite(posterior$lfc_km) &
      is.finite(posterior$lfc_B),
    ,
    drop = FALSE
  ]

  if (nrow(posterior) < 10) {
    stop(
      "Too few finite log2FC draws.",
      call. = FALSE
    )
  }

  list(
    posterior = posterior,
    draws1 = draws1,
    draws2 = draws2,
    paired1 = paired1,
    paired2 = paired2,
    idx1 = idx1,
    idx2 = idx2
  )
}

#' Classify km/B log2 fold-change draws into star regions
#'
#' Assigns each joint `(log2FC(km), log2FC(B))` draw to a practical-change
#' region: no meaningful change, frequency, size, or both. Directional
#' labels are also returned.
#'
#' @param lfc_km Numeric vector of log2 fold changes in burst frequency `km`.
#' @param lfc_B Numeric vector of log2 fold changes in burst size `B`.
#' @param radius_cutoff Radius around the origin considered practically
#'   unchanged.
#' @param frequency_angle Upper angle, in degrees, for frequency-dominated
#'   changes.
#' @param size_angle Lower angle, in degrees, for size-dominated changes.
#'
#' @return A list with `radius`, `angle`, `mechanism`, and
#'   `directional_region` for each draw.
#' @keywords internal
#' @noRd
.classify_star_kmB <- function(
    lfc_km,
    lfc_B,
    radius_cutoff = 0.3,
    frequency_angle = 30,
    size_angle = 60
) {

  radius <- sqrt(
    lfc_km^2 +
      lfc_B^2
  )

  angle <- atan2(
    abs(lfc_B),
    abs(lfc_km)
  ) * 180 / pi

  mechanism <- dplyr::case_when(
    radius <= radius_cutoff ~
      "No meaningful change",

    angle < frequency_angle ~
      "Frequency",

    angle > size_angle ~
      "Size",

    TRUE ~
      "Both"
  )

  directional_region <- dplyr::case_when(
    mechanism == "No meaningful change" ~
      "No meaningful change",

    mechanism == "Frequency" &
      lfc_km >= 0 ~
      "Frequency up",

    mechanism == "Frequency" &
      lfc_km < 0 ~
      "Frequency down",

    mechanism == "Size" &
      lfc_B >= 0 ~
      "Size up",

    mechanism == "Size" &
      lfc_B < 0 ~
      "Size down",

    mechanism == "Both" &
      lfc_km >= 0 &
      lfc_B >= 0 ~
      "Both up",

    mechanism == "Both" &
      lfc_km < 0 &
      lfc_B < 0 ~
      "Both down",

    mechanism == "Both" ~
      "Both opposing",

    TRUE ~
      "Unclassified"
  )

  list(
    radius = radius,
    angle = angle,
    mechanism = mechanism,
    directional_region = directional_region
  )
}

# Posterior comparison ----

#' Compare two km/B posterior distributions with a gated star classifier
#'
#' Compares posterior samples from two fitted conditions in joint
#' `log2FC(km)` / `log2FC(B)` space. The procedure combines three ideas:
#'
#' 1. a Mahalanobis gate testing whether the posterior center is separated
#'    from zero in the joint two-dimensional posterior;
#' 2. a practical-effect gate requiring sufficient posterior mass outside a
#'    user-defined radius around zero; and
#' 3. a star-region classifier that attributes supported changes to burst
#'    frequency, burst size, both, or an opposing trade-off.
#'
#' Posterior draws from the two conditions are sampled independently, while
#' `km` and `B` remain paired within each condition.
#'
#' @param res1 First fitted result. Posterior samples must be stored in
#'   `attr(res1, "samples")` and contain numeric vectors `km` and `B`.
#' @param res2 Second fitted result with the same posterior-sample structure
#'   as `res1`.
#' @param gene Optional label used as the plot title.
#' @param radius_cutoff Practical-change radius in joint log2FC space.
#'   Draws with radius less than or equal to this value are classified as
#'   no meaningful change. Default is `0.3`.
#' @param evidence_level Credible-ellipse level used to define the
#'   chi-squared Mahalanobis cutoff. Default is `0.95`.
#' @param min_outside_mass Minimum posterior probability required outside
#'   `radius_cutoff` to support a meaningful change. Default is `0.90`.
#' @param mechanism_prob_cutoff Minimum conditional posterior mass among
#'   changed draws required to support one mechanism. Default is `0.60`.
#' @param frequency_angle Angle in degrees below which a changed draw is
#'   classified as frequency dominated. Default is `30`.
#' @param size_angle Angle in degrees above which a changed draw is classified
#'   as size dominated. Default is `60`.
#' @param dominance_cutoff For opposing `km`/`B` changes, minimum ratio of
#'   absolute log2FC magnitudes required to call one component dominant.
#'   Default is `1.20`.
#' @param dominance_prob_cutoff Minimum posterior agreement among opposing
#'   draws required to support the dominant component. Default is `0.70`.
#' @param n Optional number of posterior comparison draws. If `NULL`, uses
#'   the smaller number of valid draws across conditions.
#' @param seed Optional random seed for reproducible independent sampling.
#' @param eps Small positive threshold used to exclude non-positive `km` and
#'   `B` draws before taking logarithms.
#' @param ci_level Credible-interval level used for the marginal log2FC
#'   summaries of `km` and `B`. Default is `0.95`.
#' @param lfc_ci_width_veto Optional maximum credible-interval width for
#'   mechanistic attribution. The default `Inf` disables this veto.
#' @param plot Logical; if `TRUE`, construct the km/B star plot.
#' @param limits Optional numeric vector of length two giving symmetric plot
#'   limits. If `NULL`, limits are chosen from the posterior draws.
#' @param max_plot_points Maximum number of posterior draws displayed as
#'   points. This affects plotting only.
#' @param grid_n Number of grid values per axis used for the plot background.
#'
#' @return A named list containing:
#' \describe{
#'   \item{summary}{One-row data frame with the main comparison results.}
#'   \item{call_change}{Logical indicating whether both the joint-evidence
#'     and practical-magnitude gates are passed.}
#'   \item{interpretation}{Concise overall interpretation.}
#'   \item{interpretation_detail}{Diagnostic text describing the gate and
#'     mechanism probabilities.}
#'   \item{attrib_ok}{Logical indicating whether mechanistic attribution
#'     passes the optional credible-interval-width veto.}
#'   \item{region_mass}{Unconditional posterior mass in each star region.}
#'   \item{mechanism_mass}{Mechanism probabilities conditional on meaningful
#'     change.}
#'   \item{directional_mass}{Conditional masses for directional regions.}
#'   \item{km_summary,B_summary}{Credible-interval summaries of the two
#'     log2 fold changes.}
#'   \item{posterior}{Joint posterior draw table used for classification.}
#'   \item{plot}{A `ggplot` object when `plot = TRUE`, otherwise `NULL`.}
#' }
#'
#' @importFrom rlang .data
#' @export
compare_posteriors_gate_star_kmB <- function(
    res1,
    res2,
    gene = NULL,

    radius_cutoff = 0.3,

    # Zero is outside this posterior credible ellipse
    evidence_level = 0.95,

    # Required probability of exceeding the practical radius
    min_outside_mass = 0.90,

    # Required conditional posterior mass for one mechanism
    mechanism_prob_cutoff = 0.60,

    frequency_angle = 30,
    size_angle = 60,

    # For draws in the opposing B-km ridge, one component must be
    # at least this many times larger to be called dominant.
    dominance_cutoff = 1.20,

    # Required posterior agreement among opposing draws.
    dominance_prob_cutoff = 0.70,

    n = NULL,
    seed = NULL,
    eps = 1e-12,

    # Optional attribution reliability threshold
    ci_level = 0.95,
    lfc_ci_width_veto = Inf,

    plot = TRUE,
    limits = NULL,
    max_plot_points = 3000,
    grid_n = 400
) {

  # A. Correctly construct the joint LFC posterior
  samples1 <- attr(res1, "samples")
  samples2 <- attr(res2, "samples")
  comparison <- .make_lfc_kmB_independent(
    samples1 = samples1,
    samples2 = samples2,
    n = n,
    seed = seed,
    eps = eps
  )

  posterior <- comparison$posterior

  # B. Joint Mahalanobis gate
  #
  # This asks whether zero lies outside the posterior
  # credible ellipse.

  posterior_matrix <- as.matrix(
    posterior[
      ,
      c("lfc_km", "lfc_B"),
      drop = FALSE
    ]
  )

  posterior_center <- colMeans(
    posterior_matrix
  )

  posterior_cov <- stats::cov(
    posterior_matrix
  )

  # Small ridge term for numerical stability.
  covariance_regularization <- max(
    mean(diag(posterior_cov)) * 1e-8,
    1e-10
  )

  posterior_cov_regularized <-
    posterior_cov +
    diag(
      covariance_regularization,
      2
    )

  mahalanobis_d2 <- as.numeric(
    stats::mahalanobis(
      x = matrix(
        posterior_center,
        nrow = 1
      ),
      center = c(0, 0),
      cov = posterior_cov_regularized
    )
  )

  mahalanobis_cutoff <- stats::qchisq(
    evidence_level,
    df = 2
  )

  joint_evidence <-
    mahalanobis_d2 >
    mahalanobis_cutoff

  # C. Practical-effect circle and star regions

  star <- .classify_star_kmB(
    lfc_km = posterior$lfc_km,
    lfc_B = posterior$lfc_B,
    radius_cutoff = radius_cutoff,
    frequency_angle = frequency_angle,
    size_angle = size_angle
  )

  posterior <- posterior |>
    dplyr::mutate(
      radius = as.numeric(star$radius),
      angle = as.numeric(star$angle),
      outside = .data$radius > radius_cutoff,
      mechanism = as.character(star$mechanism),
      directional_region = as.character(star$directional_region),

      # Opposing means km and B move in different directions.
      opposing =
        .data$outside &
        (.data$lfc_km * .data$lfc_B < 0),

      # Among opposing draws, ask which absolute LFC is larger.
      opposing_dominance = dplyr::case_when(
        !.data$opposing ~ NA_character_,
        abs(.data$lfc_km) >=
          dominance_cutoff * abs(.data$lfc_B) ~
          "Frequency dominated",
        abs(.data$lfc_B) >=
          dominance_cutoff * abs(.data$lfc_km) ~
          "Size dominated",
        TRUE ~ "Balanced trade-off"
      )
    )

  outside_mass <- mean(
    posterior$outside
  )

  magnitude_supported <-
    outside_mass >=
    min_outside_mass

  # D. Region probabilities over all posterior draws

  region_mass <- posterior |>
    dplyr::count(
      .data$mechanism,
      name = "n"
    ) |>
    dplyr::mutate(
      posterior_mass = .data$n /
        nrow(posterior)
    ) |>
    tidyr::complete(
      mechanism = c(
        "No meaningful change",
        "Frequency",
        "Size",
        "Both"
      ),
      fill = list(
        n = 0,
        posterior_mass = 0
      )
    )

  # E. Mechanism probabilities conditional on meaningful change

  changed_draws <- posterior[
    posterior$outside %in% TRUE,
    ,
    drop = FALSE
  ]

  if (nrow(changed_draws) > 0) {

    mechanism_mass <- changed_draws |>
      dplyr::count(
        .data$mechanism,
        name = "n"
      ) |>
      tidyr::complete(
        mechanism = c(
          "Frequency",
          "Size",
          "Both"
        ),
        fill = list(n = 0)
      ) |>
      dplyr::mutate(
        conditional_mass =
          .data$n / sum(.data$n)
      )

    directional_mass <- changed_draws |>
      dplyr::count(
        .data$directional_region,
        name = "n"
      ) |>
      dplyr::mutate(
        conditional_mass =
          .data$n / sum(.data$n)
      )

    opposing_draws <- changed_draws[
      changed_draws$opposing %in% TRUE,
      ,
      drop = FALSE
    ]

    if (nrow(opposing_draws) > 0L) {
      opposing_dominance_mass <- opposing_draws |>
        dplyr::count(
          .data$opposing_dominance,
          name = "n"
        ) |>
        tidyr::complete(
          opposing_dominance = c(
            "Frequency dominated",
            "Size dominated",
            "Balanced trade-off"
          ),
          fill = list(n = 0L)
        ) |>
        dplyr::mutate(
          conditional_mass =
            .data$n / sum(.data$n)
        )
    } else {
      opposing_dominance_mass <- tibble::tibble(
        opposing_dominance = c(
          "Frequency dominated",
          "Size dominated",
          "Balanced trade-off"
        ),
        n = 0L,
        conditional_mass = 0
      )
    }

  } else {

    mechanism_mass <- tibble::tibble(
      mechanism = c(
        "Frequency",
        "Size",
        "Both"
      ),
      n = 0L,
      conditional_mass = 0
    )

    directional_mass <- tibble::tibble(
      directional_region = character(0),
      n = integer(0),
      conditional_mass = numeric(0)
    )

    opposing_dominance_mass <- tibble::tibble(
      opposing_dominance = c(
        "Frequency dominated",
        "Size dominated",
        "Balanced trade-off"
      ),
      n = 0L,
      conditional_mass = 0
    )
  }

  top_mechanism <- mechanism_mass |>
    dplyr::arrange(
      dplyr::desc(.data$conditional_mass)
    ) |>
    dplyr::slice(1)

  top_directional_region <- if (
    nrow(directional_mass) > 0
  ) {
    directional_mass |>
      dplyr::arrange(
        dplyr::desc(.data$conditional_mass)
      ) |>
      dplyr::slice(1)
  } else {
    tibble::tibble(
      directional_region = NA_character_,
      n = 0L,
      conditional_mass = NA_real_
    )
  }

  mechanism_supported <-
    top_mechanism$conditional_mass >=
    mechanism_prob_cutoff

  opposing_mass <- if (nrow(changed_draws) > 0L) {
    mean(changed_draws$opposing %in% TRUE)
  } else {
    0
  }

  opposing_supported <-
    identical(
      as.character(
        top_directional_region$directional_region[[1]]
      ),
      "Both opposing"
    ) &&
    opposing_mass >= mechanism_prob_cutoff

  top_opposing_dominance <- opposing_dominance_mass |>
    dplyr::arrange(
      dplyr::desc(.data$conditional_mass)
    ) |>
    dplyr::slice(1L)

  dominance_supported <-
    opposing_supported &&
    top_opposing_dominance$conditional_mass[[1]] >=
    dominance_prob_cutoff

  # F. Credible intervals and attribution reliability

  alpha <- (
    1 -
      ci_level
  ) / 2

  summarize_lfc <- function(x) {

    qs <- stats::quantile(
      x,
      probs = c(
        alpha,
        0.5,
        1 - alpha
      ),
      names = FALSE,
      na.rm = TRUE
    )

    tibble::tibble(
      lower = qs[1],
      median = qs[2],
      upper = qs[3],
      width = qs[3] - qs[1]
    )
  }

  km_summary <- summarize_lfc(
    posterior$lfc_km
  )

  B_summary <- summarize_lfc(
    posterior$lfc_B
  )

  attribution_unreliable <- any(
    c(
      km_summary$width,
      B_summary$width
    ) > lfc_ci_width_veto
  )

  # G. Final interpretation

  interpretation <- dplyr::case_when(
    !joint_evidence ~
      "No joint evidence",

    !magnitude_supported ~
      "No supported meaningful change",

    attribution_unreliable ~
      "Change supported; attribution unreliable",

    opposing_supported &&
      dominance_supported &&
      top_opposing_dominance$opposing_dominance[[1]] ==
      "Frequency dominated" ~
      "Frequency supported; opposing B shift",

    opposing_supported &&
      dominance_supported &&
      top_opposing_dominance$opposing_dominance[[1]] ==
      "Size dominated" ~
      "Size supported; opposing km shift",

    opposing_supported ~
      "Change supported; balanced opposing trade-off",

    !mechanism_supported ~
      "Change supported; mechanism unresolved",

    top_mechanism$mechanism ==
      "Frequency" ~
      "Frequency supported",

    top_mechanism$mechanism ==
      "Size" ~
      "Size supported",

    top_mechanism$mechanism ==
      "Both" ~
      "Both concordant supported",

    TRUE ~
      "Unresolved"
  )

  interpretation_detail <- paste0(
    "D\u00B2 = ",
    sprintf("%.2f", mahalanobis_d2),
    " versus cutoff ",
    sprintf("%.2f", mahalanobis_cutoff),
    "; P(radius > ",
    radius_cutoff,
    ") = ",
    sprintf("%.3f", outside_mass),
    "; top conditional mechanism mass = ",
    sprintf(
      "%.3f",
      top_mechanism$conditional_mass
    ),
    " for ",
    top_mechanism$mechanism,
    "; opposing mass = ",
    sprintf("%.3f", opposing_mass),
    "; top opposing dominance = ",
    top_opposing_dominance$opposing_dominance[[1]],
    " (",
    sprintf(
      "%.3f",
      top_opposing_dominance$conditional_mass[[1]]
    ),
    ")."
  )

  # H. Optional star plot

  plt <- NULL

  if (isTRUE(plot)) {

    if (is.null(limits)) {

      lim <- max(
        stats::quantile(
          abs(c(
            posterior$lfc_km,
            posterior$lfc_B
          )),
          probs = 0.995,
          na.rm = TRUE
        ),
        radius_cutoff,
        na.rm = TRUE
      )

      limits <- c(
        -1.08 * lim,
        1.08 * lim
      )
    }

    if (
      length(limits) != 2 ||
      !all(is.finite(limits)) ||
      limits[1] >= limits[2]
    ) {
      stop(
        "limits must contain two increasing finite values.",
        call. = FALSE
      )
    }

    grid_values <- seq(
      limits[1],
      limits[2],
      length.out = grid_n
    )

    background <- tidyr::expand_grid(
      lfc_km = grid_values,
      lfc_B = grid_values
    )

    background_star <- .classify_star_kmB(
      lfc_km = background$lfc_km,
      lfc_B = background$lfc_B,
      radius_cutoff = radius_cutoff,
      frequency_angle = frequency_angle,
      size_angle = size_angle
    )

    background$region <-
      background_star$mechanism

    circle_angles <- seq(
      0,
      2 * pi,
      length.out = 500
    )

    circle_df <- tibble::tibble(
      x = radius_cutoff *
        cos(circle_angles),
      y = radius_cutoff *
        sin(circle_angles)
    )

    boundary_angles <- c(
      frequency_angle,
      size_angle,
      180 - size_angle,
      180 - frequency_angle,
      180 + frequency_angle,
      180 + size_angle,
      360 - size_angle,
      360 - frequency_angle
    )

    plot_radius <- max(
      abs(limits)
    ) * sqrt(2)

    wedge_df <- tibble::tibble(
      angle = boundary_angles *
        pi / 180
    ) |>
      dplyr::mutate(
        x = plot_radius * cos(.data$angle),
        y = plot_radius * sin(.data$angle)
      )

    set.seed(
      if (is.null(seed)) 1L else
        as.integer(seed) + 1L
    )

    shown_draws <- posterior |>
      dplyr::slice_sample(
        n = min(
          max_plot_points,
          nrow(posterior)
        )
      )

    region_colors <- c(
      "No meaningful change" = "#ECECEC",
      "Frequency" = "#F4A3A3",
      "Size" = "#8FD19E",
      "Both" = "#F7C59F"
    )

    mass_caption <- mechanism_mass |>
      dplyr::mutate(
        text = paste0(
          .data$mechanism,
          ": ",
          sprintf(
            "%.1f%%",
            100 * .data$conditional_mass
          )
        )
      ) |>
      dplyr::pull("text") |>
      paste(
        collapse = "   "
      )

    plt <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = background,
        ggplot2::aes(
          x = .data$lfc_km,
          y = .data$lfc_B,
          fill = .data$region
        ),
        alpha = 0.48
      ) +
      ggplot2::geom_segment(
        data = wedge_df,
        ggplot2::aes(
          x = 0,
          y = 0,
          xend = .data$x,
          yend = .data$y
        ),
        inherit.aes = FALSE,
        linewidth = 0.5,
        linetype = 2,
        color = "#555555"
      ) +
      ggplot2::geom_path(
        data = circle_df,
        ggplot2::aes(
          x = .data$x,
          y = .data$y
        ),
        inherit.aes = FALSE,
        linewidth = 0.75,
        linetype = 2,
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.3,
        color = "#555555"
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linewidth = 0.3,
        color = "#555555"
      ) +
      ggplot2::stat_density_2d(
        data = posterior,
        ggplot2::aes(
          x = .data$lfc_km,
          y = .data$lfc_B
        ),
        bins = 7,
        linewidth = 0.55,
        color = "#222222"
      ) +
      ggplot2::geom_point(
        data = shown_draws,
        ggplot2::aes(
          x = .data$lfc_km,
          y = .data$lfc_B
        ),
        size = 0.45,
        alpha = 0.08,
        color = "#222222"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = km_summary$median,
          y = B_summary$median
        ),
        inherit.aes = FALSE,
        shape = 4,
        size = 4,
        stroke = 1.3,
        color = "black"
      ) +
      ggplot2::scale_fill_manual(
        values = region_colors,
        breaks = c(
          "No meaningful change",
          "Frequency",
          "Size",
          "Both"
        ),
        name = "Joint region"
      ) +
      ggplot2::coord_fixed(
        xlim = limits,
        ylim = limits,
        expand = FALSE
      ) +
      ggplot2::labs(
        title = if (
          is.null(gene)
        ) {
          "Joint kinetic classification"
        } else {
          gene
        },
        subtitle = interpretation,
        x = expression(
          "Burst frequency " *
            log[2] * "FC (" * k[m] * ")"
        ),
        y = expression(
          "Burst size " *
            log[2] * "FC (" * B * ")"
        ),
        caption = paste0(
          interpretation_detail,
          "\nConditional mechanism mass: ",
          mass_caption,
          "\nCondition posteriors were sampled independently; ",
          "km and B remained paired within each condition."
        )
      ) +
      ggplot2::theme_classic(
        base_size = 12
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 15
        ),
        plot.subtitle = ggplot2::element_text(
          face = "bold"
        ),
        plot.caption = ggplot2::element_text(
          hjust = 0
        ),
        legend.position = "bottom"
      )
  }

  # I. Return everything needed downstream
  get_mass <- function(
    tab,
    key,
    key_col,
    value_col
  ) {
    if (
      is.null(tab) ||
      !all(c(key_col, value_col) %in% names(tab))
    ) {
      return(NA_real_)
    }

    hit <- tab[
      tab[[key_col]] == key,
      value_col,
      drop = TRUE
    ]

    if (length(hit) == 0) {
      return(0)
    }

    as.numeric(hit[[1]])
  }

  summary_row <- data.frame(
    call_change = joint_evidence &&
      magnitude_supported,

    interpretation = interpretation,
    interpretation_detail = interpretation_detail,

    attrib_ok = !attribution_unreliable,

    # Gate
    mahalanobis_d2 = mahalanobis_d2,
    mahalanobis_cutoff = mahalanobis_cutoff,
    joint_evidence = joint_evidence,
    outside_mass = outside_mass,
    magnitude_supported = magnitude_supported,

    # Mechanism
    mechanism_supported = mechanism_supported,
    top_mechanism = as.character(
      top_mechanism$mechanism[[1]]
    ),
    top_mechanism_mass = as.numeric(
      top_mechanism$conditional_mass[[1]]
    ),

    top_directional_region = if (
      nrow(top_directional_region) > 0
    ) {
      as.character(
        top_directional_region$directional_region[[1]]
      )
    } else {
      NA_character_
    },

    opposing_mass = opposing_mass,
    opposing_supported = opposing_supported,
    dominance_supported = dominance_supported,
    top_opposing_dominance = as.character(
      top_opposing_dominance$opposing_dominance[[1]]
    ),
    top_opposing_dominance_mass = as.numeric(
      top_opposing_dominance$conditional_mass[[1]]
    ),
    conditional_opposing_frequency = get_mass(
      opposing_dominance_mass,
      "Frequency dominated",
      "opposing_dominance",
      "conditional_mass"
    ),
    conditional_opposing_size = get_mass(
      opposing_dominance_mass,
      "Size dominated",
      "opposing_dominance",
      "conditional_mass"
    ),
    conditional_opposing_balanced = get_mass(
      opposing_dominance_mass,
      "Balanced trade-off",
      "opposing_dominance",
      "conditional_mass"
    ),

    # Unconditional posterior mass
    mass_no_change = get_mass(
      region_mass,
      "No meaningful change",
      "mechanism",
      "posterior_mass"
    ),

    mass_frequency = get_mass(
      region_mass,
      "Frequency",
      "mechanism",
      "posterior_mass"
    ),

    mass_size = get_mass(
      region_mass,
      "Size",
      "mechanism",
      "posterior_mass"
    ),

    mass_both = get_mass(
      region_mass,
      "Both",
      "mechanism",
      "posterior_mass"
    ),

    # Conditional mass among changed draws
    conditional_frequency = get_mass(
      mechanism_mass,
      "Frequency",
      "mechanism",
      "conditional_mass"
    ),

    conditional_size = get_mass(
      mechanism_mass,
      "Size",
      "mechanism",
      "conditional_mass"
    ),

    conditional_both = get_mass(
      mechanism_mass,
      "Both",
      "mechanism",
      "conditional_mass"
    ),

    # LFC summaries
    median_lfc_km = km_summary$median[[1]],
    lfc_km_lo = km_summary$lower[[1]],
    lfc_km_hi = km_summary$upper[[1]],
    lfc_km_width = km_summary$width[[1]],

    median_lfc_B = B_summary$median[[1]],
    lfc_B_lo = B_summary$lower[[1]],
    lfc_B_hi = B_summary$upper[[1]],
    lfc_B_width = B_summary$width[[1]],

    # Pairing and draw counts
    pairing = "independent condition indices",
    n_posterior = nrow(posterior),
    valid_draws_cond1 = nrow(comparison$draws1),
    valid_draws_cond2 = nrow(comparison$draws2),

    # Settings
    radius_cutoff = radius_cutoff,
    evidence_level = evidence_level,
    min_outside_mass = min_outside_mass,
    mechanism_prob_cutoff = mechanism_prob_cutoff,
    dominance_cutoff = dominance_cutoff,
    dominance_prob_cutoff = dominance_prob_cutoff,
    frequency_angle = frequency_angle,
    size_angle = size_angle,

    stringsAsFactors = FALSE
  )
  list(
    summary = summary_row,

    call_change =
      joint_evidence &&
      magnitude_supported,

    interpretation = interpretation,
    interpretation_detail =
      interpretation_detail,

    attrib_ok =
      !attribution_unreliable,

    # Gate
    mahalanobis_d2 =
      mahalanobis_d2,

    mahalanobis_cutoff =
      mahalanobis_cutoff,

    joint_evidence =
      joint_evidence,

    outside_mass =
      outside_mass,

    magnitude_supported =
      magnitude_supported,

    # Mechanism
    mechanism_supported =
      mechanism_supported,

    top_mechanism =
      top_mechanism$mechanism,

    top_mechanism_mass =
      top_mechanism$conditional_mass,

    top_directional_region =
      top_directional_region$directional_region,

    opposing_mass = opposing_mass,
    opposing_supported = opposing_supported,
    dominance_supported = dominance_supported,
    top_opposing_dominance =
      top_opposing_dominance$opposing_dominance,
    top_opposing_dominance_mass =
      top_opposing_dominance$conditional_mass,
    opposing_dominance_mass =
      opposing_dominance_mass,

    region_mass =
      region_mass,

    mechanism_mass =
      mechanism_mass,

    directional_mass =
      directional_mass,

    # LFC summaries
    km_summary =
      km_summary,

    B_summary =
      B_summary,

    # Posterior draws
    posterior =
      posterior,

    # Pairing diagnostics
    pairing =
      "independent condition indices",

    idx1 =
      comparison$idx1,

    idx2 =
      comparison$idx2,

    valid_draws_cond1 =
      nrow(comparison$draws1),

    valid_draws_cond2 =
      nrow(comparison$draws2),

    # Settings
    settings = list(
      radius_cutoff =
        radius_cutoff,

      evidence_level =
        evidence_level,

      min_outside_mass =
        min_outside_mass,

      mechanism_prob_cutoff =
        mechanism_prob_cutoff,

      dominance_cutoff =
        dominance_cutoff,

      dominance_prob_cutoff =
        dominance_prob_cutoff,

      frequency_angle =
        frequency_angle,

      size_angle =
        size_angle,

      n =
        nrow(posterior),

      seed =
        seed
    ),

    plot = plt
  )
}
