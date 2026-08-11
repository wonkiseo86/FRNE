############################################################
## FRNE simulation table generator
##
## This script loads saved Monte Carlo output files and writes
## all requested LaTeX tables in a prespecified order.
##
## IMPORTANT:
## 1. Run all simulation scripts before running this file (see README-SUPP).
## 2. Keep the final output files in the directory specified by
##    `result_dirs`.
## 3. The current alternative-DGP file names do not necessarily
##    encode decfac2. Remove or relocate obsolete alternative-DGP
##    files before generating the tables.
##
## Note: Assuming that outputs from other simulation codes are stored
##       in the 'results' subfolder, this script loads these pre-saved 
##       results to reproduce the tables in the Supplementary Material. 
##       If the save directory was modified in the simulation scripts, 
##       root_dir and result_dirs must be updated accordingly. 
##       The authors utilized ChatGPT-5.5 to generate this summary code.
############################################################

rm(list = ls())

############################################################
## 1. User settings
############################################################

root_dir <- "Path, e.g., C:/Users/R_code"
setwd(root_dir)

## Search directories in order. The first matching file is used.
## It is safest to keep only the final results directory here.
result_dirs <- file.path(root_dir, c("results"))

output_tex <- file.path(root_dir, "FRNE_simulation_tables.tex")
output_manifest <- file.path(root_dir, "FRNE_simulation_table_manifest.csv")

## Print the complete LaTeX output to the console as well as to the file.
print_tables_to_console <- TRUE

## Wrap revised material in \commRV{...}.
use_revision_markup <- TRUE

## Table formatting.
table_position <- "p"
array_stretch <- 1
table_font_command <- "\\small"
T_vals <- c(100, 200, 400, 800, 1600)
kappa_vals <- c(0, 1, 2)
rho_list <- c(0, 0.2, 0.4)

## File-name parameters used by the simulation scripts.
Kscheme <- 1
Kval_cut <- 0.75
cuteigen <- 0.5
decfac1 <- 0.2
alternative_decfac2 <- 0.2

## Error-magnitude lists.
report_esig150 <- FALSE
baseline_esignals <- if (report_esig150) c(0, 50, 100, 150) else c(0, 50, 100)
cumulative_esignals <- baseline_esignals
alternative_hs_esignals <- c(0, 50, 100)

## Error magnitudes reported for the alternative-DGP CI tables.
alternative_ci_esignals <- c(0, 50, 100)

## Growth-rate values.
baseline_bdw <- 0.30
sensitivity_bdw <- c(0.25, 0.35)
alternative_bdw <- c(0.25, 0.30, 0.35)

## Switches controlling the output sequence.
include_baseline_main <- TRUE
include_cumulative_main <- TRUE
include_rate_sensitivity <- TRUE
include_alternative_dgp <- TRUE
include_selected_K_tables <- FALSE

## Append the size and power tables for the residual-based test.
include_residual_test_tables <- TRUE
residual_test_dn <- 2
residual_test_bdw <- 0.30
residual_test_esignals <- c(0, 50, 100)
residual_test_rhos <- c(0, 0.2, 0.4)
residual_size_ranintmax <- 1

## Keep this FALSE to reproduce only the tables currently reported in
## the Supplementary Material. If TRUE, additional d_N=3 sensitivity
## tables are appended after the reported tables.
include_unreported_dn3_tables <- FALSE

############################################################
## 2. Table plan
##
## The entries below are listed in exactly the same order as
## the simulation tables currently appear in the Supplementary
## Material. Hence, the generated LaTeX file can be read from
## top to bottom in the same order as the Supplement.
############################################################

make_table_spec <- function(family, metric, dn, bdw, esignals,
                            show_bdw = FALSE, caption = NULL,
                            label = NULL) {
  list(
    family = family,
    metric = metric,
    dn = dn,
    bdw = bdw,
    esignals = esignals,
    show_bdw = show_bdw,
    caption = caption,
    label = label
  )
}

table_plan <- list()

## ---------------------------------------------------------
## A. Baseline simulation results in Section Sec_monte_carlo2
##
## Supplement order:
##   HS, d_N=2
##   CI, d_N=2
##   HS, d_N=3
##   CI, d_N=3
## ---------------------------------------------------------

if (include_baseline_main) {
  table_plan <- c(
    table_plan,
    list(
      make_table_spec(
        family = "baseline", metric = "HS", dn = 2,
        bdw = baseline_bdw, esignals = baseline_esignals,
        label = "tabsim_supp0_dn2"
      ),
      make_table_spec(
        family = "baseline", metric = "CI", dn = 2,
        bdw = baseline_bdw, esignals = baseline_esignals,
        label = "tabsimci_supp0_dn2"
      ),
      make_table_spec(
        family = "baseline", metric = "HS", dn = 3,
        bdw = baseline_bdw, esignals = baseline_esignals,
        label = "tabsim_supp0_dn3"
      ),
      make_table_spec(
        family = "baseline", metric = "CI", dn = 3,
        bdw = baseline_bdw, esignals = baseline_esignals,
        label = "tabsimci_supp0_dn3"
      )
    )
  )
}

## ---------------------------------------------------------
## B. Cumulative-variation criterion in Section sec_addsim1
##
## Supplement order:
##   HS, d_N=2
##   CI, d_N=2
## ---------------------------------------------------------

if (include_cumulative_main) {
  table_plan <- c(
    table_plan,
    list(
      make_table_spec(
        family = "cumulative", metric = "HS", dn = 2,
        bdw = baseline_bdw, esignals = cumulative_esignals,
        label = "tabsim_supp1_dn2"
      ),
      make_table_spec(
        family = "cumulative", metric = "CI", dn = 2,
        bdw = baseline_bdw, esignals = cumulative_esignals,
        label = "tabsimci_supp1_dn2"
      )
    )
  )
}

## ---------------------------------------------------------
## C. Threshold-growth sensitivity in Section sec_addsim2
##
## Supplement order:
##   HS, c_2=0.35
##   HS, c_2=0.25
##   CI, c_2=0.35
##   CI, c_2=0.25
## ---------------------------------------------------------

if (include_rate_sensitivity) {
  table_plan <- c(
    table_plan,
    list(
      make_table_spec(
        family = "baseline", metric = "HS", dn = 2,
        bdw = 0.35, esignals = baseline_esignals,
        show_bdw = TRUE,
        label = "tabsim_supp0_dn2_faster"
      ),
      make_table_spec(
        family = "baseline", metric = "HS", dn = 2,
        bdw = 0.25, esignals = baseline_esignals,
        show_bdw = TRUE,
        label = "tabsim_supp0_dn2_slower"
      ),
      make_table_spec(
        family = "baseline", metric = "CI", dn = 2,
        bdw = 0.35, esignals = baseline_esignals,
        show_bdw = TRUE,
        label = "tabsimci_supp0_dn2_faster"
      ),
      make_table_spec(
        family = "baseline", metric = "CI", dn = 2,
        bdw = 0.25, esignals = baseline_esignals,
        show_bdw = TRUE,
        label = "tabsimci_supp0_dn2_slower"
      )
    )
  )
}

## ---------------------------------------------------------
## D. Alternative DGP in Section sec_addsim3
##
## Supplement order:
##   HS, c_2=0.30
##   HS, c_2=0.35
##   HS, c_2=0.25
##   CI, c_2=0.30
##   CI, c_2=0.35
##   CI, c_2=0.25
## ---------------------------------------------------------

if (include_alternative_dgp) {
  table_plan <- c(
    table_plan,
    list(
      make_table_spec(
        family = "alternative", metric = "HS", dn = 2,
        bdw = 0.30, esignals = alternative_hs_esignals,
        show_bdw = TRUE,
        label = "tabsim_addadd_supp0_dn2"
      ),
      make_table_spec(
        family = "alternative", metric = "HS", dn = 2,
        bdw = 0.35, esignals = alternative_hs_esignals,
        show_bdw = TRUE,
        label = "tabsim_addadd_supp0_dn2a"
      ),
      make_table_spec(
        family = "alternative", metric = "HS", dn = 2,
        bdw = 0.25, esignals = alternative_hs_esignals,
        show_bdw = TRUE,
        label = "tabsim_addadd_supp0_dn2b"
      ),
      make_table_spec(
        family = "alternative", metric = "CI", dn = 2,
        bdw = 0.30, esignals = alternative_ci_esignals,
        show_bdw = TRUE,
        label = "tabsimci_addadd_supp0_dn2"
      ),
      make_table_spec(
        family = "alternative", metric = "CI", dn = 2,
        bdw = 0.35, esignals = alternative_ci_esignals,
        show_bdw = TRUE,
        label = "tabsimci_supp0_dn2a"
      ),
      make_table_spec(
        family = "alternative", metric = "CI", dn = 2,
        bdw = 0.25, esignals = alternative_ci_esignals,
        show_bdw = TRUE,
        label = "tabsimci_supp0_dn2b"
      )
    )
  )
}

## ---------------------------------------------------------
## E. Optional tables not currently reported in the Supplement
##
## These are appended only when explicitly requested, so they
## do not interrupt the Supplement's reported table sequence.
## ---------------------------------------------------------

if (include_unreported_dn3_tables) {
  extra_plan <- list()

  if (include_cumulative_main) {
    extra_plan <- c(
      extra_plan,
      list(
        make_table_spec(
          family = "cumulative", metric = "HS", dn = 3,
          bdw = baseline_bdw, esignals = cumulative_esignals,
          label = "tabsim_supp1_dn3"
        ),
        make_table_spec(
          family = "cumulative", metric = "CI", dn = 3,
          bdw = baseline_bdw, esignals = cumulative_esignals,
          label = "tabsimci_supp1_dn3"
        )
      )
    )
  }

  if (include_rate_sensitivity) {
    for (bdw in c(0.35, 0.25)) {
      extra_plan <- c(
        extra_plan,
        list(
          make_table_spec(
            family = "baseline", metric = "HS", dn = 3,
            bdw = bdw, esignals = baseline_esignals,
            show_bdw = TRUE
          )
        )
      )
    }

    for (bdw in c(0.35, 0.25)) {
      extra_plan <- c(
        extra_plan,
        list(
          make_table_spec(
            family = "baseline", metric = "CI", dn = 3,
            bdw = bdw, esignals = baseline_esignals,
            show_bdw = TRUE
          )
        )
      )
    }
  }

  if (include_alternative_dgp) {
    for (bdw in c(0.30, 0.35, 0.25)) {
      extra_plan <- c(
        extra_plan,
        list(
          make_table_spec(
            family = "alternative", metric = "HS", dn = 3,
            bdw = bdw, esignals = alternative_hs_esignals,
            show_bdw = TRUE
          )
        )
      )
    }

    for (bdw in c(0.30, 0.35, 0.25)) {
      extra_plan <- c(
        extra_plan,
        list(
          make_table_spec(
            family = "alternative", metric = "CI", dn = 3,
            bdw = bdw, esignals = alternative_ci_esignals,
            show_bdw = TRUE
          )
        )
      )
    }
  }

  table_plan <- c(table_plan, extra_plan)
}

## Optional selected-K tables are placed at the end because they
## are not standalone tables in the current Supplement.
if (include_selected_K_tables) {
  selected_K_plan <- list()

  for (dn in c(2, 3)) {
    for (bdw in c(0.30, 0.35, 0.25)) {
      selected_K_plan <- c(
        selected_K_plan,
        list(
          make_table_spec(
            family = "baseline", metric = "K", dn = dn,
            bdw = bdw, esignals = baseline_esignals,
            show_bdw = TRUE
          )
        )
      )
    }
  }

  for (dn in c(2, 3)) {
    selected_K_plan <- c(
      selected_K_plan,
      list(
        make_table_spec(
          family = "cumulative", metric = "K", dn = dn,
          bdw = baseline_bdw, esignals = cumulative_esignals,
          show_bdw = TRUE
        )
      )
    )
  }

  for (dn in c(2, 3)) {
    for (bdw in c(0.30, 0.35, 0.25)) {
      selected_K_plan <- c(
        selected_K_plan,
        list(
          make_table_spec(
            family = "alternative", metric = "K", dn = dn,
            bdw = bdw, esignals = alternative_hs_esignals,
            show_bdw = TRUE
          )
        )
      )
    }
  }

  table_plan <- c(table_plan, selected_K_plan)
}

############################################################
## 3. Utility functions
############################################################

number_token <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, digits = 15)
}

bdw_tag <- function(x) {
  gsub("\\.", "", sprintf("%.2f", x))
}

latex_escape_text <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x
}

comm_wrap <- function(x) {
  if (use_revision_markup) {
    paste0("\\commRV{", x, "}")
  } else {
    x
  }
}

metric_file_code <- function(metric) {
  if (metric == "K") "HS" else metric
}

family_stem <- function(family) {
  switch(
    family,
    baseline = "FRNENEWBLOCKKAPP",
    alternative = "FRNENEWBLOCKSMALLKAPP",
    cumulative = "FRNENEW_KsuppadjBLOCKKAPP",
    stop("Unknown family: ", family)
  )
}

build_file_basenames <- function(spec, esig, rho) {
  family <- spec$family
  metric_code <- metric_file_code(spec$metric)
  bdw <- number_token(spec$bdw)
  esig <- number_token(esig)
  rho <- number_token(rho)
  dn <- number_token(spec$dn)

  if (family == "cumulative") {
    core <- paste0(
      family_stem(family),
      "_esignal", esig,
      "_rho", rho,
      "_Kscheme", number_token(Kscheme),
      "_Kcutmin", number_token(Kval_cut),
      "_dn", dn,
      "_", metric_code,
      "_bdw", bdw
    )

    return(unique(c(
      paste0(core, "_defac1", number_token(decfac1), ".RData"),
      paste0(core, ".RData")
    )))
  }

  core <- paste0(
    family_stem(family),
    "_esignal", esig,
    "_rho", rho,
    "_dn", dn,
    "_", metric_code,
    "_bdw", bdw,
    "_cuteigen", number_token(cuteigen)
  )

  candidates <- c(
    paste0(core, "_defac1", number_token(decfac1), ".RData"),
    paste0(core, ".RData")
  )

  ## Optional future naming convention that records decfac2 explicitly.
  if (family == "alternative") {
    candidates <- c(
      paste0(
        core,
        "_defac1", number_token(decfac1),
        "_defac2", number_token(alternative_decfac2),
        ".RData"
      ),
      candidates
    )
  }

  unique(candidates)
}

resolve_result_file <- function(spec, esig, rho) {
  basenames <- build_file_basenames(spec, esig, rho)
  candidates <- unlist(
    lapply(result_dirs, function(directory) file.path(directory, basenames)),
    use.names = FALSE
  )
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    return(list(path = NA_character_, duplicates = character(0)))
  }

  if (length(existing) > 1L) {
    warning(
      "Multiple matching files found; using the first one:\n",
      paste(existing, collapse = "\n")
    )
  }

  list(path = existing[1L], duplicates = existing[-1L])
}

required_object_names <- function(metric) {
  prefix <- if (metric == "K") "RESULTK" else "RESULT"
  c(
    paste0(prefix, seq_along(T_vals)),
    paste0(prefix, seq_along(T_vals), "_sp")
  )
}

load_result_matrices <- function(path, metric) {
  result_env <- new.env(parent = emptyenv())
  load(path, envir = result_env)

  required <- required_object_names(metric)
  missing <- required[
    !vapply(required, exists, logical(1), envir = result_env, inherits = FALSE)
  ]

  if (length(missing) > 0L) {
    stop(
      "The following objects are missing from ", basename(path), ": ",
      paste(missing, collapse = ", ")
    )
  }

  prefix <- if (metric == "K") "RESULTK" else "RESULT"

  exponential <- lapply(
    seq_along(T_vals),
    function(j) get(paste0(prefix, j), envir = result_env, inherits = FALSE)
  )

  sparse <- lapply(
    seq_along(T_vals),
    function(j) get(paste0(prefix, j, "_sp"), envir = result_env, inherits = FALSE)
  )

  all_matrices <- c(exponential, sparse)

  if (!all(vapply(all_matrices, is.matrix, logical(1)))) {
    stop("At least one required object is not a matrix in ", basename(path))
  }

  if (any(vapply(all_matrices, ncol, integer(1)) < length(kappa_vals))) {
    stop(
      "At least one matrix has fewer than ", length(kappa_vals),
      " kappa columns in ", basename(path)
    )
  }

  list(
    exponential = exponential,
    sparse = sparse,
    n_iter = min(vapply(all_matrices, nrow, integer(1)))
  )
}

summarize_vector <- function(x, metric) {
  if (metric == "HS") {
    return(mean(sqrt(x), na.rm = TRUE))
  }

  if (metric %in% c("CI", "K")) {
    return(mean(x, na.rm = TRUE))
  }

  stop("Unknown metric: ", metric)
}

extract_table_block <- function(path, metric) {
  matrices <- load_result_matrices(path, metric)
  out <- matrix(
    NA_real_,
    nrow = length(kappa_vals),
    ncol = 2L * length(T_vals)
  )

  for (k in seq_along(kappa_vals)) {
    out[k, seq_along(T_vals)] <- vapply(
      matrices$exponential,
      function(mat) summarize_vector(mat[, k], metric),
      numeric(1)
    )

    out[k, length(T_vals) + seq_along(T_vals)] <- vapply(
      matrices$sparse,
      function(mat) summarize_vector(mat[, k], metric),
      numeric(1)
    )
  }

  list(values = out, n_iter = matrices$n_iter)
}

format_entry <- function(x, metric) {
  if (is.na(x) || !is.finite(x)) {
    return("--")
  }

  digits <- if (metric == "K") 2L else 3L
  sprintf(paste0("%.", digits, "f"), x)
}

default_caption <- function(spec) {
  metric_text <- switch(
    spec$metric,
    HS = "Hilbert--Schmidt norms of $\\hat{f}_{\\kappa}-f$",
    CI = "Empirical coverage probabilities",
    K = "Average selected values of $\\KK$"
  )

  family_text <- switch(
    spec$family,
    baseline = "",
    cumulative = " with $\\KK$ selected by the cumulative-variation criterion",
    alternative = " under the alternative DGP"
  )

  bdw_text <- if (isTRUE(spec$show_bdw)) {
    paste0(" with $c_2=", number_token(spec$bdw), "$")
  } else {
    ""
  }

  section_text <- if (spec$family == "cumulative") {
    " (for Section~\\ref{sec_addsim1})"
  } else if (spec$family == "alternative") {
    " (for Section~\\ref{sec_addsim3})"
  } else if (isTRUE(spec$show_bdw)) {
    " (for Section~\\ref{sec_addsim2})"
  } else {
    " (for Section~\\ref{Sec_monte_carlo2})"
  }

  paste0(
    metric_text,
    family_text,
    bdw_text,
    " with $d_N=", spec$dn, "$",
    section_text
  )
}

default_label <- function(spec) {
  if (spec$family == "baseline" &&
      abs(spec$bdw - baseline_bdw) < 1e-12 &&
      !isTRUE(spec$show_bdw)) {
    prefix <- if (spec$metric == "CI") "tabsimci" else if (spec$metric == "K") "tabsimK" else "tabsim"
    return(paste0(prefix, "_supp0_dn", spec$dn))
  }

  if (spec$family == "cumulative" &&
      abs(spec$bdw - baseline_bdw) < 1e-12 &&
      !isTRUE(spec$show_bdw)) {
    prefix <- if (spec$metric == "CI") "tabsimci" else if (spec$metric == "K") "tabsimK" else "tabsim"
    return(paste0(prefix, "_supp1_dn", spec$dn))
  }

  metric_tag <- switch(spec$metric, HS = "hs", CI = "ci", K = "K")
  paste0(
    "tabsim_", metric_tag,
    "_", spec$family,
    "_bdw_", bdw_tag(spec$bdw),
    "_dn_", spec$dn
  )
}

table_note <- function(metric) {
  error_note <- paste0(
    "The scale of measurement errors is given by $\\mathfrak{a}_{e}\\%$ ",
    "when $\\sigma_e$ is chosen so that the nuclear norm of the covariance ",
    "operator of $e_t$ matches $\\mathfrak{a}_{e}\\%$ of that of ",
    "$\\mathcal E_t^x=(\\Delta\\PP^N x_t,\\PP^S x_t)$."
  )

  switch(
    metric,
    HS = paste0(
      "Notes: The average Hilbert--Schmidt norm of ",
      "$\\hat{f}_{\\kappa}-f$ is computed from Monte Carlo replications. ",
      error_note
    ),
    CI = paste0(
      "Notes: Empirical coverage probabilities of the nominal $95\\%$ ",
      "confidence intervals are computed from Monte Carlo replications. ",
      error_note
    ),
    K = paste0(
      "Notes: The entries report Monte Carlo averages of the selected ",
      "truncation parameter $\\KK$. ",
      error_note
    )
  )
}

############################################################
## 4. Output helpers
############################################################

tex_connection <- file(output_tex, open = "wt", encoding = "UTF-8")

emit <- function(...) {
  text <- paste0(...)
  cat(text, file = tex_connection)
  if (print_tables_to_console) {
    cat(text)
  }
}

manifest <- data.frame(
  table_number = integer(0),
  label = character(0),
  family = character(0),
  metric = character(0),
  dn = integer(0),
  bdw = numeric(0),
  esignal = numeric(0),
  rho = numeric(0),
  file = character(0),
  n_iter = integer(0),
  status = character(0),
  stringsAsFactors = FALSE
)

append_manifest <- function(table_number, label, spec, esig, rho,
                            file, n_iter, status) {
  manifest <<- rbind(
    manifest,
    data.frame(
      table_number = table_number,
      label = label,
      family = spec$family,
      metric = spec$metric,
      dn = spec$dn,
      bdw = spec$bdw,
      esignal = esig,
      rho = rho,
      file = ifelse(is.na(file), "", normalizePath(file, winslash = "/", mustWork = FALSE)),
      n_iter = n_iter,
      status = status,
      stringsAsFactors = FALSE
    )
  )
}

############################################################
## 5. LaTeX table writer
############################################################

write_latex_table <- function(spec, table_number) {
  caption <- if (is.null(spec$caption)) default_caption(spec) else spec$caption
  label <- if (is.null(spec$label)) default_label(spec) else spec$label

  emit("\n% ==========================================================\n")
  emit("% Table ", table_number, ": ", label, "\n")
  emit("% ==========================================================\n")
  emit("\\begin{table}[", table_position, "]\n")
  emit("\t\\renewcommand{\\arraystretch}{", array_stretch, "}\n")
  emit("\t\\caption{", comm_wrap(caption), "}\n")
  emit("\t\\label{", label, "}\n")
  emit("\t\\vskip -5pt\n")

  if (nzchar(table_font_command)) {
    emit("\t", table_font_command, "\n")
  }

  if (use_revision_markup) {
    emit("\t\\commRV{\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}} ccc rrrrr rrrrr }\n")
  } else {
    emit("\t\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}} ccc rrrrr rrrrr }\n")
  }

  emit("\t\t\\toprule\n")
  emit("\t\t& & & \\multicolumn{5}{c}{Exponential design} & ",
       "\\multicolumn{5}{c}{Sparse design} \\\\\n")
  emit("\t\t\\cmidrule(lr){4-8} \\cmidrule(l){9-13}\n")
  emit("\t\t\\multirow{2}{*}{Error} & \\multirow{2}{*}{$\\rho$} & ",
       "\\multirow{2}{*}{$\\kappa$} & \\multicolumn{5}{c}{Sample size ($T$)} & ",
       "\\multicolumn{5}{c}{Sample size ($T$)} \\\\\n")
  emit("\t\t\\cmidrule(lr){4-8} \\cmidrule(l){9-13}\n")

  T_header <- paste(T_vals, collapse = " & ")
  emit("\t\t& & & ", T_header, " & ", T_header, " \\\\\n")
  emit("\t\t\\midrule\n")

  for (e_idx in seq_along(spec$esignals)) {
    esig <- spec$esignals[e_idx]

    for (r_idx in seq_along(rho_list)) {
      rho <- rho_list[r_idx]
      resolved <- resolve_result_file(spec, esig, rho)

      if (is.na(resolved$path)) {
        block <- matrix(
          "--",
          nrow = length(kappa_vals),
          ncol = 2L * length(T_vals)
        )

        append_manifest(
          table_number, label, spec, esig, rho,
          NA_character_, 0L, "missing"
        )

        message(
          "[Missing] family=", spec$family,
          ", metric=", spec$metric,
          ", dn=", spec$dn,
          ", bdw=", number_token(spec$bdw),
          ", esignal=", esig,
          ", rho=", rho
        )
      } else {
        extracted <- tryCatch(
          extract_table_block(resolved$path, spec$metric),
          error = function(e) e
        )

        if (inherits(extracted, "error")) {
          block <- matrix(
            "--",
            nrow = length(kappa_vals),
            ncol = 2L * length(T_vals)
          )

          append_manifest(
            table_number, label, spec, esig, rho,
            resolved$path, 0L,
            paste0("error: ", conditionMessage(extracted))
          )

          warning(
            "Could not process ", resolved$path, ": ",
            conditionMessage(extracted)
          )
        } else {
          block <- apply(
            extracted$values,
            c(1, 2),
            format_entry,
            metric = spec$metric
          )

          append_manifest(
            table_number, label, spec, esig, rho,
            resolved$path, extracted$n_iter, "loaded"
          )

          message(
            "[Loaded] ", basename(resolved$path),
            " | iterations=", extracted$n_iter
          )
        }
      }

      for (k_idx in seq_along(kappa_vals)) {
        prefix <- ""

        if (r_idx == 1L && k_idx == 1L) {
          prefix <- paste0(
            "\\multirow{", length(rho_list) * length(kappa_vals),
            "}{*}{$", esig, "\\%$} & ",
            "\\multirow{", length(kappa_vals),
            "}{*}{$", number_token(rho), "$} & "
          )
        } else if (k_idx == 1L) {
          prefix <- paste0(
            "& \\multirow{", length(kappa_vals),
            "}{*}{$", number_token(rho), "$} & "
          )
        } else {
          prefix <- "& & "
        }

        row_values <- paste(block[k_idx, ], collapse = " & ")

        emit(
          "\t\t", prefix,
          "$", kappa_vals[k_idx], "$ & ",
          row_values,
          " \\\\\n"
        )
      }

      if (r_idx < length(rho_list)) {
        emit("\t\t\\cmidrule{2-13}\n")
      }
    }

    if (e_idx < length(spec$esignals)) {
      emit("\t\t\\midrule\n")
    }
  }

  emit("\t\t\\bottomrule\n")

  if (use_revision_markup) {
    emit("\t\\end{tabular*}}\n")
  } else {
    emit("\t\\end{tabular*}\n")
  }

  emit("\t\\vskip 4pt\n")
  emit("\t{\\footnotesize ", comm_wrap(table_note(spec$metric)), "}\n")
  emit("\\end{table}\n")
}

############################################################
## 6. Residual-based test tables
##
## These functions are separate from the general HS/CI table
## functions because the residual-based test uses only kappa=0
## and stores rejection indicators in the first result column.
############################################################

build_residual_test_basenames <- function(test_type, esig, rho) {
  stem <- switch(
    test_type,
    size = "FRNENEWBLOCKTESTSIZE",
    power = "FRNENEWBLOCKTESTPOWER",
    stop("Unknown residual-test type: ", test_type)
  )

  core <- paste0(
    stem,
    "_esignal", number_token(esig),
    "_rho", number_token(rho),
    "_dn", number_token(residual_test_dn),
    "_HS_bdw", number_token(residual_test_bdw),
    "_cuteigen", number_token(cuteigen),
    "_defac1", number_token(decfac1)
  )

  if (test_type == "size") {
    return(paste0(
      core,
      "_ranintmax", number_token(residual_size_ranintmax),
      ".RData"
    ))
  }

  paste0(core, ".RData")
}

resolve_residual_test_file <- function(test_type, esig, rho) {
  basenames <- build_residual_test_basenames(test_type, esig, rho)
  candidates <- unlist(
    lapply(result_dirs, function(directory) file.path(directory, basenames)),
    use.names = FALSE
  )
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    return(list(path = NA_character_, duplicates = character(0)))
  }

  if (length(existing) > 1L) {
    warning(
      "Multiple matching residual-test files found; using the first one:\n",
      paste(existing, collapse = "\n")
    )
  }

  list(path = existing[1L], duplicates = existing[-1L])
}

load_residual_test_matrices <- function(path) {
  result_env <- new.env(parent = emptyenv())
  load(path, envir = result_env)

  required <- c(
    paste0("RESULT", seq_along(T_vals)),
    paste0("RESULT", seq_along(T_vals), "_sp")
  )

  missing <- required[
    !vapply(required, exists, logical(1), envir = result_env, inherits = FALSE)
  ]

  if (length(missing) > 0L) {
    stop(
      "The following objects are missing from ", basename(path), ": ",
      paste(missing, collapse = ", ")
    )
  }

  exponential <- lapply(
    seq_along(T_vals),
    function(j) get(paste0("RESULT", j), envir = result_env, inherits = FALSE)
  )

  sparse <- lapply(
    seq_along(T_vals),
    function(j) get(paste0("RESULT", j, "_sp"), envir = result_env, inherits = FALSE)
  )

  all_matrices <- c(exponential, sparse)

  if (!all(vapply(all_matrices, is.matrix, logical(1)))) {
    stop("At least one required residual-test object is not a matrix in ",
         basename(path))
  }

  if (any(vapply(all_matrices, ncol, integer(1)) < 1L)) {
    stop("At least one residual-test matrix has no kappa=0 column in ",
         basename(path))
  }

  list(
    exponential = exponential,
    sparse = sparse,
    n_iter = min(vapply(all_matrices, nrow, integer(1)))
  )
}

extract_residual_test_block <- function(path) {
  matrices <- load_residual_test_matrices(path)

  values <- c(
    vapply(
      matrices$exponential,
      function(mat) mean(mat[, 1L], na.rm = TRUE),
      numeric(1)
    ),
    vapply(
      matrices$sparse,
      function(mat) mean(mat[, 1L], na.rm = TRUE),
      numeric(1)
    )
  )

  list(values = values, n_iter = matrices$n_iter)
}

residual_test_caption <- function(test_type) {
  switch(
    test_type,
    size = "Size of the residual-based test",
    power = "Power of the residual-based test",
    stop("Unknown residual-test type: ", test_type)
  )
}

residual_test_label <- function(test_type) {
  switch(
    test_type,
    size = "tabsim_supp1_dn2_residual1",
    power = "tabsim_supp1_dn2_residual2",
    stop("Unknown residual-test type: ", test_type)
  )
}

residual_test_note <- paste0(
  "Notes: Rejection rates of the residual-based test at the $5\\%$ ",
  "nominal level over 3,000 Monte Carlo replications. The scale of ",
  "measurement errors is given by $\\mathfrak{a}_e\\%$ when $\\sigma_e$ is chosen ",
  "so that the nuclear norm of the covariance operator of $e_t$ matches ",
  "$\\mathfrak{a}_e\\%$ of that of ",
  "$\\mathcal E_t^x=(\\Delta\\PP^N x_t,\\PP^S x_t)$."
)

write_residual_test_table <- function(test_type, table_number) {
  caption <- residual_test_caption(test_type)
  label <- residual_test_label(test_type)

  emit("\n% ==========================================================\n")
  emit("% Table ", table_number, ": ", label, "\n")
  emit("% ==========================================================\n")
  emit("\\begin{table}[H]\n")
  emit("\t\\renewcommand{\\arraystretch}{0.85}\n")
  emit("\t\\caption{", comm_wrap(caption), "}\n")
  emit("\t\\label{", label, "}\n")
  emit("\t\\vskip -5pt\n")
  emit("\t\\small\n")

  if (use_revision_markup) {
    emit("\t\\commRV{\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}} ccc rrrrr rrrrr }\n")
  } else {
    emit("\t\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}} ccc rrrrr rrrrr }\n")
  }

  emit("\t\t\\toprule\n")
  emit("\t\t& & & \\multicolumn{5}{c}{Exponential design} & ",
       "\\multicolumn{5}{c}{Sparse design} \\\\\n")
  emit("\t\t\\cmidrule(lr){4-8} \\cmidrule(l){9-13}\n")
  emit("\t\t\\multirow{2}{*}{Error} & \\multirow{2}{*}{$\\rho$} & ",
       "\\multirow{2}{*}{$\\kappa$} & \\multicolumn{5}{c}{Sample size ($T$)} & ",
       "\\multicolumn{5}{c}{Sample size ($T$)} \\\\\n")
  emit("\t\t\\cmidrule(lr){4-8} \\cmidrule(l){9-13}\n")

  T_header <- paste(T_vals, collapse = " & ")
  emit("\t\t& & & ", T_header, " & ", T_header, " \\\\\n")
  emit("\t\t\\midrule\n")

  manifest_spec <- list(
    family = "residual_test",
    metric = toupper(test_type),
    dn = residual_test_dn,
    bdw = residual_test_bdw
  )

  for (e_idx in seq_along(residual_test_esignals)) {
    esig <- residual_test_esignals[e_idx]

    for (r_idx in seq_along(residual_test_rhos)) {
      rho <- residual_test_rhos[r_idx]
      resolved <- resolve_residual_test_file(test_type, esig, rho)

      if (is.na(resolved$path)) {
        formatted <- rep("--", 2L * length(T_vals))

        append_manifest(
          table_number, label, manifest_spec, esig, rho,
          NA_character_, 0L, "missing"
        )

        message(
          "[Missing] residual test=", test_type,
          ", dn=", residual_test_dn,
          ", bdw=", number_token(residual_test_bdw),
          ", esignal=", esig,
          ", rho=", rho
        )
      } else {
        extracted <- tryCatch(
          extract_residual_test_block(resolved$path),
          error = function(e) e
        )

        if (inherits(extracted, "error")) {
          formatted <- rep("--", 2L * length(T_vals))

          append_manifest(
            table_number, label, manifest_spec, esig, rho,
            resolved$path, 0L,
            paste0("error: ", conditionMessage(extracted))
          )

          warning(
            "Could not process ", resolved$path, ": ",
            conditionMessage(extracted)
          )
        } else {
          formatted <- sprintf("%.3f", extracted$values)

          append_manifest(
            table_number, label, manifest_spec, esig, rho,
            resolved$path, extracted$n_iter, "loaded"
          )

          message(
            "[Loaded] ", basename(resolved$path),
            " | iterations=", extracted$n_iter
          )
        }
      }

      if (r_idx == 1L) {
        prefix <- paste0(
          "\\multirow{", length(residual_test_rhos),
          "}{*}{$", esig, "\\%$} & ",
          "$", number_token(rho), "$ & $0$ & "
        )
      } else {
        prefix <- paste0(
          "& $", number_token(rho), "$ & $0$ & "
        )
      }

      emit(
        "\t\t", prefix,
        paste(formatted, collapse = " & "),
        " \\\\\n"
      )

      if (r_idx < length(residual_test_rhos)) {
        emit("\t\t\\cmidrule{2-13}\n")
      }
    }

    if (e_idx < length(residual_test_esignals)) {
      emit("\t\t\\midrule\n")
    }
  }

  emit("\t\t\\bottomrule\n")

  if (use_revision_markup) {
    emit("\t\\end{tabular*}}\n")
  } else {
    emit("\t\\end{tabular*}\n")
  }

  emit("\t\\vskip 4pt\n")

  if (use_revision_markup) {
    emit("\t\\commRV{\\footnotesize ", residual_test_note, "}\n")
  } else {
    emit("\t{\\footnotesize ", residual_test_note, "}\n")
  }

  emit("\\end{table}\n")
}

############################################################
## 7. Generate all tables
############################################################

cat("\n============================================================\n")
cat("FRNE simulation table generation\n")
cat("============================================================\n")
cat("Output file: ", output_tex, "\n", sep = "")
n_residual_tables <- if (include_residual_test_tables) 2L else 0L
total_requested_tables <- length(table_plan) + n_residual_tables

cat("Number of requested tables: ", total_requested_tables, "\n", sep = "")

for (i in seq_along(table_plan)) {
  write_latex_table(table_plan[[i]], i)
}

if (include_residual_test_tables) {
  residual_start <- length(table_plan) + 1L
  write_residual_test_table("size", residual_start)
  write_residual_test_table("power", residual_start + 1L)
}

close(tex_connection)
write.csv(manifest, output_manifest, row.names = FALSE, na = "")

cat("\n============================================================\n")
cat("Completed\n")
cat("============================================================\n")
cat("LaTeX tables: ", output_tex, "\n", sep = "")
cat("Source manifest: ", output_manifest, "\n", sep = "")
cat("Loaded files: ", sum(manifest$status == "loaded"), "\n", sep = "")
cat("Missing files: ", sum(manifest$status == "missing"), "\n", sep = "")
cat("Processing errors: ", sum(grepl("^error:", manifest$status)), "\n", sep = "")



