# Replication code for Main Paper Table 1.

setwd("C:/Users/user/Dropbox/Academic_Life(201808-Present)/2_Main_Research/0_Journal_Pub/2026_Submission/202602_Sector_CC/202509_P1_Code/Program/202608_JASA_Revision_Code")  

get_script_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }

  call_frames <- sys.frames()
  for (frame_index in rev(seq_along(call_frames))) {
    source_path <- call_frames[[frame_index]]$ofile
    if (!is.null(source_path)) {
      return(dirname(normalizePath(source_path, winslash = "/", mustWork = TRUE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

required_packages <- c("fda", "geigen", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L))
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required R package(s) before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

code_dir <- get_script_dir()
data_file <- file.path(code_dir, "Data_FCRGT.xlsx")
reference_file <- file.path(code_dir, "Table1_MonteCarlo_Reference.R")
output_file <- file.path(code_dir, "Table1_Replication_Output.csv")

if (!file.exists(data_file)) {
  stop("Replication data file not found: ", data_file)
}
if (!file.exists(reference_file)) {
  stop("Monte Carlo reference file not found: ", reference_file)
}

temperature_grid <- unlist(
  readxl::read_excel(
    data_file,
    sheet = "LTemp",
    range = "B1:CRB1",
    col_names = FALSE,
    col_types = "numeric",
    .name_repair = "minimal"
  ),
  use.names = FALSE
)

temperature_density <- t(as.matrix(readxl::read_excel(
  data_file,
  sheet = "LTemp",
  range = "B2:CRB70",
  col_names = FALSE,
  col_types = "numeric",
  .name_repair = "minimal"
)))

stopifnot(
  length(temperature_grid) == 2497L,
  identical(dim(temperature_density), c(2497L, 69L)),
  all(is.finite(temperature_grid)),
  all(diff(temperature_grid) > 0),
  all(is.finite(temperature_density)),
  all(temperature_density > 0)
)

smax <- 5L
numbasis <- 50L
nobs <- ncol(temperature_density)

log_density <- log(temperature_density)
x_mat <- sweep(log_density, 2L, colMeans(log_density), FUN = "-")

basis_function <- fda::create.bspline.basis(
  rangeval = range(temperature_grid),
  nbasis = numbasis
)

fd_xx <- fda::Data2fd(
  y = x_mat - rowMeans(x_mat),
  argvals = temperature_grid,
  basisobj = basis_function
)
hhtau <- eigen(crossprod(t(fd_xx$coefs)), symmetric = TRUE)

fd_x <- fda::Data2fd(
  y = x_mat,
  argvals = temperature_grid,
  basisobj = basis_function
)
h2kmat <- t(fd_x$coefs - rowMeans(fd_x$coefs))
fd_z <- t(h2kmat %*% hhtau$vectors[, seq_len(numbasis)])
fd_z <- fd_z - rowMeans(fd_z)

kmat <- t(fd_z[seq_len(smax), , drop = FALSE])
cmat <- crossprod(kmat)
partial_sums <- apply(kmat, 2L, cumsum)
smat <- crossprod(partial_sums)
tau <- (nobs^2) * geigen::geigen(
  cmat,
  smat,
  symmetric = TRUE,
  only.values = TRUE
)$values

statistics_ascending <- vapply(
  seq_len(smax),
  function(d0) sum(tau[seq_len(d0)]),
  numeric(1L)
)

monte_carlo <- new.env(parent = baseenv())
sys.source(reference_file, envir = monte_carlo)
stopifnot(
  monte_carlo$table1_monte_carlo_draws == 100000L,
  monte_carlo$table1_brownian_grid_points == 1000L,
  length(monte_carlo$grid) == 1001L,
  identical(dim(monte_carlo$CV), c(5L, 1001L))
)

p_value_labels_ascending <- vapply(seq_len(smax), function(d0) {
  first_quantile <- which(
    monte_carlo$CV[d0, ] >= statistics_ascending[d0]
  )[1L]

  if (is.na(first_quantile) || monte_carlo$grid[first_quantile] >= 0.999) {
    return("<0.1")
  }

  sprintf("%.1f", 100 * (1 - monte_carlo$grid[first_quantile]))
}, character(1L))

table_1 <- data.frame(
  d0 = rev(seq_len(smax)),
  test_statistic = rev(round(statistics_ascending, 2L)),
  p_value_percent = rev(p_value_labels_ascending),
  check.names = FALSE
)

published_statistics <- c(7247.24, 3216.90, 1214.36, 177.39, 11.73)
published_p_value_labels <- c("<0.1", "<0.1", "0.3", "26.1", "87.3")

statistics_match <- identical(
  formatC(table_1$test_statistic, format = "f", digits = 2L),
  formatC(published_statistics, format = "f", digits = 2L)
)
p_values_match <- identical(
  table_1$p_value_percent,
  published_p_value_labels
)

if (!statistics_match || !p_values_match) {
  stop(
    "Table 1 replication check failed. Computed results do not match the ",
    "published statistics and p-values."
  )
}

print(table_1, row.names = FALSE)
write.csv(table_1, output_file, row.names = FALSE, quote = FALSE)
cat("\nSaved: ", output_file, "\n", sep = "")
