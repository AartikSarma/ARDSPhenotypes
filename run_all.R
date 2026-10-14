# run_all.R
#
# Runs the AJRCCM 2025 analysis scripts in order, each in a fresh R session,
# stopping at the first failure. Run from the repository root:
#
#   uvr run run_all.R            run scripts 0 to 7
#   uvr run run_all.R -- 4       resume from script 4
#
# Each script runs with the same R as this one; the project's .Rprofile puts
# the uvr package library first. Script 3 (MAST) takes over an hour. Scripts
# 1, 2 and 5 stop until GEO provides age; see README.md.

script_directory <- file.path("scripts", "ajrccm_2025")
script_numbers <- 0:7

arguments <- commandArgs(trailingOnly = TRUE)
first_script <- if (length(arguments) > 0) as.integer(arguments[1]) else 0L
if (is.na(first_script) || !first_script %in% script_numbers) {
  stop("The first script must be a number from 0 to 7.")
}

rscript <- file.path(R.home("bin"), "Rscript")

for (script_number in script_numbers[script_numbers >= first_script]) {
  script_path <- list.files(script_directory,
                            pattern = paste0("^ARDS Phenotypes - ", script_number, " - .*\\.R$"),
                            full.names = TRUE)
  stopifnot(length(script_path) == 1)
  message("=== ", basename(script_path), " (", format(Sys.time(), "%H:%M:%S"), ")")
  exit_status <- system2(rscript, shQuote(script_path))
  if (exit_status != 0) {
    stop(basename(script_path), " failed (exit status ", exit_status, "); later scripts were not run.")
  }
}

message("=== All scripts finished (", format(Sys.time(), "%H:%M:%S"), ")")
