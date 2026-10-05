args <- commandArgs()
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- dirname(normalizePath(script))
trailing <- commandArgs(trailingOnly = TRUE)
output <- if (length(trailing)) trailing[1] else file.path(root, "results")
expected <- c("metrics.csv", "overview.png", "predictions.csv.gz", "selected_features.csv.gz", "splits.csv.gz", "tuning.csv.gz", "summary.md", "run_info.txt")
files <- file.path(output, expected)
stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
m <- read.csv(file.path(output, "metrics.csv"))
stopifnot(all(c("repeat_id","method","balanced_accuracy") %in% names(m)),
          nrow(m) > 0, all(is.finite(m$balanced_accuracy)),
          all(m$balanced_accuracy >= 0 & m$balanced_accuracy <= 1))
message("Verified benchmark outputs and bounded balanced accuracy. Custom runs may differ from the bundled reference.")
