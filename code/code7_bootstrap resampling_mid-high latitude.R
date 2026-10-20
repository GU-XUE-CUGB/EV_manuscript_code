# Final site-balanced latitude analysis; both latitude groups are analysed together.
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1])))
python <- Sys.getenv("KPG_PYTHON", unset="python3")
file_path <- file.path(here, "modules/Code07_latitude/Code07_latitude_analysis.py")
status <- system2(python, c(shQuote(file_path), "--no-plots"))
stopifnot(status == 0)
