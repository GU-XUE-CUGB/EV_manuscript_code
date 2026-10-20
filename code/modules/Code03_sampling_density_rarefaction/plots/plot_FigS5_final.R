#!/usr/bin/env Rscript
# Evaluate the exact final S5 build/drawing code from the complete preserved
# unified formatter, stopping before its S6/S7 uncertainty input section.
argv <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
source_path <- file.path(dirname(normalizePath(script)),"plot_FigS4_S7_final_unified_vector.R")
source <- readLines(source_path,warn=FALSE)
end <- grep("^inventory_dir <-",source)[1]-1L
stopifnot(length(end)==1L,is.finite(end),end>1L)
eval(parse(text=source[seq_len(end)]),envir=.GlobalEnv)
export_vector_set("Fig_S5_final",fig_s5,width=10.35,height=12.20)
writeLines("Plot-only final S5: historical formal 10000-iteration A1/A2 summaries; no rarefaction rerun.",
 file.path(output_dir,"PLOT_INPUT_STATUS.txt"))
