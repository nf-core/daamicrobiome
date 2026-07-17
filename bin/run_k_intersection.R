#!/usr/bin/env Rscript
# ======================================================================
# K-intersection consensus for real DA results (Path A output)
#
# For each k = 1..n_tools, compute taxa detected by at least k tools
# in the same direction. Output an XLSX workbook (one sheet per k)
# and individual CSV files.
# ======================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
})

read_any <- function(p) {
  ext <- tolower(tools::file_ext(p))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  utils::read.table(p, header = TRUE, sep = sep, check.names = FALSE,
                    quote = "", comment.char = "")
}

main <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "Text file: one *_da.tsv path per line"),
    make_option(c("-o", "--outdir"), type = "character", default = ".",
                help = "Output directory [default %default]"),
    make_option(c("-a", "--alpha"), type = "numeric", default = 0.05,
                help = "Adj. p-value threshold [default %default]"),
    make_option(c("-l", "--lfc_min"), type = "numeric", default = 0,
                help = "Minimum |log2 fold-change| [default %default]")
  )
  args <- parse_args(OptionParser(option_list = option_list))
  if (is.null(args$input)) stop("--input is required.")

  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

  lines <- readLines(args$input, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !grepl("^\\s*#", lines)]
  lines <- sub("^['\"](.*)['\"]$", "\\1", lines)
  files <- lines
  if (length(files) == 0) stop("No files listed in --input.")

  req <- c("taxon_id", "lfc", "pvalue", "adj_pval", "tool")
  all_rows <- list()
  for (f in files) {
    if (!file.exists(f)) stop("File not found: ", f)
    d <- read_any(f)
    miss <- setdiff(req, names(d))
    if (length(miss)) stop("File ", basename(f), " missing columns: ", paste(miss, collapse = ", "))
    d$taxon_id <- as.character(d$taxon_id)
    suppressWarnings({
      d$lfc      <- as.numeric(d$lfc)
      d$adj_pval <- as.numeric(d$adj_pval)
      d$tool     <- as.character(d$tool)
    })
    all_rows[[length(all_rows) + 1]] <- d
  }
  all_df <- do.call(rbind, all_rows)

  if (is.null(all_df) || nrow(all_df) == 0) {
    message("No data rows across input files. Writing empty outputs.")
    wb <- createWorkbook()
    addWorksheet(wb, "k1")
    saveWorkbook(wb, file.path(args$outdir, "k_intersection.xlsx"), overwrite = TRUE)
    return(invisible(NULL))
  }

  all_df$called <- !is.na(all_df$adj_pval) & all_df$adj_pval < args$alpha &
    !is.na(all_df$lfc) & abs(all_df$lfc) >= args$lfc_min

  all_df$direction <- ifelse(is.na(all_df$lfc), NA_character_,
                             ifelse(all_df$lfc > 0, "up", "down"))

  called_df <- all_df[all_df$called & !is.na(all_df$direction), ]

  tool_names <- sort(unique(all_df$tool))
  n_tools <- length(tool_names)
  cat("[k_intersection] Tools:", paste(tool_names, collapse = ", "), "\n")
  cat("[k_intersection] n_tools:", n_tools, "\n")

  # For each taxon: count tools calling it in each direction (strict = same direction)
  pos_agg <- aggregate(tool ~ taxon_id, data = called_df[called_df$direction == "up", ],
                       FUN = function(x) list(sort(unique(x))))
  neg_agg <- aggregate(tool ~ taxon_id, data = called_df[called_df$direction == "down", ],
                       FUN = function(x) list(sort(unique(x))))
  names(pos_agg) <- c("taxon_id", "tools_up")
  names(neg_agg) <- c("taxon_id", "tools_dn")

  all_taxa <- sort(unique(called_df$taxon_id))
  taxon_tbl <- data.frame(taxon_id = all_taxa, stringsAsFactors = FALSE)

  taxon_tbl <- merge(taxon_tbl, pos_agg, by = "taxon_id", all.x = TRUE)
  taxon_tbl <- merge(taxon_tbl, neg_agg, by = "taxon_id", all.x = TRUE)

  taxon_tbl$n_up <- sapply(taxon_tbl$tools_up, function(x) if (is.null(x)) 0L else length(x))
  taxon_tbl$n_dn <- sapply(taxon_tbl$tools_dn, function(x) if (is.null(x)) 0L else length(x))

  # Majority-direction: assign each taxon to the direction supported by more tools.
  # Ties (equal up and down) are dropped as ambiguous.
  taxon_tbl$direction <- ifelse(
    taxon_tbl$n_up > taxon_tbl$n_dn, "up",
    ifelse(taxon_tbl$n_dn > taxon_tbl$n_up, "down", "ambiguous")
  )
  taxon_tbl$support <- ifelse(taxon_tbl$direction == "up", taxon_tbl$n_up,
                              ifelse(taxon_tbl$direction == "down", taxon_tbl$n_dn, 0L))

  n_ambig <- sum(taxon_tbl$direction == "ambiguous")
  if (n_ambig > 0)
    cat("[k_intersection]", n_ambig, "taxa dropped (equal support in both directions)\n")

  taxon_tbl$tools_list_up <- sapply(taxon_tbl$tools_up, function(x) {
    if (is.null(x)) "" else paste(x, collapse = ";")
  })
  taxon_tbl$tools_list_dn <- sapply(taxon_tbl$tools_dn, function(x) {
    if (is.null(x)) "" else paste(x, collapse = ";")
  })

  # Build per-tool detection columns
  for (t in tool_names) {
    taxon_tbl[[paste0("det_", t)]] <- sapply(seq_len(nrow(taxon_tbl)), function(i) {
      up_tools <- if (!is.null(taxon_tbl$tools_up[[i]])) taxon_tbl$tools_up[[i]] else character(0)
      dn_tools <- if (!is.null(taxon_tbl$tools_dn[[i]])) taxon_tbl$tools_dn[[i]] else character(0)
      as.integer(t %in% c(up_tools, dn_tools))
    })
  }

  # Build XLSX and CSVs for each k
  wb <- createWorkbook()
  for (k in seq_len(n_tools)) {
    sheet_df <- taxon_tbl[taxon_tbl$direction %in% c("up", "down") & taxon_tbl$support >= k, ]
    out_df <- data.frame(
      taxon_id = sheet_df$taxon_id,
      direction = sheet_df$direction,
      n_tools_agreeing = sheet_df$support,
      tools_up = sheet_df$tools_list_up,
      tools_down = sheet_df$tools_list_dn,
      stringsAsFactors = FALSE
    )
    det_cols <- grep("^det_", names(sheet_df), value = TRUE)
    if (length(det_cols) > 0) out_df <- cbind(out_df, sheet_df[, det_cols, drop = FALSE])
    out_df <- out_df[order(-out_df$n_tools_agreeing, out_df$taxon_id), ]

    sheet_name <- paste0("k", k)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, out_df)

    csv_path <- file.path(args$outdir, paste0("k", k, "_intersection.csv"))
    utils::write.csv(out_df, csv_path, row.names = FALSE)
    cat("[k_intersection] k=", k, ": ", nrow(out_df), " taxa\n", sep = "")
  }

  xlsx_path <- file.path(args$outdir, "k_intersection.xlsx")
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("[k_intersection] Wrote:", xlsx_path, "\n")
}

if (sys.nframe() == 0) main()
