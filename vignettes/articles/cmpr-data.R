CMPR_RECORD <- "https://reshare.ukdataservice.ac.uk/851454/"
CMPR_ARCHIVE_URL <- paste0(CMPR_RECORD, "2/CMPR_Data.zip")
CMPR_ARCHIVE_SHA256 <- "741029dda3ea1bcda757b9d408c20d72a933e8a21676f47396389df8d2f619b4"
CMPR_DTA_SHA256 <- "9d7886a0366465cee4f054608f0258bc3b84bcf7eab8052f40ae885f37492f90"

.cmpr_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The website build requires the '", package, "' package.", call. = FALSE)
  }
}

.cmpr_sha256 <- function(path) {
  .cmpr_require("digest")
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

.cmpr_verify <- function(path, expected, label) {
  actual <- .cmpr_sha256(path)
  if (!identical(actual, expected)) {
    stop(
      label, " checksum mismatch. Expected ", expected, ", received ", actual,
      ". Delete the file and obtain a fresh copy from the ReShare record.",
      call. = FALSE
    )
  }
  invisible(path)
}

.cmpr_find_dta <- function(path) {
  if (dir.exists(path)) {
    candidates <- list.files(
      path,
      pattern = "^CMPR[.]dta$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(candidates) != 1L) {
      stop("The CMPR directory must contain exactly one CMPR.dta file.", call. = FALSE)
    }
    return(candidates)
  }

  if (!file.exists(path)) {
    stop("CMPR source path does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.cmpr_extract_dta <- function(archive, exdir) {
  listing <- utils::unzip(archive, list = TRUE)$Name
  entry <- listing[tolower(basename(listing)) == "cmpr.dta"]
  if (length(entry) != 1L) {
    stop("The CMPR archive must contain exactly one CMPR.dta file.", call. = FALSE)
  }
  utils::unzip(archive, files = entry, exdir = exdir, junkpaths = TRUE)
  file.path(exdir, basename(entry))
}

cmpr_source_file <- function(
    path = Sys.getenv("CSDM_CMPR_PATH", unset = ""),
    cache_dir = Sys.getenv(
      "CSDM_CMPR_CACHE",
      unset = tools::R_user_dir("csdm", which = "cache")
    )) {
  supplied <- nzchar(path)
  if (supplied) {
    source <- .cmpr_find_dta(path)
  } else {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    source <- file.path(cache_dir, "CMPR_Data.zip")
    if (!file.exists(source)) {
      utils::download.file(CMPR_ARCHIVE_URL, source, mode = "wb", quiet = TRUE)
    }
  }

  if (tolower(tools::file_ext(source)) == "zip") {
    tryCatch(
      .cmpr_verify(source, CMPR_ARCHIVE_SHA256, "CMPR archive"),
      error = function(error) {
        if (!supplied) unlink(source)
        stop(error)
      }
    )
    extract_root <- if (supplied) tempdir() else cache_dir
    extract_dir <- file.path(extract_root, "cmpr-extracted")
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    source <- .cmpr_extract_dta(source, extract_dir)
  }

  if (tolower(tools::file_ext(source)) != "dta") {
    stop("CSDM_CMPR_PATH must identify CMPR.dta, CMPR_Data.zip, or its directory.",
      call. = FALSE)
  }
  .cmpr_verify(source, CMPR_DTA_SHA256, "CMPR.dta")
  source
}

read_cmpr <- function(path = Sys.getenv("CSDM_CMPR_PATH", unset = "")) {
  .cmpr_require("haven")
  data <- as.data.frame(haven::read_dta(cmpr_source_file(path)))
  required <- c("gdebt", "year", "gdp", "cpi", "ccode")
  if (!identical(names(data), required) || nrow(data) != 1840L ||
      length(unique(data$ccode)) != 40L || !identical(range(data$year), c(1965, 2010))) {
    stop("CMPR.dta does not have the expected 40-country, 1965-2010 structure.",
      call. = FALSE)
  }
  data
}

.cmpr_panel_difference <- function(x, id, time) {
  out <- rep(NA_real_, length(x))
  groups <- split(seq_along(x), id)
  for (rows in groups) {
    rows <- rows[order(time[rows])]
    consecutive <- diff(time[rows]) == 1
    change <- diff(x[rows])
    change[!consecutive] <- NA_real_
    out[rows[-1L]] <- change
  }
  out
}

prepare_cmpr <- function(data = read_cmpr()) {
  data <- data[order(data$ccode, data$year), , drop = FALSE]
  rownames(data) <- NULL
  data$y <- log(data$gdp)
  data$dy <- .cmpr_panel_difference(data$y, data$ccode, data$year)
  data$p <- log(data$cpi)
  data$dp <- .cmpr_panel_difference(data$p, data$ccode, data$year)
  data$gd <- log(data$gdebt)
  data$dgd <- .cmpr_panel_difference(data$gd, data$ccode, data$year)
  data
}
