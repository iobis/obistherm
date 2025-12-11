#' Check Copernicus user
#'
#' @param user Copernicus user
#' @param pwd Copernicus password
#'
#' @return nothing, user and password set at Global env
#' @export
#' 
#' @details
#' User and password can be set as environmental variables also through
#' COPERNICUS_USER and COPERNICUS_PWD
#'
#' @examples
#' \dontrun{
#' check_user()
#' }
#'
check_user <- function(user = NULL, pwd = NULL) {
    if (Sys.getenv("COPERNICUS_USER") != "") {
        .user <- Sys.getenv("COPERNICUS_USER")
        .pwd <- Sys.getenv("COPERNICUS_PWD")
    } else {
        if (interactive()) {
            .user <- readline("Enter your user\n")
            .pwd <- readline("Enter your password\n")
        } else {
            if (any(is.null(c(user, pwd)))) {
                stop("You need to supply the copernicus user and password in settings.yml")
            } else {
                .user <- user
                .pwd <- pwd
            }
        }
    }
    .user <<- .user
    .pwd <<- .pwd
    return(invisible())
}

#' Get OBIS dataset from AWS bucket
#'
#' @param obis_source optional path to where to save the files (or where files are)
#'
#' @return file path, with files downloaded
#' @export
#' 
#' @details
#' If path exists and contains file, it will only sync
#'
#' @examples
#' \dontrun{
#' get_obis()
#' }
#'
get_obis <- function(obis_source = NULL) {
    message("Downloading/updating OBIS dataset...")

    if (is.null(obis_source)) {
        f <- "data/obis_data"
    } else {
        f <- obis_source
    }
    fs::dir_create(f)
    # Test if aws is available
    if (!nzchar(Sys.which("aws"))) {
        stop("aws is not available on the system; install it to download OBIS data")
    }
    comm <- paste0("aws s3 sync --no-sign-request s3://obis-open-data/occurrence/ ", f)
    system(comm)

    writeLines(paste("obis_data dataset downloaded on", format(Sys.time()), "using the command", comm,
                    "\n\nCite as: Ocean Biodiversity Information System (OBIS) (25 March 2025) OBIS Occurrence Data. Intergovernmental Oceanographic Commission of UNESCO. https://doi.org/10.25607/obis.occurrence.b89117cd."),
            "data/obis-download.txt")
    return(f)
}

#' Partition OBIS dataset by year
#'
#' @param obis_source path to the OBIS dataset as downloaded by [get_obis()]
#' @param min_year minimum year to keep, the rest will be ignored
#' @param skip if TRUE, it will only return the folder (if you already divided and want to avoid doing again)
#'
#' @return saved files in the folder using Hive format
#' @export
#' 
#' @details
#' Files are saves as .parquet using hive naming schema
#'
#' @examples
#' \dontrun{
#' partition_by_year()
#' }
#'
partition_by_year <- function(obis_source, min_year = 1982, skip = FALSE) {

    message("Partitioning dataset by year")
    require(duckdb)

    f <- dirname(obis_source) |>
        file.path("obis_partitioned")
    
    if (skip) {
        message("Skipping data preparing and using old version at ", f)
        return(f)
    }

    if (dir.exists(f)) {
        message("Cleaning folder...")
        fs::dir_delete(f)
    }

    con <- DBI::dbConnect(duckdb())

    DBI::dbSendQuery(con, paste0("
    COPY(
        SELECT _id, dataset_id, source.occurrenceID as occurrenceID, source.datasetID as datasetID,
            interpreted.aphiaid as AphiaID, interpreted.scientificName as scientificName,
            interpreted.species as species, interpreted.genus as genus, interpreted.family as family,
            interpreted.order as order, interpreted.class as class, interpreted.phylum as phylum,
            interpreted.kingdom as kingdom,
            interpreted.eventDate as eventDate, interpreted.date_start as date_start,
            interpreted.date_mid as date_mid, interpreted.date_end as date_end, interpreted.date_year as date_year,
            interpreted.decimalLongitude as decimalLongitude, interpreted.decimalLatitude as decimalLatitude,
            interpreted.coordinatePrecision as coordinatePrecision, interpreted.coordinateUncertaintyInMeters as coordinateUncertaintyInMeters,
            interpreted.minimumDepthInMeters as minimumDepthInMeters, interpreted.maximumDepthInMeters as maximumDepthInMeters,
            dropped, absence, flags
        FROM read_parquet('/Volumes/OBIS2/data/obis_data/*')
        WHERE dropped IS NOT TRUE AND date_year >= ", min_year, "
    ) TO '/Volumes/OBIS2/data/obis_partitioned'
    (FORMAT parquet, PARTITION_BY (date_year));
    "))

    DBI::dbDisconnect(con)

    return(f)
}

#' Start a DuckDB connection with the OBIS dataset
#'
#' @param source_folder folder containing the OBIS dataset divided by year (using [partition_by_year()])
#'
#' @return A DBI connection
#' @export
#' 
#' @details
#' A temporary view of the files is built in memory
#'
#' @examples
#' \dontrun{
#' db_con <- open_db()
#' }
#'
open_db <- function(source_folder) {
    require(duckdb)

    con <- DBI::dbConnect(duckdb())

    DBI::dbSendQuery(con, paste0(
    "
    CREATE VIEW obis AS SELECT * FROM read_parquet('", file.path(source_folder, "*/*"), "');
    "
    ))

    return(con)
}

#' Check depth difference between original and new
#'
#' @param depth_original vector of original depths
#' @param depth_new vector of new depths
#' @param limit maximum limit to return TRUE
#'
#' @return what_return
#' @export
#' 
#' @details
#' if NA, it will return FALSE
#'
#' @examples
#' \dontrun{
#' diff_depths <- check_depth_diff(original_depth, new_depth)
#' }
#'
check_depth_diff <- function(depth_original, depth_new, limit = 5) {
    diff_depth <- abs(depth_original - depth_new)
    diff_depth <- ifelse(diff_depth > limit, TRUE, FALSE)
    ifelse(is.na(diff_depth), FALSE, diff_depth)
}

#' Start Dask connection
#'
#' @param browse if TRUE, open the browser for monitoring
#'
#' @return the Dask client
#' @export
#' 
#' @details
#' Dask enables parallel processing of large datasets (through Python)
#'
#' @examples
#' \dontrun{
#' start_dask()
#' }
#'
start_dask <- function(browse = TRUE) {
    da <- import("dask")
    dd <- import("dask.distributed")
    os <- import("os")
    os$environ["BOKEH_SESSION_TOKEN_EXPIRATION"] <- "86400"
    client <- dd$Client()
    if (browse) {
        browseURL("http://localhost:8787/status")
    } else {
        cat("For browsing, access http://localhost:8787/status\n")
    }
    return(client)
}

#' Decode flags
#'
#' @param flag vector of flags
#'
#' @return decoded flags
#' @export
#' 
#' @details
#' Convert numeric flags to text. Those are the flags:
#' "1"  = "date is range",
#' "2"  = "GLORYS coordinate is approximated",
#' "4"  = "Minimum depth differs >5m from true value",
#' "8"  = "Maximum depth differs >5m from true value",
#' "16" = "CoralTempSST coordinate is approximated",
#' "32" = "MUR SST coordinate is approximated",
#' "64" = "OSTIA SST coordinate is approximated"
#'
#' @examples
#' \dontrun{
#' text_flags <- decode_flag(c(1,3,12))
#' }
#'
decode_flag <- function(flag, collapse = TRUE) {

    .decode <- function(flag, collapse) {
        if (!is.na(flag)) {
            flag_text <- c(
                "1"  = "date is range",
                "2"  = "GLORYS coordinate is approximated",
                "4"  = "Minimum depth differs >5m from true value",
                "8"  = "Maximum depth differs >5m from true value",
                "16" = "CoralTempSST coordinate is approximated",
                "32" = "MUR SST coordinate is approximated",
                "64" = "OSTIA SST coordinate is approximated"
            )

            flag_id <- as.numeric(names(flag_text))

            active <- flag_id[bitwAnd(flag, flag_id) > 0]

            if (length(active) == 0) {
                return(NA)
            }

            result <- flag_text[as.character(active)]

            if (collapse) {
                paste(result, collapse = "; ")
            } else {
                result
            }
        } else {
            return(NA)
        }
    }

    unlist(
        lapply(flag, .decode, collapse = collapse),
        use.names = FALSE
    )
}

#' Check if Python virtual environment is active
#'
#' @return an error if not being used
#' @export
#'
#' @examples
#' \dontrun{
#' check_venv()
#' }
#'
check_venv <- function() {
    use_python(".venv/bin/python", required = TRUE)
    if (!grepl("\\.venv/bin/python", py_config()$python)) stop("Check virtual environment")
    return(invisible(NULL))
}

#' Check if the expected number of files was returned for the dataset
#'
#' @param files files vector, as returned by copernicusmarine function
#' @param sel_year the year
#' @param sel_month the month
#' @param start_chr start character of the file
#' @param single or two values indicating the target. If two values, it considers minimum/max
#'
#' @return the files or an error if there is a problem
#' @export
#'
#' @examples
#' \dontrun{
#' check_ifdate(files, sel_year, sel_month, target = c(28, 31))
#' }
#'
check_ifdate <- function(files, sel_year, sel_month, start_chr = "/", target = NULL) {
    files_sb <- unlist(lapply(files, as.character), recursive = T)
    files <- files[grepl(paste0(start_chr, sel_year, sprintf("%02d", sel_month)), files_sb)]
    if (length(files) > 31) {
        stop("More days than maximum in a month")
    } else if (length(files) < 1) {
        stop("No valid files")
    } else {
        if (!is.null(target)) {
            if (length(target) == 1) {
                if (length(files) != target) stop("Number of files different than target")
            } else {
                if (length(files) < target[1]) {
                    stop("Less files than target[1]")
                }
                if (length(files) > target[2]) {
                    stop("More files than target[2]")
                }
            }
        }
        return(files)
    }
}

#' Get MUR dataset list (files and sizes)
#'
#' @return a data.frame
#' @export
#' 
#' @details
#' All MUR files and sizes according to https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/
#'
#' @examples
#' \dontrun{
#' mur_files <- get_mur_ds()
#' }
#'
get_mur_ds <- function() {
    tf <- tempfile(fileext = ".csv")
    download.file("https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/.csv", tf)
    read.csv(tf)
}

#' title
#'
#' @param name what
#'
#' @return what_return
#' @export
#' 
#' @details
#' details
#'
#' @examples
#' \dontrun{
#' example
#' }
#'
safe_download_mur <- function(sel_year, sel_month, destfile, mur_info, ...) {

    url <- c(
        paste0(
            "https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/",
            paste0(
            sel_year, sprintf("%02d", sel_month), "01", sel_year, sprintf("%02d", sel_month),
            lubridate::days_in_month(
                as.Date(paste0(sel_year, sprintf("%02d", sel_month), "01"),
                format = "%Y%m%d"
                )
            )
            ),
            "-GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc"
        )
    )

    target_size <- mur_info[mur_info$Name == basename(url), "Size"]
    if (is.null(target_size) || length(target_size) < 1) target_size <- 500*1000*1000
    td <- try(download.file(url, destfile, ...), silent = T)
    if (inherits(td, "try-error")) {
        fsiz <- 0
    } else {
        fsiz <- file.size(destfile)
    }
    if (fsiz < target_size) {
        message("Download failed. Retrying...")
        max_retry <- 2
        downloaded <- FALSE
        tried <- 0
        while (!downloaded & tried < max_retry) {
            tried <- tried + 1
            td <- try(download.file(url, destfile, ...), silent = T)
            if (!inherits(td, "try-error")) {
                fsiz <- file.size(destfile)
                if (fsiz >= target_size) downloaded <- TRUE
            }
        }
    } else {
        downloaded <- TRUE
    }
    if (!downloaded) message("Problem when downloading file for MUR")
    return(td)
}

#' Set log file status
#'
#' @param log_df log data.frame
#' @param sel_year the year
#' @param sel_month the month
#' @param which_col target column to add the information
#'
#' @return the log_df data.frame, modified
#' @export
#'
#' @examples
#' \dontrun{
#' log_df <- log_df |> set_failed(2010, 3, "status_ostia")
#' }
#'
set_failed <- function(log_df, sel_year, sel_month, which_col) {
    log_df[log_df$year == sel_year & log_df$month == sel_month, which_col] <- "failed_extract"
    return(log_df)
}

#' @rdname set_failed
#' @export
set_success <- function(log_df, sel_year, sel_month, which_col) {
    log_df[log_df$year == sel_year & log_df$month == sel_month, which_col] <- "concluded"
    return(log_df)
}

#' @rdname set_failed
#' @export
set_unavailable <- function(log_df, sel_year, sel_month, which_col) {
    log_df[log_df$year == sel_year & log_df$month == sel_month, which_col] <- "unavailable"
    return(log_df)
}

#' Print message with one additional line
#'
#' @param ... objects to print in the message
#'
#' @return what_return
#' @export
#' 
#' @details
#' This is a wrapper around `cat` that adds a new line ("\n")
#'
#' @examples
#' \dontrun{
#' catn("my message", "with more info")
#' catg("my message", "with more info", "with green background")
#' }
#'
catn <- function(...) {
    cat(..., "\n")
}

#' @rdname catn
#' @export
catg <- function(...) {
    cat("\033[42m", ..., "\033[49m\n")
}
