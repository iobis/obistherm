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

    require(duckdb)

    f <- dirname(obis_source) |>
        file.path("obis_partitioned")

    if (skip) {
        message("Skipping data preparing and using old version at ", f)
        return(f)
    }

    message("Partitioning dataset by year")

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
        FROM read_parquet('/Volumes/OBIS2/data/obis_data/*.parquet')
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
    CREATE VIEW obis AS SELECT * FROM read_parquet('", file.path(source_folder, "*/*.parquet"), "');
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
#' @param n_workers maximum number of workers. If NULL it uses number of cores - 2
#' @param memory_limit Memory limit per worker. If NULL uses default. Should be written as "4GB"
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
start_dask <- function(browse = TRUE, n_workers = NULL, memory_limit = NULL) {
    da <- import("dask")
    dd <- import("dask.distributed")
    os <- import("os")
    os$environ["BOKEH_SESSION_TOKEN_EXPIRATION"] <- "86400"
    if (is.null(n_workers)) {
        n_workers <- max(1L, parallel::detectCores() - 2L)
    }
    on_mac <- Sys.info()[["sysname"]] == "Darwin"
    if (on_mac) {
        # On macOS, process-based workers communicate over TCP and are prone to
        # "Connection reset by peer" (errno 54) crashes when the OS kills a worker
        # under memory pressure. Thread-based workers share memory without TCP,
        # eliminating these crashes. numpy/xarray release the GIL so parallelism
        # is preserved.
        cluster <- dd$LocalCluster(
            n_workers = as.integer(n_workers),
            threads_per_worker = 2L,
            memory_limit = memory_limit,
            processes = FALSE
        )
        client <- dd$Client(cluster)
    } else {
        client <- dd$Client(
            n_workers = as.integer(n_workers),
            threads_per_worker = 2L,
            memory_limit = memory_limit
        )
    }
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

#' Safely download MUR data
#'
#' @param sel_year year
#' @param sel_month month
#' @param destfile destination file
#' @param mur_info data.frame containing details, as returned by [get_mur_ds()]
#' @param ... additional parameter passed to [download.file()]
#'
#' @return downloaded file
#' @export
#'
#' @examples
#' \dontrun{
#' safe_download_mur(2010, 10, "tempmur.nc", mur_info)
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


#' Get CoralTemp dataset list (files and sizes)
#'
#' @return a data.frame
#' @export
#' 
#' @details
#' All CoralTemp files and sizes according to https://coastwatch.pfeg.noaa.gov/erddap/files/NOAA_DHW_monthly/
#'
#' @examples
#' \dontrun{
#' ctemp_files <- get_ctemp_ds()
#' }
#'
get_ctemp_ds <- function() {
    tf <- tempfile(fileext = ".csv")
    download.file("https://coastwatch.pfeg.noaa.gov/erddap/files/NOAA_DHW_monthly/.csv", tf)
    read.csv(tf)
}

#' Safely download CoralTemp data
#'
#' @param sel_year year
#' @param sel_month month
#' @param destfile destination file
#' @param ctemp_info data.frame containing details, as returned by [get_ctemp_ds()]
#' @param ... additional parameter passed to [download.file()]
#'
#' @return downloaded file
#' @export
#'
#' @examples
#' \dontrun{
#' safe_download_ctemp(2010, 10, "tempctemp.nc", ctemp_info)
#' }
#'
safe_download_ctemp <- function(sel_year, sel_month, destfile, ctemp_info, ...) {

    url <- paste0(
        "https://coastwatch.pfeg.noaa.gov/erddap/files/NOAA_DHW_monthly/ct5km_sst_ssta_monthly_v31_",
        sel_year, sprintf("%02d", sel_month), ".nc"
    )

    target_size <- ctemp_info[ctemp_info$Name == basename(url), "Size"]
    if (is.null(target_size) || length(target_size) < 1) target_size <- 18*1000*1000
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
    if (!downloaded) message("Problem when downloading file for CoralTemp")
    return(td)
}

#' Backup files
#'
#' @param file file to backup
#' @param what name of the product
#' @param file file to backup
#' @param fyear year of the product
#' @param fmonth month of the product
#' @param backup_folder the folder to save files
#' @param backup logical, if TRUE do the backup
#'
#' @return the file name
#' @export
#'
#' @examples
#' \dontrun{
#' backup_it("cttemp.nc", "coraltemp", 2020, 1, "backup")
#' }
#'
backup_it <- function(file, what, fyear, fmonth, backup_folder, backup = TRUE) {
  if (backup) {
    fs::file_copy(
        file,
        file.path(backup_folder, paste0(what, "_", fyear, "_", fmonth, ".nc"))
    )
  }
  return(file)
}

#' Bypass download
#'
#' @param what name of the product
#' @param folder folder where files are
#' @param fyear year of the product
#' @param fmonth month of the product
#'
#' @return the file name or NULL if file does not exists
#' @export
#'
#' @examples
#' \dontrun{
#' bypass_try("coraltemp", "backup", 2020, 1)
#' }
#'
bypass_try <- function(what, folder, fyear, fmonth, bypass = FALSE) {
    if (bypass) {
        tf <- file.path(folder, paste0(what, "_", fyear, "_", fmonth, ".nc"))
        if (file.exists(tf)) {
            cat("Bypassing download for file", tf, "\n")
            return(tf)
        } else {
            cat("Bypass: file", tf, "not found, proceeding with download.\n")
            return(NULL)
        }
    } else {
        return(NULL)
    }
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

#' Add invalid occurrence IDs to an RDS file
#'
#' @param invalid_ids a vector of invalid occurrence IDs
#' @param sel_year selected year
#' @param sel_month selected month
#' @param rds_path path to save the RDS file
#'
#' @return what_return
#' @export
#' 
#' @details
#' Some occurrences are not matched with temperature (e.g. fall on land). In those cases
#' we record the occurrence IDs, so that when updating the dataset it skips those IDs.
#' When a dataset is updated on OBIS, any changes on the record will give it a new ID.
#' Thus, using this method will not affect getting temperature for records that were
#' eventually updated in terms of their coordinates.
#' 
#' The RDS file is structured as a list, with each combination sel_year + sel_month
#' containing the invalid IDs. E.g.
#' `rds_results[["year=2010_month=10"]]`
#'
#' @examples
#' \dontrun{
#' add_invalid_ids(void_ids, sel_year, sel_month)
#' }
#'
add_invalid_ids <- function(invalid_ids, sel_year, sel_month, rds_path = "_void-ids.rds") {
    ymcomb <- paste0("year=", sel_year, "_month=", sel_month)

    if (!file.exists(rds_path)) {
        message("RDS file not available, creating and adding...")
        rds_list <- list(invalid_ids)
        names(rds_list) <- ymcomb
        saveRDS(rds_list, file = rds_path)
    } else {
        message("Adding invalid _id's to RDS file")
        existent <- readRDS(rds_path)

        if (is.null(existent[[ymcomb]])) {
            existent[[ymcomb]] <- invalid_ids
        } else {
            existent[[ymcomb]] <- c(
                existent[[ymcomb]], invalid_ids
            )
            existent[[ymcomb]] <- existent[[ymcomb]][!duplicated(existent[[ymcomb]])]
        }

        saveRDS(existent, file = rds_path)
    }

    return(invisible(NULL))
}

#' Fetch AphiaIDs from WoRMS by attribute key ID
#'
#' @param attribute_ids integer vector of WoRMS attribute key IDs to query.
#'   Default is c(4L, 70L) for Functional Group and Zonation respectively.
#' @param sleep seconds to wait between requests
#' @param output_dir directory to save the final CSV
#' @param storr_path path for the storr cache directory used to checkpoint progress
#'
#' @return path to the saved CSV file
#' @export
#'
#' @details
#' Each page is written to a storr cache as it is fetched. If the function is
#' interrupted, re-running it resumes from the last completed page. On successful
#' completion the storr is destroyed and results are saved to a CSV in output_dir.
#'
#' @examples
#' \dontrun{
#' path <- get_worms_attributes_by_key()
#' }
#'
get_worms_attributes_by_key <- function(
  attribute_ids = c(4L, 70L),
  sleep = 0.5,
  output_dir = "data",
  storr_path = "worms_attr_storr"
) {
    base_url <- "https://www.marinespecies.org/rest/AphiaIDsByAttributeKeyID"

    st <- storr::storr_rds(storr_path)

    null_chr <- function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)

    parse_page <- function(page, attr_id) {
        rows <- list()
        for (item in page) {
            for (attr in item$Attributes) {
                base_row <- data.frame(
                    AphiaID           = as.integer(item$AphiaID),
                    AphiaID_Inherited = suppressWarnings(as.integer(null_chr(attr$AphiaID_Inherited))),
                    measurementTypeID = as.integer(null_chr(attr$measurementTypeID)),
                    measurementType   = null_chr(attr$measurementType),
                    measurementValue  = null_chr(attr$measurementValue),
                    qualitystatus     = null_chr(attr$qualitystatus),
                    stringsAsFactors  = FALSE
                )
                children <- attr$children
                if (length(children) > 0) {
                    for (child in children) {
                        row <- base_row
                        row$child_measurementTypeID <- as.integer(null_chr(child$measurementTypeID))
                        row$child_measurementType   <- null_chr(child$measurementType)
                        row$child_measurementValue  <- null_chr(child$measurementValue)
                        rows[[length(rows) + 1L]] <- row
                    }
                } else {
                    base_row$child_measurementTypeID <- NA_integer_
                    base_row$child_measurementType   <- NA_character_
                    base_row$child_measurementValue  <- NA_character_
                    rows[[length(rows) + 1L]] <- base_row
                }
            }
        }
        if (length(rows) == 0) {
            return(NULL)
        }
        dplyr::bind_rows(rows) |> dplyr::mutate(attributeKeyID = as.integer(attr_id))
    }

    collect_storr_pages <- function(attr_id) {
        all_keys <- st$list()
        page_keys <- all_keys[grepl(paste0("^page_", attr_id, "_"), all_keys)]
        if (length(page_keys) == 0) {
            return(NULL)
        }
        offsets <- as.integer(sub(paste0("^page_", attr_id, "_"), "", page_keys))
        page_keys <- page_keys[order(offsets)]
        dplyr::bind_rows(lapply(page_keys, function(k) st$get(k)))
    }

    fetch_attr <- function(attr_id) {
        done_key     <- paste0("done_",     attr_id)
        progress_key <- paste0("progress_", attr_id)

        if (st$exists(done_key)) {
            message("Attribute ID ", attr_id, " already complete, loading from storr...")
            return(collect_storr_pages(attr_id))
        }

        offset <- if (st$exists(progress_key)) {
            o <- st$get(progress_key)
            message("Resuming attribute ID ", attr_id, " from offset ", o)
            o
        } else {
            1L
        }

        repeat {
            url <- paste0(base_url, "/", attr_id, "?offset=", offset)
            message("Fetching attribute ID ", attr_id, ", offset ", offset)

            resp <- tryCatch(
                httr::GET(url, httr::timeout(60)),
                error = function(e) {
                    warning("Request failed (attr_id=", attr_id, ", offset=", offset, "): ", conditionMessage(e))
                    NULL
                }
            )

            if (is.null(resp)) break

            sc <- httr::status_code(resp)

            if (sc == 204L) {
                st$set(done_key, TRUE)
                break
            }

            if (httr::http_error(resp)) {
                warning("Request failed (attr_id=", attr_id, ", offset=", offset, "): HTTP ", sc)
                break
            }

            page <- tryCatch(
                jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                                   simplifyDataFrame = FALSE),
                error = function(e) {
                    warning("Parse failed (attr_id=", attr_id, ", offset=", offset, "): ", conditionMessage(e))
                    NULL
                }
            )

            if (is.null(page)) break

            if (length(page) == 0) {
                st$set(done_key, TRUE)
                break
            }

            page_rows <- parse_page(page, attr_id)
            if (!is.null(page_rows)) {
                st$set(paste0("page_", attr_id, "_", offset), page_rows)
            }

            offset <- offset + length(page)
            st$set(progress_key, offset)
            Sys.sleep(sleep)
        }

        Sys.sleep(sleep) 
        collect_storr_pages(attr_id)
    }

    results <- dplyr::bind_rows(lapply(attribute_ids, fetch_attr))

    outfile <- file.path(output_dir, paste0("WoRMS_attributes_", Sys.Date(), ".parquet"))
    arrow::write_parquet(results, outfile)
    message("Saved WoRMS attributes to ", outfile)

    all_done <- all(vapply(attribute_ids, function(id) st$exists(paste0("done_", id)), logical(1)))
    if (all_done) {
        st$destroy()
        message("Storr cache cleared.")
    } else {
        incomplete <- attribute_ids[!vapply(attribute_ids, function(id) st$exists(paste0("done_", id)), logical(1))]
        warning("Some attribute IDs did not complete successfully (", paste(incomplete, collapse = ", "),
                "). Storr cache preserved at: ", storr_path)
    }

    outfile
}

#' Build species functional group attributes
#'
#' Assembles a table of functional group (benthos/pelagic) classifications for
#' all species present in the final dataset, drawing from WoRMS attributes and,
#' optionally, FishBase and SeaLifeBase for species not covered by WoRMS.
#'
#' @param final_dataset path to the Arrow/Parquet dataset produced by the
#'   pipeline (must contain `AphiaID` and `species` columns).
#' @param local_worms optional path to a `WoRMS_attributes_<date>.parquet` file.
#'   If `NULL`, the function searches `output_dir` and picks the most recently
#'   modified matching file.
#' @param output_dir directory used both to search for the WoRMS file (when
#'   `local_worms` is `NULL`) and to write the output parquet. Defaults to
#'   `"data"`.
#' @param get_fishbase logical; if `TRUE`, species without a WoRMS functional
#'   group are looked up in FishBase via \pkg{rfishbase}. Defaults to `FALSE`.
#' @param get_sealifebase logical; if `TRUE`, species without a WoRMS functional
#'   group are looked up in SeaLifeBase via \pkg{rfishbase}. Defaults to
#'   `FALSE`.
#' @param force logical; if TRUE it will generate a new table even if one is
#'  already available.
#'
#' @return `NULL` (invisibly). Writes
#'   `species_functional_attr_<Sys.Date()>.parquet` to `output_dir` as a
#'   side effect.  The file contains columns `AphiaID`, `functional_group`,
#'   `functional_group_original`, `is_inherited`, `source`, and `flag`.
#' @export
#'
#' @details
#' Functional groups are resolved in priority order:
#' \enumerate{
#'   \item WoRMS trait data (measurement type 4, adult life stage preferred).
#'   \item FishBase/SeaLifeBase `DemersPelag` field and ecology table (only for
#'     species missing from WoRMS, and only when the respective `get_*`
#'     argument is `TRUE`).
#' }
#' The `flag` column records the data source and, where applicable, whether the
#' functional group was inherited from a higher taxonomic level.
#'
#' @examples
#' \dontrun{
#' build_species_info(
#'     final_dataset = "data/obis_final",
#'     get_fishbase   = TRUE,
#'     get_sealifebase = TRUE
#' )
#' }
#'
build_species_info <- function(
  final_dataset,
  local_worms = NULL,
  output_dir = "data",
  get_fishbase = FALSE,
  get_sealifebase = FALSE,
  force = FALSE
) {
    find_latest <- function(dir, pattern = "WoRMS_attributes_\\d{4}-\\d{2}-\\d{2}\\.parquet$") {
        files <- list.files(
            dir,
            pattern = pattern,
            full.names = TRUE
        )
        if (length(files) == 0) {
            return(NULL)
        }
        files[order(file.mtime(files), decreasing = TRUE)][1]
    }

    existing_data <- find_latest(output_dir, "species_functional_attr_\\d{4}-\\d{2}-\\d{2}\\.parquet$")

    if (!is.null(existing_data) & !force) {
        return(existing_data)
    } else {
        message("Non existent functional attributes table or `force = TRUE`. Generating new table.")
    }

    if (is.null(local_worms)) {
        local_worms <- find_latest(output_dir)
    }

    if (is.null(local_worms) || !file.exists(local_worms)) {
        stop("No WoRMS file found. Build it first with `get_worms_attributes_by_key()`")
    }

    local_worms <- arrow::read_parquet(local_worms)

    ds <- arrow::open_dataset(final_dataset)

    species_attr <- ds |>
        dplyr::select(AphiaID, species) |>
        dplyr::distinct() |>
        dplyr::collect() |>
        dplyr::filter(!is.na(AphiaID)) |>
        dplyr::mutate(AphiaID = as.integer(AphiaID)) |>
        dplyr::group_by(AphiaID) |>
        dplyr::slice(1) |>
        dplyr::ungroup()

    rm(ds)

    worms_filtered <- local_worms |>
        mutate(is_inherited = ifelse(AphiaID != AphiaID_Inherited, TRUE, FALSE)) |>
        mutate(source = "WoRMS") |>
        filter(measurementTypeID == 4) |>
        mutate(child_measurementValue = ifelse(
            is.na(child_measurementValue), "not available", child_measurementValue
        )) |>
        filter(child_measurementValue %in% c("adult", "Adult", "not available")) |>
        mutate(life_stage = tolower(child_measurementValue)) |>
        mutate(functional_group = case_when(
            grepl("bentho", tolower(measurementValue)) ~ "benthos",
            grepl("pelagic", tolower(measurementValue)) ~ "pelagic",
            grepl("plankton", tolower(measurementValue)) ~ "pelagic",
            grepl("euston", tolower(measurementValue)) ~ "pelagic",
            grepl("nekton", tolower(measurementValue)) ~ "pelagic",
            .default = NA
        )) |>
        select(AphiaID, functional_group,
            functional_group_original = measurementValue,
            life_stage, is_inherited, source
        ) |>
        group_by(AphiaID) |>
        arrange(life_stage == "not available", .by_group = TRUE) |>
        slice(1)

    covered <- worms_filtered |>
        filter(AphiaID %in% species_attr$AphiaID)

    non_covered <- species_attr |>
        filter(!AphiaID %in% covered$AphiaID)

    if (nrow(non_covered) > 0 && (get_fishbase || get_sealifebase)) {
        to_lookup <- non_covered |> dplyr::filter(!is.na(species))
        message("Looking up ", nrow(to_lookup), " species not found in WoRMS using FishBase/SeaLifeBase...")

        pelagic_cols <- c("Epipelagic", "Mesopelagic", "Bathypelagic", "Abyssopelagic", "Hadopelagic", "Pelagic")
        benthic_dp <- c("demersal", "bathydemersal", "reef-associated", "benthopelagic")
        pelagic_dp <- c("pelagic", "pelagic-oceanic", "pelagic-neritic")

        lookup_server <- function(spp, server, source_label) {
            sp_df <- rfishbase::species(spp, server = server) |>
                dplyr::select(SpecCode, Species, DemersPelag) |>
                dplyr::distinct()

            eco_df <- rfishbase::ecology(sp_df$Species, server = server)

            if (!is.null(eco_df) && nrow(eco_df) > 0) {
                eco_summary <- eco_df |>
                    dplyr::select(SpecCode, dplyr::any_of(c("Benthic", pelagic_cols))) |>
                    dplyr::group_by(SpecCode) |>
                    dplyr::slice(1) |>
                    dplyr::ungroup()
                sp_df <- dplyr::left_join(sp_df, eco_summary, by = "SpecCode")
            }

            sp_df |>
                dplyr::mutate(
                    functional_group = dplyr::case_when(
                        tolower(DemersPelag) %in% benthic_dp ~ "benthos",
                        tolower(DemersPelag) %in% pelagic_dp ~ "pelagic",
                        !is.na(Benthic) & Benthic == 1 ~ "benthos",
                        dplyr::if_any(dplyr::any_of(pelagic_cols), ~ !is.na(.) & . == 1) ~ "pelagic",
                        .default = NA
                    ),
                    functional_group_original = DemersPelag,
                    is_inherited = FALSE,
                    source = source_label,
                    life_stage = "not available"
                ) |>
                dplyr::filter(!is.na(functional_group)) |>
                dplyr::select(Species, functional_group, functional_group_original, is_inherited, source)
        }

        fb_result <- dplyr::bind_rows(
            if (get_fishbase) lookup_server(to_lookup$species, "fishbase", "FishBase"),
            if (get_sealifebase) lookup_server(to_lookup$species, "sealifebase", "SeaLifeBase")
        ) |>
            dplyr::distinct(Species, .keep_all = TRUE) |>
            dplyr::inner_join(to_lookup, by = c("Species" = "species")) |>
            dplyr::select(AphiaID, functional_group, functional_group_original, is_inherited, source)

        covered <- dplyr::bind_rows(covered, fb_result)
    }

    final_list <- covered |>
        mutate(flag = case_when(
            source == "WoRMS" & is_inherited ~ "FG source: WoRMS; FG inherited from higher taxonomic level",
            source == "WoRMS" & !is_inherited ~ "FG source: WoRMS",
            source == "SeaLifeBase" ~ "FG source: SeaLifeBase",
            source == "FishBase" ~ "FG source: FishBase"
        )) |>
        mutate(flag = ifelse(
            is.na(flag), flag, ifelse(
                source == "WoRMS" & life_stage == "adult",
                paste(flag, "FG life stage: adult", sep = "; "),
                paste(flag, "FG life stage: not available", sep = "; ")
            )
        )) |>
        mutate(flag = ifelse(
            is.na(flag), flag, ifelse(
                functional_group != functional_group_original,
                paste(paste0("FG original value: ", functional_group_original), flag, sep = "; "),
                flag
            )
        ))

    outf <- file.path(output_dir, paste0("species_functional_attr_", Sys.Date(), ".parquet"))
    arrow::write_parquet(final_list, outf)

    return(outf)
}
