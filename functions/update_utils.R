#' Add previous dataset to a DuckDB connection with the OBIS dataset
#'
#' @param connection previous connection opened using [open_db()]
#' @param s3_path path to the S3 bucket where the dataset is stored, without the last "/"
#'
#' @return A DBI connection
#' @export
#' 
#' @details
#' A temporary view of the files is built in memory
#'
#' @examples
#' \dontrun{
#' db_con <- open_db() |> open_done("s3://obis-products/obistherm")
#' }
#'
open_done <- function(connection, s3_path) {
    message("Installing extensions")
    dbSendQuery(connection, "install httpfs; load httpfs;")
    dbSendQuery(connection, "install icu; load icu;")

    message("Reading previous S3 dataset")
    dbSendQuery(connection, glue::glue(
        "create table obistherm as
            select _id, decimalLongitude, decimalLatitude, year, month
            from read_parquet('{s3_path}/*/*.parquet')"
    ))

    return(connection)
}

#' Get records not matched with previously done dataset for year/month
#'
#' @param connection previous connection opened using [open_db()]
#' @param sel_year target year
#' @param sel_month target_month
#'
#' @return data.frame containing records to be processed
#' @export
#'
#' @examples
#' \dontrun{
#' to_do <- previous_done(obis_ds, 2010, 10)
#' }
#'
previous_done <- function(connection, sel_year, sel_month) {
    pd <- dbGetQuery(connection, glue::glue("
        with current_slice as (
            select 
                o.*,
                extract(year from to_timestamp(date_mid / 1000.0)) as year,
                extract(month from to_timestamp(date_mid / 1000.0)) as month
            from obis o
        )

        select o.*
        from current_slice o
        left join obistherm t
            on o._id = t._id
        where
            o.date_year = {sel_year}
            and o.month = {sel_month}
            and (
                t._id is NULL
                or o.decimalLongitude is distinct from t.decimalLongitude
                or o.decimalLatitude is distinct from t.decimalLatitude
                or o.month is distinct from t.month
                or o.date_year is distinct from t.year
            );
        "
    ))
    pd <- pd |>
        mutate(
            extractedDateStart = get_date(date_start),
            extractedDateMid = get_date(date_mid),
            extractedDateEnd = get_date(date_end)
        ) |>
        filter(!is.na(extractedDateMid)) |>
        mutate(
            extractedDateYear = lubridate::year(extractedDateMid),
            extractedDateMonth = lubridate::month(extractedDateMid)
        ) |>
        mutate(flagDate = check_date(extractedDateStart, extractedDateEnd)) |>
        select(-year, -month)
    pd$temp_ID <- seq_len(nrow(pd))
    return(pd)
}


invalid_ids <- function(dataset, rds_path = "_void-ids.rds") {
    if (file.exists(rds_path)) {
        message("Checking invalid _id's")
        existent <- readRDS(rds_path)
        existent$toremove <- 1
        prev <- nrow(dataset)
        dataset <- left_join(dataset, existent)
        dataset <- dataset[is.na(dataset$toremove), -length(dataset)]
        after <- nrow(dataset)
        message(prev-after, " records removed.")
    } else {
        message("RDS file not available, ignoring...")
    }
    return(dataset)
}

# Ok


cleanup_repeated <- function(outfolder_final, sel_year, sel_month, current_ids) {
    message("Cleaning up previous processed files")
    previous_file <- list.files(outfolder_final)
    previous_file <- previous_file[grepl(
        paste0(
            paste0("year=", sel_year, "_month=", sel_month, ".parquet"), "|",
            paste0("year=", sel_year, "_month=", sel_month, "_part")
        ),
        previous_file
    )]
    r <- lapply(previous_file, \(f) {
        pf_ds <- read_parquet(file.path(outfolder_final, f))
        stnr <- nrow(pf_ds)
        pf_ds <- pf_ds[!pf_ds$`_id` %in% current_ids, ]
        stnr <- stnr - nrow(pf_ds)
        write_parquet(pf_ds, file.path(outfolder_final, f))
        if (stnr > 0) message("Removed ", stnr, " records from ", f)
        return(invisible(NULL))
    })
    return(invisible(NULL))
}
