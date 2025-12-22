#' Load data for a year
#'
#' @param sel_year target year
#' @param obis_filt path to the the filtered dataset (by year)
#' @param outfolder if not NULL, an intermediate file will be saved (recommended), 
#' otherwise will keep it in memory.
#'
#' @return data.frame in memory or open Arrow dataset
#' @export
#'
#' @examples
#' \dontrun{
#' load_data_year(2010, "obis_ds", "temp")
#' }
#'
load_data_year <- function(sel_year, obis_filt, outfolder = NULL) {

    obis_sel <- dbGetQuery(obis_filt,
                           paste0("SELECT * FROM obis WHERE date_year = ", sel_year, ";"))
    obis_sel <- tibble::tibble(obis_sel)

    if (nrow(obis_sel) > 0) {
        obis_sel <- obis_sel %>%
            mutate(
                extractedDateStart = get_date(date_start),
                extractedDateMid = get_date(date_mid),
                extractedDateEnd = get_date(date_end)
            ) %>%
            filter(!is.na(extractedDateMid)) %>%
            mutate(
                extractedDateYear = lubridate::year(extractedDateMid),
                extractedDateMonth = lubridate::month(extractedDateMid)
            ) %>%
            mutate(flagDate = check_date(extractedDateStart, extractedDateEnd))

        if (!is.null(outfolder)) {
            outf <- file.path(outfolder, paste0("obisdata_year=", sel_year, ".parquet"))
            arrow::write_parquet(obis_sel, outf)
            obis_sel <- arrow::open_dataset(outf)
        }
    }

    return(obis_sel)
}

#' Filter data for a month
#'
#' @param obis_sel data.frame or Arrow dataset
#' @param sel_month target month
#'
#' @return filtered data
#' @export
#' 
#' @examples
#' \dontrun{
#' filter_data_month(obis_df, 1)
#' }
#'
filter_data_month <- function(obis_sel, sel_month) {
    
    if (nrow(obis_sel) > 0) {
      obis_sel_month <- obis_sel %>%
        filter(extractedDateYear == sel_year) %>%
        filter(extractedDateMonth == sel_month) %>%
        collect() %>%
        mutate(temp_ID = seq_len(nrow(.)))
    } else {
      obis_sel_month <- obis_sel
    }

    return(obis_sel_month)
}
