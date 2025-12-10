load_data_year <- function(sel_year, obis_filt) {

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
    }

    return(obis_sel)
}

filter_data_month <- function(obis_sel, sel_month) {
    
    if (nrow(obis_sel) > 0) {
      obis_sel_month <- obis_sel %>%
        filter(extractedDateYear == sel_year) %>%
        filter(extractedDateMonth == sel_month) %>%
        mutate(temp_ID = seq_len(nrow(.)))
    } else {
      obis_sel_month <- obis_sel
    }

    return(obis_sel_month)
}
