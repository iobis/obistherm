######################### OBIS sea temperature dataset ########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
################################ Data aggregation ##############################

library(arrow)
library(dplyr)
library(reticulate)
source("functions/utils.R")
reticulate::source_python("functions/convert_geoarrow.py")
settings <- yaml::read_yaml("settings.yml", readLines.warn = FALSE)
input_folder <- settings$outfolder_final
output_folder <- settings$outfolder_aggregated
fs::dir_create(output_folder)

species_info <- build_species_info(input_folder) |>
  arrow::read_parquet()

aggregate_data <- function(input_folder, output_folder,
                           h3_resolutions = c(7L),
                           species_attrs) {
  
  cat("Aggregating as Parquet dataset (hive format)\n")

  file_schema <- schema(read_parquet(list.files(input_folder, full.names = T)[1]))
  file_schema$AphiaID <- arrow::int32()
  file_schema$minimumDepthInMeters <- double()
  file_schema$maximumDepthInMeters <- double()
  file_schema$surfaceTemperature <- double()
  file_schema$midTemperature <- double()
  file_schema$deepTemperature <- double()
  file_schema$bottomTemperature <- double()
  file_schema$midDepth <- double()
  file_schema$deepDepth <- double()
  file_schema$minimumDepthTemperature <- double()
  file_schema$maximumDepthTemperature <- double()
  file_schema$minimumDepthClosestDepth <- double()
  file_schema$maximumDepthClosestDepth <- double()
  file_schema$coraltempSST <- double()
  file_schema$murSST <- double()
  file_schema$ostiaSST <- double()
  
  ds <- open_dataset(input_folder,
                     schema = file_schema)
  
  ds |>
    select(-extractedDateStart, -extractedDateMid, -extractedDateEnd, -date_year, -dropped) |>
    rename(year = extractedDateYear, month = extractedDateMonth) |>
    mutate(year = as.integer(year), 
           month = as.integer(month),
           surfaceTemperature = round(surfaceTemperature, 2),
           midTemperature = round(midTemperature, 2),
           deepTemperature = round(deepTemperature, 2),
           bottomTemperature = round(bottomTemperature, 2),
           minimumDepthTemperature = round(minimumDepthTemperature, 2),
           maximumDepthTemperature = round(maximumDepthTemperature, 2),
           coraltempSST = round(coraltempSST, 2),
           murSST = round(murSST, 2),
           ostiaSST = round(ostiaSST, 2)) |>
    group_by(year) |>
    write_dataset(path = output_folder)

  rm(ds)
  
  cat("Converting to GeoArrow and adding H3\n")

  species_attrs <- species_attrs |>
    select(AphiaID, adultFunctionalGroup = functional_group, fg_flag = flag) |>
    filter(grepl("FG life stage: adult", fg_flag)) |>
    mutate(fg_flag = gsub("; FG life stage: adult|FG life stage: adult", "", fg_flag))

  tf <- list.files(output_folder, recursive = T, full.names = T)

  pb <- progress::progress_bar$new(total = length(tf))
  
  con <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbSendQuery(con, "install h3 from community; load h3;")
  duckdb::duckdb_register(con, "species_attrs_db", species_attrs)

  for (id in seq_along(tf)) {
    pb$tick()

    file_path <- tf[id]
    h3_cols <- paste(
      sprintf("h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, %d) AS h3_%d",
              h3_resolutions, h3_resolutions),
      collapse = ",\n        "
    )

    de <- DBI::dbExecute(con, glue::glue("
      COPY (
        WITH decoded AS (
          SELECT
            * EXCLUDE (obistherm_flags),
            CASE
              WHEN obistherm_flags IS NULL THEN NULL
              ELSE nullif(array_to_string(list_filter([
                CASE WHEN (obistherm_flags::INTEGER &  1) > 0 THEN 'Date is range' END,
                CASE WHEN (obistherm_flags::INTEGER &  2) > 0 THEN 'GLORYS coordinate is approximated' END,
                CASE WHEN (obistherm_flags::INTEGER &  4) > 0 THEN 'Minimum depth differs >5m from true value' END,
                CASE WHEN (obistherm_flags::INTEGER &  8) > 0 THEN 'Maximum depth differs >5m from true value' END,
                CASE WHEN (obistherm_flags::INTEGER & 16) > 0 THEN 'CoralTempSST coordinate is approximated' END,
                CASE WHEN (obistherm_flags::INTEGER & 32) > 0 THEN 'MUR SST coordinate is approximated' END,
                CASE WHEN (obistherm_flags::INTEGER & 64) > 0 THEN 'OSTIA SST coordinate is approximated' END
              ], x -> x IS NOT NULL), '; '), '')
            END AS obistherm_flags
          FROM read_parquet('{file_path}')
        ),
        joined AS (
          SELECT d.*, sa.adultFunctionalGroup, sa.fg_flag
          FROM decoded d
          LEFT JOIN species_attrs_db sa ON d.AphiaID::INTEGER = sa.AphiaID::INTEGER
        ),
        merged_flags AS (
          SELECT
            * EXCLUDE (obistherm_flags, fg_flag),
            CASE
              WHEN fg_flag IS NULL THEN obistherm_flags
              WHEN obistherm_flags IS NULL THEN fg_flag
              ELSE obistherm_flags || '; ' || fg_flag
            END AS obistherm_flags
          FROM joined
        )
        SELECT
          * EXCLUDE (obistherm_flags, AphiaID),
          CASE
            WHEN obistherm_flags IS NULL THEN NULL
            WHEN minimumDepthInMeters IS NULL AND maximumDepthInMeters IS NULL
              THEN obistherm_flags || '; ' || 'No depth information available for the record'
            WHEN minimumDepthInMeters IS NULL OR maximumDepthInMeters IS NULL
              THEN obistherm_flags
            WHEN abs(maximumDepthInMeters - minimumDepthInMeters) > 100
              THEN obistherm_flags || '; ' || 'Difference between minimumDepthInMeters and maximumDepthInMeters is higher than 100m'
            ELSE obistherm_flags
          END AS obistherm_flags,
          AphiaID::INTEGER AS AphiaID,
          {h3_cols}
        FROM merged_flags
      ) TO '{file_path}' (FORMAT PARQUET)
    "))

    proc_res <- process_file(tf[id])

    if (!proc_res) stop("Failed processing file ", paste0("\033[47m", tf[id], "\033[49m"))
  }

  DBI::dbDisconnect(con)
  
  cat("Aggregation concluded \n")
    
  return(invisible(NULL))
  
}

# Aggregate
aggregate_data(input_folder, output_folder, species_attrs = species_info)

### END