# Dataset structure

<!-- 
library(arrow)
ds <- open_dataset("agg2")
ds_s <- schema(ds)
obj <- ds$schema$ToString(truncate = FALSE)
obj <- strsplit(obj, "\n")
obj <- obj[[1]][1:length(colnames(ds))]
obj <- strsplit(obj, ": ")
df <- data.frame(
    column = unlist(lapply(obj, \(x) x[1])),
    class = unlist(lapply(obj, \(x) x[2]))
)
if (file.exists("_columns.csv")) {
    cols_expl <- read.csv("_columns.csv", sep = ";")
    df <- dplyr::left_join(df, cols_expl[,c("column", "explanation")])
    write.csv(df, "_columns.csv", row.names = FALSE)
} else {
    write.csv(df, "_columns.csv", row.names = FALSE)
}
knitr::kable(df)
-->

The dataset is structured as this:

|column                        |class        |explanation                                                                              |
|:-----------------------------|:------------|:----------------------------------------------------------------------------------------|
|_id                           |string       |Globally unique identifier assigned by OBIS                                              |
|dataset_id                    |string       |Internal dataset identifier assigned by OBIS                                             |
|occurrenceID                  |string       |Occurrence ID (DwC)                                                                      |
|datasetID                     |string       |Dataset ID (DwC)                                                                         |
|AphiaID                       |int32        |WoRMS AphiaID                                                                            |
|scientificName                |string       |Scientific name                                                                          |
|species                       |string       |Species (from WoRMS)                                                                     |
|genus                         |string       |Genus (from WoRMS)                                                                       |
|family                        |string       |Family (from WoRMS)                                                                      |
|order                         |string       |Order (from WoRMS)                                                                       |
|class                         |string       |Class (from WoRMS)                                                                       |
|phylum                        |string       |Phylum (from WoRMS)                                                                      |
|kingdom                       |string       |Kingdom (from WoRMS)                                                                     |
|eventDate                     |string       |Event date (DwC)                                                                         |
|date_start                    |double       |Unix timestamp based on eventDate (start)                                                |
|date_mid                      |double       |Unix timestamp based on eventDate (middle)                                               |
|date_end                      |double       |Unix timestamp based on eventDate (end)                                                  |
|decimalLongitude              |double       |Parsed and validated by OBIS                                                             |
|decimalLatitude               |double       |Parsed and validated by OBIS                                                             |
|coordinatePrecision           |string       |Precision of the coordinates                                                             |
|coordinateUncertaintyInMeters |double       |Uncertainty of the coordinates in meters                                                 |
|minimumDepthInMeters          |double       |Maximum depth in meters                                                                  |
|maximumDepthInMeters          |double       |Minimum depth in meters                                                                  |
|absence                       |bool         |If TRUE, is an absence                                                                   |
|flags                         |list<element |OBIS QC flags                                                                            |
|month                         |int32        |Month                                                                                    |
|surfaceTemperature            |double       |GLORYS surface temperature                                                               |
|midTemperature                |double       |GLORYS temperature at the mid-range of the water column                                  |
|deepTemperature               |double       |GLORYS temperature at the maximum range of depth                                         |
|bottomTemperature             |double       |GLORYS bottom temperature                                                                |
|midDepth                      |double       |Depth equivalent to the mid-range                                                        |
|deepDepth                     |double       |Depth equivalent to maximum depth                                                        |
|minimumDepthTemperature       |double       |If “minimumDepthInMeters” is available, temperature at that depth                        |
|maximumDepthTemperature       |double       |If “maximumDepthInMeters” is available, temperature at that depth                        |
|minimumDepthClosestDepth      |double       |If the exact depth of “minimumDepthInMeters” is not available, the closest depth matched |
|maximumDepthClosestDepth      |double       |If the exact depth of “maximumDepthInMeters” is not available, the closest depth matched |
|coraltempSST                  |double       |CoralTemp sea surface temperature                                                        |
|murSST                        |double       |MUR sea surface temperature                                                              |
|ostiaSST                      |double       |OSTIA sea surface temperature                                                            |
|ostiaProduct                  |string       |Which OSTIA product was used. REP = reprocessed, NRT = near real time                    |
|obistherm_flags               |double       |obistherm flags                                                                          |
|h3_7                          |string       |H3 grid cell at resolution 7                                                             |
|geometry                      |binary       |Geometry in WKB format                                                                   |
|year                          |int32        |Year (read from the dataset structure)                                                   | 

The columns `_id` and `dataset_id` enable you to link and join this dataset with the OBIS database.

The `geometry` column is a binary geometry format (see more about the GeoParquet [format here](https://geoparquet.org/)). The `h3_7` column contains the H3 grid code at the resolution 7. You can use this to easily aggregate data. Because Uber's H3 system is hierarchical, you can also aggregate in coarser resolutions. See more about the H3 system [here](https://h3geo.org/) and the resolutions table [here](https://h3geo.org/docs/core-library/restable).

## Temperature columns

* `surfaceTemperature`: this is the GLORYS surface temperature
* `midTemperature`: this is the GLORYS temperature for the mid depth (that is, the mid point between the maximum depth with valid values and the surface depth)
* `deepTemperature`: this is the GLORYS temperature for the maximum depth with valid values
* `bottomTemperature`: this is the GLORYS temperature for the bottom (i.e. the 'Sea water potential temperature at sea floor' variable)

For both `midTemperature` and `deepTemperature` you should look at `midDepth` and `deepDepth` to see which is the depth used.

* `minimumDepthTemperature`: this is the GLORYS temperature for the depth recorded on the column `minimumDepthInMeters`, that is, the minimum depth given by the original data
* `maximumDepthTemperature`: this is the GLORYS temperature for the depth recorded on the column `maximumDepthInMeters`, that is, the maximum depth given by the original data

For both `minimumDepthTemperature` and `maximumDepthTemperature` you should look at `minimumDepthClosestDepth` and `maximumDepthClosestDepth` to see which is the depth that was actually used.

* `coraltempSST`: SST according to the CoralTemp
* `murSST`: SST according to the MUR
* `ostiaSST`: SST according to the OSTIA

For OSTIA, there is an additional column:

* `ostiaProduct`: Which OSTIA product was used. REP = reprocessed, NRT = near real time

## Flags

* NA = no problem identified  
* date is range (i.e. the date_start and date_end are different). Note that the date that is used to retrieve the data is the date_mid/date_year column  
* GLORYS coordinate is approximated (i.e., the target cell had no value - was NA - and we searched for the nearest valid point in the 25 nearest cells)  
* Minimum depth closest value is more than 5 meters different than the true value  
* Maximum depth closest value is more than 5 meters different than the true value  
* CoralTempSST coordinate is approximated   
* MUR SST coordinate is approximated  
* OSTIA SST coordinate is approximated  

Flags can be combined. E.g. for a record, date is range and surfaceTemperature and medianTemperature are approximated. Flags are separated by semicolons: "date is range; surfaceTemperature coordinate is approximated; medianTemperature coordinate is approximated". On `R` you can split the flags with `strsplit(flags, "; ")`.

## Sample of the dataset

<!--
r <- read_parquet("agg2/year=1982/part-0.parquet")
r <- r |> filter(class != "Aves") |> select(-dropped, -geometry) 
knitr::kable(r[1:5,])
rm(r)
-->

|_id                                  |dataset_id                           |occurrenceID                                               |datasetID                                    | AphiaID|scientificName     |species            |genus       |family          |order          |class          |phylum        |kingdom  |eventDate           |   date_start|     date_mid|     date_end| decimalLongitude| decimalLatitude|coordinatePrecision | coordinateUncertaintyInMeters| minimumDepthInMeters| maximumDepthInMeters|absence |flags    | month| surfaceTemperature| midTemperature| deepTemperature| bottomTemperature| midDepth| deepDepth| minimumDepthTemperature| maximumDepthTemperature| minimumDepthClosestDepth| maximumDepthClosestDepth| coraltempSST| murSST| ostiaSST|ostiaProduct | obistherm_flags|h3_7            |
|:------------------------------------|:------------------------------------|:----------------------------------------------------------|:--------------------------------------------|-------:|:------------------|:------------------|:-----------|:---------------|:--------------|:--------------|:-------------|:--------|:-------------------|------------:|------------:|------------:|----------------:|---------------:|:-------------------|-----------------------------:|--------------------:|--------------------:|:-------|:--------|-----:|------------------:|--------------:|---------------:|-----------------:|--------:|---------:|-----------------------:|-----------------------:|------------------------:|------------------------:|------------:|------:|--------:|:------------|---------------:|:---------------|
|e107dd2a-72bd-443b-b84f-1c4993f5164a |0264be1a-9d3f-495d-afcb-ac22718a70ce |2827733.00                                                 |NA                                           |  293541|Anguilla australis |Anguilla australis |Anguilla    |Anguillidae     |Anguilliformes |Teleostei      |Chordata      |Animalia |1982-01-01          | 378691200000| 378691200000| 378691200000|        148.41790|       -38.41510|0.0806              |                          9000|                   NA|                   NA|FALSE   |NO_DEPTH |     1|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|    17.75|REP          |               0|87bf5b0e3ffffff |
|f293c0db-bfea-4301-b272-3b4510e69525 |031dba90-d480-419a-a562-be9a34bc5c49 |QLD-Wildnet-4975932                                        |NA                                           |  345834|Auricularia        |NA                 |Auricularia |Auriculariaceae |Auriculariales |Agaricomycetes |Basidiomycota |Fungi    |1982-01-12          | 379641600000| 379641600000| 379641600000|        153.05939|       -26.04007|NA                  |                          2000|                    0|                  0.0|FALSE   |ON_LAND  |     1|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|    26.29|REP          |              64|87be88765ffffff |
|08222b67-79fd-4074-a14e-70118edc687d |04017518-c9c6-4801-a1c9-8019618c9b0d |Station_154_Date_27JAN1982:14:45:00.000_Gobiosoma_robustum |TPWD_HARC_Texas_Upper_Laguna_Madre_Bag_Seine |  276514|Gobiosoma robustum |Gobiosoma robustum |Gobiosoma   |Gobiidae        |Gobiiformes    |Teleostei      |Chordata      |Animalia |1982-01-27 14:45:00 | 380937600000| 380937600000| 380937600000|        -97.66917|        27.30972|NA                  |                           100|                    0|                  0.5|TRUE    |ON_LAND  |     1|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|    18.69|REP          |              64|8748958b2ffffff |
|f427d2eb-6c0f-4131-9c94-22c13d74b247 |04017518-c9c6-4801-a1c9-8019618c9b0d |Station_39_Date_22JAN1982:14:55:00.000_Gobiosoma_robustum  |TPWD_HARC_Texas_Upper_Laguna_Madre_Bag_Seine |  276514|Gobiosoma robustum |Gobiosoma robustum |Gobiosoma   |Gobiidae        |Gobiiformes    |Teleostei      |Chordata      |Animalia |1982-01-22 14:55:00 | 380505600000| 380505600000| 380505600000|        -97.30417|        27.58417|NA                  |                           100|                    0|                  0.5|TRUE    |ON_LAND  |     1|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|    18.02|REP          |               0|874895572ffffff |
|991ae291-318b-4dee-858e-ae79fbb85de5 |04017518-c9c6-4801-a1c9-8019618c9b0d |Station_42_Date_22JAN1982:15:40:00.000_Gobiosoma_robustum  |TPWD_HARC_Texas_Upper_Laguna_Madre_Bag_Seine |  276514|Gobiosoma robustum |Gobiosoma robustum |Gobiosoma   |Gobiidae        |Gobiiformes    |Teleostei      |Chordata      |Animalia |1982-01-22 15:40:00 | 380505600000| 380505600000| 380505600000|        -97.26389|        27.58778|NA                  |                           100|                    0|                  0.3|TRUE    |ON_LAND  |     1|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|    17.96|REP          |               0|874895509ffffff |

Note: in the table above, we removed the `geometry` column.
