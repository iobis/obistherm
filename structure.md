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
|obistherm_flags               |string       |obistherm flags                                                                          |
|adultFunctionalGroup          |string       |Adult functional group (benthic or pelagic) obtained from WoRMS                          |
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
* Date is range (i.e. the date_start and date_end are different). Note that the date that is used to retrieve the data is the date_mid/date_year column  
* GLORYS coordinate is approximated (i.e., the target cell had no value - was NA - and we searched for the nearest valid point in the 25 nearest cells)  
* Minimum depth closest value is more than 5 meters different than the true value  
* Maximum depth closest value is more than 5 meters different than the true value  
* CoralTempSST coordinate is approximated   
* MUR SST coordinate is approximated  
* OSTIA SST coordinate is approximated  
* Difference between minimumDepthInMeters and maximumDepthInMeters is higher than 100m
* No depth information available for the record

And there are also three flags associated with the `adultFunctionalGroup` field:

* FG source: WoRMS - source from the functional group information (for now, only WoRMS)
* FG inherited from higher taxonomic level - the functional group information is inherited from a higher taxonomic level
* FG original value - the original value for the function group (before conversion to benthic/pelagic)

Flags can be combined. E.g. for a record, date is range and surfaceTemperature and medianTemperature are approximated. Flags are separated by semicolons: "date is range; surfaceTemperature coordinate is approximated; medianTemperature coordinate is approximated". On `R` you can split the flags with `strsplit(flags, "; ")`.

## Sample of the dataset

<!--
r <- read_parquet("/Volumes/OBIS2/results_final/year=1982/part-0.parquet")
r <- r |> filter(class != "Aves") |> select(-geometry) 
knitr::kable(r[1:5,])
clipr::write_clip(knitr::kable(r[1:5,]))
rm(r)
-->

|_id                                  |dataset_id                           |occurrenceID                                       |datasetID                              |scientificName        |species               |genus      |family         |order             |class     |phylum   |kingdom  |eventDate           |   date_start|     date_mid|     date_end| decimalLongitude| decimalLatitude|coordinatePrecision | coordinateUncertaintyInMeters| minimumDepthInMeters| maximumDepthInMeters|absence |flags | month| surfaceTemperature| midTemperature| deepTemperature| bottomTemperature| midDepth| deepDepth| minimumDepthTemperature| maximumDepthTemperature| minimumDepthClosestDepth| maximumDepthClosestDepth| coraltempSST| murSST| ostiaSST|ostiaProduct | year|adultFunctionalGroup |obistherm_flags | AphiaID|h3_7            |
|:------------------------------------|:------------------------------------|:--------------------------------------------------|:--------------------------------------|:---------------------|:---------------------|:----------|:--------------|:-----------------|:---------|:--------|:--------|:-------------------|------------:|------------:|------------:|----------------:|---------------:|:-------------------|-----------------------------:|--------------------:|--------------------:|:-------|:-----|-----:|------------------:|--------------:|---------------:|-----------------:|--------:|---------:|-----------------------:|-----------------------:|------------------------:|------------------------:|------------:|------:|--------:|:------------|----:|:--------------------|:---------------|-------:|:---------------|
|3abe74ad-f625-4377-9fa5-1edd0880439d |d17e4cae-baf5-4f8b-980f-aa6806f070e2 |349009_158879_43_F_1_2.1111_W2A_26_43_-9_1_2_-9_34 |https://marineinfo.org/id/dataset/8239 |Myzopsetta ferruginea |Myzopsetta ferruginea |Myzopsetta |Pleuronectidae |Pleuronectiformes |Teleostei |Chordata |Animalia |1982-10-01T22:55:00 | 402278400000| 402278400000| 402278400000|           -58.15|         44.3833|NA                  |                            NA|                    0|                   29|FALSE   |      |    10|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|     11.9|REP          | 1982|NA                   |NA              | 1571969|871a96941ffffff |
|adab43b9-5774-47b7-8674-b3b55acb5d6d |d17e4cae-baf5-4f8b-980f-aa6806f070e2 |349009_158879_43_F_1_2.1111_W2A_27_43_-9_1_2_-9_34 |https://marineinfo.org/id/dataset/8239 |Myzopsetta ferruginea |Myzopsetta ferruginea |Myzopsetta |Pleuronectidae |Pleuronectiformes |Teleostei |Chordata |Animalia |1982-10-01T22:55:00 | 402278400000| 402278400000| 402278400000|           -58.15|         44.3833|NA                  |                            NA|                    0|                   29|FALSE   |      |    10|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|     11.9|REP          | 1982|NA                   |NA              | 1571969|871a96941ffffff |
|24e8014e-ba59-40fd-8fec-4ca7b860aa3e |d17e4cae-baf5-4f8b-980f-aa6806f070e2 |349009_158879_43_F_1_2.1111_W2A_28_43_-9_1_2_-9_34 |https://marineinfo.org/id/dataset/8239 |Myzopsetta ferruginea |Myzopsetta ferruginea |Myzopsetta |Pleuronectidae |Pleuronectiformes |Teleostei |Chordata |Animalia |1982-10-01T22:55:00 | 402278400000| 402278400000| 402278400000|           -58.15|         44.3833|NA                  |                            NA|                    0|                   29|FALSE   |      |    10|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|     11.9|REP          | 1982|NA                   |NA              | 1571969|871a96941ffffff |
|178927c2-b4aa-4525-804a-6dc19b864d85 |d17e4cae-baf5-4f8b-980f-aa6806f070e2 |349009_158879_43_F_1_2.1111_W2A_29_43_-9_1_2_-9_34 |https://marineinfo.org/id/dataset/8239 |Myzopsetta ferruginea |Myzopsetta ferruginea |Myzopsetta |Pleuronectidae |Pleuronectiformes |Teleostei |Chordata |Animalia |1982-10-01T22:55:00 | 402278400000| 402278400000| 402278400000|           -58.15|         44.3833|NA                  |                            NA|                    0|                   29|FALSE   |      |    10|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|     11.9|REP          | 1982|NA                   |NA              | 1571969|871a96941ffffff |
|499eb2e9-39a4-433c-939a-7100e44bb9e8 |d17e4cae-baf5-4f8b-980f-aa6806f070e2 |349009_158879_43_F_1_2.1111_W2A_30_43_-9_1_2_-9_34 |https://marineinfo.org/id/dataset/8239 |Myzopsetta ferruginea |Myzopsetta ferruginea |Myzopsetta |Pleuronectidae |Pleuronectiformes |Teleostei |Chordata |Animalia |1982-10-01T22:55:00 | 402278400000| 402278400000| 402278400000|           -58.15|         44.3833|NA                  |                            NA|                    0|                   29|FALSE   |      |    10|                 NA|             NA|              NA|                NA|       NA|        NA|                      NA|                      NA|                       NA|                       NA|           NA|     NA|     11.9|REP          | 1982|NA                   |NA              | 1571969|871a96941ffffff |

Note: in the table above, we removed the `geometry` column.
