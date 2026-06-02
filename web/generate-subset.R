library(duckdb)
library(glue)
library(dplyr)
library(ggplot2)

con <- dbConnect(duckdb())

ds_path <- "/Volumes/OBIS2/results_final_2/"

dbSendQuery(con, glue("
create or replace view obis as 
    select *
    from read_parquet('{ds_path}*/*.parquet');
"))

dbSendQuery(con, "install h3 from community; load h3;")


# Example 1: Acanthuridae family
acan_data <- dbGetQuery(con, glue(
"WITH top_species AS (
  SELECT species
  FROM obis
  WHERE family = 'Acanthuridae'
    AND absence IS NOT TRUE
    AND species IS NOT NULL
    AND deepTemperature IS NOT NULL
  GROUP BY species
  ORDER BY COUNT(*) DESC
  LIMIT 3
)
SELECT species, AphiaID, year, month, deepTemperature,
h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) as h3_5
FROM obis
WHERE family = 'Acanthuridae'
  AND absence IS NOT TRUE
  AND species IN (SELECT species FROM top_species)
  AND deepTemperature IS NOT NULL;"
))

acan_data |>
    mutate(date = as.Date(paste0(year, "-", month, "-01"))) |>
    group_by(species, date) |>
    summarise(average_temperature = mean(deepTemperature)) |>
    filter(!is.na(average_temperature)) |>
    ggplot() +
        geom_line(
            aes(x = date, y = average_temperature, color = species), alpha = .2
        ) +
        geom_point(
            aes(x = date, y = average_temperature, color = species)
        ) +
        facet_wrap(~species) +
        theme_minimal() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        )

fs::dir_create("web/static/")
arrow::write_parquet(acan_data, "web/static/acanthurus.parquet")


# Example 2: Fiddler crabs
uca <- dbGetQuery(con,
"
select species, AphiaID, year, month, surfaceTemperature as glorysSST, 
coraltempSST, murSST, ostiaSST,
h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) as h3_5
from obis
where species in ('Minuca rapax', 'Leptuca thayeri') and
absence is not true;
")

ocy_proc <- uca |>
    mutate(date = as.Date(paste(year, month, "01", sep = "-"))) |>
    tidyr::pivot_longer(cols = ends_with("SST"),
                        names_to = "SSTsource", values_to = "sst")

ggplot(ocy_proc) +
    geom_point(aes(x = date, y = sst, color = SSTsource), alpha = .5) +
    geom_smooth(aes(x = date, y = sst, group = SSTsource),
                method = "lm", color = "grey30") +
    scale_color_manual(values = c("#11b5ae", "#4046ca", "#f68512", "#de3c82")) +
    xlab(NULL) + ylab("Temperature (°C)") +
    theme_light() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(face = "italic")
    ) +
    facet_grid(SSTsource ~ species)

arrow::write_parquet(uca, "web/static/uca.parquet")

# Example 3 - Gadus morhua
gadus <- dbGetQuery(con,
"
SELECT
  h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 5) AS h3_5,
  MAKE_DATE(year, month, 1) AS date,
  AVG(surfaceTemperature) AS glorysSST,
  AVG(coraltempSST) AS coraltempSST,
  AVG(murSST) AS murSST,
  AVG(ostiaSST) AS ostiaSST,
  AVG(maximumDepthTemperature) AS glorysMaxDepthSST,
  COUNT(*) AS n
FROM obis
WHERE species = 'Gadus morhua'
  AND absence IS NOT TRUE
  AND year IS NOT NULL
  AND month IS NOT NULL
GROUP BY h3_5, date
ORDER BY date, h3_5;
")

gadus |>
    mutate(date_year = lubridate::year(as.Date(date))) |>
    tidyr::pivot_longer(cols = ends_with("SST"),
                        names_to = "SSTsource", values_to = "sst") |>
    group_by(date_year, SSTsource) |>
    filter(!is.na(sst)) |>
    summarise(average = mean(sst), max = max(sst), min = min(sst)) |>
    ggplot() +
        geom_linerange(aes(x = date_year, ymin = min, ymax = max, color = SSTsource)) +
        geom_point(
            aes(x = date_year, y = average, color = SSTsource), alpha = .2
        ) +
        facet_wrap(~ SSTsource) +
        theme_minimal() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "none",
            strip.text.x = element_text(face = "italic")
        ) +
        labs(x = NULL, y = "Temperature")

arrow::write_parquet(gadus, "web/static/gadus.parquet")
