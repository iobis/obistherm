######################### OBIS sea temperature dataset ########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
############################## Check log data.frame ############################
# Clean up log file (remove failed_extract/NA) for re-running get_temperatures.R

# Settings
library(storr)
settings <- yaml::read_yaml("settings.yml", readLines.warn = FALSE)
outfolder_final <- settings$outfolder_final
filename <- "var=thetao"
st <- storr_rds("control_storr")

# Check log data.frame
log_df <- st$get("log")

which_failed <- apply(log_df, 1, \(x) {
    any("failed_extract" %in% x[3:6])
})
which_failed <- which(which_failed)

# Clean up
cn <- colnames(log_df)[3:6]

for (k in which_failed) {
    sel_log <- log_df[k,]
    which_failed_l <- cn[which(sel_log[,3:6] %in% "failed_extract")]
    if (length(which_failed_l) > 1) {
        log_df[k, 3:7] <- NA
    } else {
        if (which_failed_l == "status_mur") {
            if (sel_log$year == 2002 & sel_log$month < 6) {
                next
            } else {
                log_df[k, 3:7] <- NA
            }
        } else {
            log_df[k, 3:7] <- NA
        }
    }
    st$del(paste0(sel_log$year, sel_log$month))
    outf <- paste0(outfolder_final, "/", filename, "_year=", sel_log$year, "_month=", sel_log$month, ".parquet")
    if (file.exists(outf)) {
        fs::file_delete(outf)
    }
}

na_entries <- which(is.na(log_df$status_general))

for (i in na_entries) {
    st$del(paste0(log_df$year[i], log_df$month[i]))
    outf <- paste0(outfolder_final, "/", filename, "_year=", log_df$year[i], "_month=", log_df$month[i], ".parquet")
    if (file.exists(outf)) {
        fs::file_delete(outf)
    }
}

# Save log file
st$set("log", log_df)

### END
