#!/bin/bash
set -euo pipefail

MAX_RETRIES=${MAX_RETRIES:-3}
RETRY_DELAY=${RETRY_DELAY:-300}

log() { echo "[$(date -u '+%Y-%m-%d %H:%M:%S UTC')] $*"; }

cd /app

# Returns 1 if any month has failed_extract status in the storr log
has_failures() {
    Rscript --vanilla -e "
suppressPackageStartupMessages(library(storr))
tryCatch({
    st <- storr_rds('control_storr')
    if (st\$exists('log')) {
        log_df <- st\$get('log')
        failed <- any(apply(log_df[, 3:6], 1, function(x) any(x == 'failed_extract', na.rm = TRUE)))
        quit(status = ifelse(failed, 1, 0))
    }
}, error = function(e) {})
quit(status = 0)
" 2>/dev/null
    return $?
}

log "Starting obistherm pipeline"

for attempt in $(seq 1 "$MAX_RETRIES"); do
    log "=== Attempt $attempt of $MAX_RETRIES ==="

    if [ "$attempt" -gt 1 ]; then
        log "Cleaning up failed storr entries before retry..."
        Rscript check_files.R 2>&1 || log "check_files.R encountered an error (continuing)"
        log "Waiting ${RETRY_DELAY}s before retry..."
        sleep "$RETRY_DELAY"
    fi

    RCODE=0
    Rscript update_temperatures.R 2>&1 || RCODE=$?

    HAS_FAIL=0
    has_failures || HAS_FAIL=$?

    if [ "$RCODE" -eq 0 ] && [ "$HAS_FAIL" -eq 0 ]; then
        log "Pipeline completed successfully on attempt $attempt"
        exit 0
    fi

    log "Attempt $attempt finished: R exit=${RCODE}, failures_in_log=${HAS_FAIL}"

    if [ "$attempt" -eq "$MAX_RETRIES" ]; then
        log "Max retries reached. Pipeline did not fully complete."
        exit 1
    fi
done
