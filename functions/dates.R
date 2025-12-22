#' Get converted date
#'
#' @param x date in POSIX format
#'
#' @return converted date
#' @export
#'
get_date <- function(x) {
  as.POSIXct(x / 1000, origin = "1970-01-01")
}

#' Check if date is range
#'
#' @param start start date
#' @param end end date
#'
#' @return flag
#' @export
#'
check_date <- function(start, end) {
  start_d <- lubridate::as_date(paste0(
    lubridate::year(start),
    "-", lubridate::month(start)
  ), format = "%Y-%m")
  
  end_d <- lubridate::as_date(paste0(
    lubridate::year(end),
    "-", lubridate::month(end)
  ), format = "%Y-%m")
  
  difference <- end_d - start_d
  
  flag <- ifelse(difference != 0, 1, 0)
  
  return(flag)
}