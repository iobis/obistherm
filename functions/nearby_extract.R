# deprecated, use nearby_from_nc.R
#' Get nearby points based on the XY coordinates
#'
#' @param coords decimal longitude and decimal latitude (in this order)
#' @param tlayer the SpatRaster
#' @param mode queen for "queen mode" or any other value for a 5x5 matrix
#' @param type to average values use "mean"
#'
#' @return nearby cells
#' @export
#'
#' @examples
#' \dontrun{
#' get_nearby(c(-5, 6), my_rast, mode = "other")
#' }
#'
get_nearby <- function(coords, tlayer, mode = "queen", type = "mean", ...) {
  if (nrow(coords) > 0) {
    tcell <- cellFromXY(tlayer, as.data.frame(coords))
    if (mode == "queen") {
      adj_m <- "queen"
    } else {
      adj_m <- matrix(c(rep(1, 12), 0, rep(1, 12)), 5, 5)
    }
    adj <- adjacent(tlayer, cells = tcell, adj_m)
    
    result <- apply(adj, 1, function(x) {
      
      ext_vals <- terra::extract(tlayer, x, ...)
      
      if (type == "mean") {
        mean(ext_vals[,1], na.rm = T)
      } else {
        vals <- na.omit(ext_vals[,1])
        higher_median <- function(x) {
          x <- sort(x)
          n <- length(x)
          if (n %% 2 == 1) {
            return(x[(n + 1) / 2])
          } else {
            return(x[n / 2 + 1])
          }
        }
        higher_median(vals)
      }
    })
    
    return(as.vector(ifelse(is.na(result), NA, result))) # to change nan for NA
  } else {
    return(NULL)
  }
}


# deprecated, use nearby_from_nc.R
#' Get nearby points based on the XY coordinates for multiple layers
#'
#' @param coords decimal longitude and decimal latitude (in this order)
#' @param tlayer the SpatRaster
#' @param mode queen for "queen mode" or any other value for a 5x5 matrix
#' @param type to average values use "mean"
#'
#' @return nearby cells
#' @export
#'
#' @examples
#' \dontrun{
#' get_nearby_mlayer(c(-5, 6), my_rast, mode = "other")
#' }
#'
get_nearby_mlayer <- function(coords, tlayer, mode = "queen", type = "mean", layer) {
  if (nrow(coords) > 0) {
    tcell <- cellFromXY(tlayer, as.data.frame(coords))
    if (mode == "queen") {
      adj_m <- "queen"
    } else {
      adj_m <- matrix(c(rep(1, 12), 0, rep(1, 12)), 5, 5)
    }
    adj <- adjacent(tlayer, cells = tcell, adj_m)
    
    result <- list()
    for (i in 1:nrow(adj)) {
      adj_sel <- adj[i,]
      xy_c <- xyFromCell(tlayer, adj_sel)
      
      ext_vals <- terra::extract(tlayer[[layer[i]]], xy_c)
      
      if (type == "mean") {
        result[[i]] <- mean(ext_vals[,1], na.rm = T)
      } else {
        vals <- na.omit(ext_vals[,1])
        higher_median <- function(x) {
          x <- sort(x)
          n <- length(x)
          if (n %% 2 == 1) {
            return(x[(n + 1) / 2])
          } else {
            return(x[n / 2 + 1])
          }
        }
        result[[i]] <- higher_median(vals)
      }
    }
    
    result <- unlist(result)
    
    return(as.vector(ifelse(is.na(result), NA, result))) # to change nan for NA
  } else {
    return(NULL)
  }
}
