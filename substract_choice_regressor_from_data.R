#' Transform data for estimation #Thao: the regressor of the corresponding chosen y will be deducted by others, return the matrix for each choice situation of each deciders less than 1
#' @description Function that transform 'data' by subtracting the regressor of the chosen alternative from the regressor matrix.
#' @param data \code{data}-object
#' @return transformed \code{data}-object \code{"data_tr"}
#' @keywords internal
#'Thao: data <- data_tr$data
substract_choice_regressor_from_data <- function(data) {
  N <- length(data)
  for (n in 1:N) { #Thao: n <- 1
    data_n <- data[[n]]
    T_n <- length(data_n$y)
    for (t in 1:T_n) { #Thao: t <- 1
      Xnt <- data_n$X[[t]]
      ynt <- as.numeric(data_n$y[t])
      Xntj <- Xnt[ynt, ]
      Xnt <- Xnt[-ynt, ,drop=FALSE] - matrix(1, dim(Xnt)[1] - 1, 1) %*% Xntj
      data_n$X[[t]] <- Xnt
    }
    data[[n]] <- data_n
  }
  return(data)
}
