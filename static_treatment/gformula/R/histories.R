#' History functions
#'
#' These functions create new columns in an input data table for covariate histories. Users must specify which covariates are to be used in the history functions. The data table is modified in place.
#'
#' \code{lagged} creates new columns for lagged versions of existing
#' variables in the dataset. The user must specify which variables are to be lagged.
#'
#' \code{cumavg} creates new columns for the cumulative average up until
#' time \eqn{t} of existing variables in the dataset.
#'
#' \code{lagavg} creates new columns for the "lagged cumulative average"
#' (cumulative average up until time t, then lagged by one time unit) up until time \eqn{t} of existing
#' variables in the dataset.
#'
#' @param t          Integer specifying the current time index.
#' @param time_name  Character string specifying the name of the time variable in \code{pool}.
#' @param pool       Data table containing all information prior to time \eqn{t} (\eqn{t} noninclusive).
#' @param histvars   Vector of character strings specifying the names of the variables for which history functions are to be applied.
#' @param id         Character string specifying the name of the ID variable in \code{pool}.
#' @param max_visits This argument is not applicable in these functions.
#' @import data.table
#' @export
lagged <- function(pool, histvars, time_name, t, id, max_visits){
  if (t == 0){
    # At first time point, set all lagged variables equal to 0
    lapply(histvars, FUN = function(histvar){
      classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
      myclass <- paste('as.', classtmp, sep = "")
      pool[pool[[time_name]] == t, (paste("lag_", histvar, sep = "")) := get(myclass)(0)]
    })
  } else {
    # Create columns for lagged variables by setting equal to the actual variable's value
    # at t-1
    current_ids <- unique(pool[pool[[time_name]] == t]$id)
    lapply(histvars, FUN = function(histvar){
      pool[pool[[time_name]] == t, (paste("lag_", histvar, sep = "")) :=
             pool[pool[[time_name]] == t - 1 & id %in% current_ids][[histvar]]]
    })
  }
}

#' @rdname lagged
#' @export
cumavg <- function(pool, histvars, time_name, t, id, max_visits){
  if (t == 0){
    # At first time point, set all cumulative averages equal to the actual value of the
    # variable
    lapply(histvars, FUN = function(histvar){
      pool[pool[[time_name]] == t, (paste("cumavg_", histvar, sep = "")) :=
             as.double(pool[pool[[time_name]] == t][[histvar]])]
    })
  } else {
    # At subsequent time points, create new column containing calculated cumulative
    # average until that time point
    current_ids <- unique(pool[pool[[time_name]] == t]$id)
    lapply(histvars, FUN = function(histvar){
      pool[pool[[time_name]] == t, (paste("cumavg_", histvar, sep = "")) :=
             as.double(tapply(pool[id %in% current_ids][[histvar]],
                              pool[id %in% current_ids]$id, FUN = sum) / (t + 1))]
    })
  }
}

#' @rdname lagged
#' @export
lagavg <- function(pool, histvars, time_name, t, id, max_visits){
  if (t == 0){
    # At first time point, set all lagged cumulative averages equal to 0
    lapply(histvars, FUN = function(histvar){
      pool[pool[[time_name]] == t, (paste("lag_cumavg_", histvar, sep = "")) := as.double(0)]
    })
  } else if (t == 1){
    # At second time point, set all lagged cumulative averages equal to the actual value
    # of the variable
    current_ids <- unique(pool[pool[[time_name]] == t]$id)
    lapply(histvars, FUN = function(histvar){
      pool[pool[[time_name]] == t, (paste("lag_cumavg_", histvar, sep = "")) :=
             as.double(pool[pool[[time_name]] == t - 1 & id %in% current_ids, ][[histvar]])]
    })
  } else {
    # At subsequent time points, create new column containing calculated lagged cumulative
    # average until that time point
    current_ids <- unique(pool[pool[[time_name]] == t]$id)
    lapply(histvars, FUN = function(histvar){
      classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
      myclass <- paste('as.', classtmp, sep = "")
      pool[pool[[time_name]] == t, (paste("lag_cumavg_", histvar, sep = "")) :=
             as.double(tapply(pool[pool[[time_name]] < t & id %in% current_ids, ][[histvar]],
                              pool[pool[[time_name]] < t & id %in% current_ids]$id, FUN = sum) / t)]
    })
  }
}

#' Create Visit Sum Covariate
#'
#' This function assists in the implementation of a visit process by creating a covariate,
#' \code{visit_sum}, that counts the number of visits in the past \code{max_visits} time points. If this
#' number is greater than 0, then the individual has not missed more than the maximum number
#' of visits.
#'
#' @param t          Integer specifying the current time index.
#' @param time_name  Character string specifying the name of the time variable in \code{pool}.
#' @param pool       Data table containing all information prior to time \eqn{t} (\eqn{t} noninclusive).
#' @param histvars   Vector of character strings specifying the names of the variables for which lagged cumulative averages are to
#'                   be created.
#' @param id         Character string specifying the name of the ID variable in \code{pool}.
#' @param max_visits A vector of one or more values denoting the maximum number of times
#'                   a binary covariate representing a visit process may be missed before
#'                   the individual is censored from the data (in the observed data) or
#'                   a visit is forced (in the simulated data). Multiple values exist in the
#'                   vector when the modeling of more than covariate is attached to a visit
#'                   process.
#' @import data.table
#' @export
visit_sum <- function(pool, histvars, time_name, t, id, max_visits){
  for (num in max_visits){
    if (t < num){
      lapply(histvars, FUN = function(histvar){
        pool[pool[[time_name]] == t, (paste("visit_sum_", num, "_", histvar, sep = "")) := 1]
      })
    } else {
      current_ids <- unique(pool[pool[[time_name]] == t]$id)
      lapply(histvars, FUN = function(histvar){
        pool[pool[[time_name]] == t, (paste("visit_sum_", num, "_", histvar, sep = "")) :=
               tapply(pool[pool[[time_name]] >= max(0, t - num) & id %in% current_ids][[histvar]],
                      pool[pool[[time_name]] >= max(0, t - num) & id %in% current_ids]$id, FUN = sum)]
      })
    }
  }
}

#' Generates Functions of History of Existing Covariates
#'
#' This internal function applies the history functions to create new columns in an input data table containing new variables that are functions
#' of the histories of existing variables in the dataset.
#'
#' @param t          Integer specifying the current time index.
#' @param time_name  Character string specifying the name of the time variable in \code{pool}.
#' @param pool       Data table containing all information prior to time \eqn{t} (\eqn{t} noninclusive).
#' @param histvars   Vector of character strings specifying the names of the variables for which history functions are to be applied.
#' @param histories  Vector of history functions to apply to the variables specified in \code{histvars}.
#' @param id         Character string specifying the name of the ID variable in \code{pool}.
#' @param max_visits A vector of one or more values denoting the maximum number of times
#'                   a binary covariate representing a visit process may be missed before
#'                   the individual is censored from the data (in the observed data) or
#'                   a visit is forced (in the simulated data). Multiple values exist in the
#'                   vector when the modeling of more than covariate is attached to a visit
#'                   process. A value of \code{NA} should be provided when there is no visit process.
#' @return           A data table containing all preexisting information at time \eqn{t}, plus columns containing
#'                   new function-of-history variables.

make_histories <- function(pool, histvars, histories, time_name, t, id, max_visits){
  if (!is.na(histvars[1]) && !is.na(histories[1])){
    lapply(1:length(histories), FUN = function(i){
      if (identical(histories[[i]], visit_sum)){
        histories[[i]](pool, histvars, time_name, t, id, max_visits[i])
      } else {
        histories[[i]](pool, histvars, time_name, t, id, max_visits)
      }
    })
  }
}
