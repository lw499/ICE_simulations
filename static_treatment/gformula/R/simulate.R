#' Simulate Binary Values
#'
#' This internal function simulates covariate values from a binomial distribution.
#'
#' @param x     Integer specifying the number of observations to be simulated.
#' @param size  Integer specifying the number of trials.
#' @param prob  Numeric vector specifying the probabilities.
#' @return      Numeric vector of simulated covariate values under the binomial distribution.
#'
predict_binomial <- function(x, size, prob){
  return (stats::rbinom(x, size, prob))
}

#' Simulate Normal Values
#'
#' This internal function simulates covariate values from a normal distribution.
#'
#' @param x       Integer specifying the number of observations to be simulated.
#' @param mean    Numeric scalar specifying the mean of the distribution.
#' @param est_sd  Numeric scalar specifying the standard deviation of the distribution.
#' @return        Numeric vector of simulated covariate values under the normal distribution.
#'
predict_normal <- function(x, mean, est_sd = NA){
  return (stats::rnorm(x, mean, est_sd))
}


#' Simulate Truncated Normal Values
#'
#' This internal function simulates covariate values from a normal distribution truncated
#' on one side.
#'
#' @param x         Integer specifying the number of observations to be simulated.
#' @param mean      Numeric scalar specifying the mean of the distribution.
#' @param est_sd    Numeric scalar specifying the standard deviation of the distribution.
#' @param a         Numeric scalar specifying the lower bound of truncation.
#' @param b         Numeric scalar specifying the upper bound of truncation.
#' @return          Numeric vector of simulated covariate values under the truncated normal
#'                  distribution.
#'
predict_trunc_normal <- function(x, mean, est_sd, a, b){
  return (truncnorm::rtruncnorm(x = x, mean = mean, sd = est_sd, a = a, b = b))
}

#' Simulate Counterfactual Outcomes Under Intervention
#'
#' This internal function simulates a new dataset containing covariates, outcome probabilities, competing event
#' probabilities (if any), outcomes, and competing events (if any) based on an observed
#' dataset and a user-specified intervention.
#'
#' @param o                       Integer specifying the index of the current intervention.
#' @param fitcov                  List of model fits for the time-varying covariates.
#' @param fitY                    Model fit for the outcome variable.
#' @param fitD                    Model fit for the competing event variable, if any.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#' @param compevent_restrictions  List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the competing event variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the competing event variable takes on the value in the
#'                                second entry.
#' @param restrictions            List of vectors. Each vector contains as its first entry
#'                                the covariate affected by the restriction; its second entry
#'                                the condition that must be \code{TRUE} for the covariate to be
#'                                modeled; its third entry a function that executes other
#'                                specific actions based on the condition; and its fourth
#'                                entry some value used by the function.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name          Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param intvars                 Vector of character strings specifying the names of the variables to be intervened
#'                                on in each round of the simulation.
#' @param interventions           List of vectors. Each vector contains a function
#'                                implementing a particular intervention, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param histvars                Vector of character strings specifying the names of the variables for which history functions
#'                                are to be applied.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}.
#' @param comprisk                Logical scalar indicating the presence of a competing event.
#' @param ranges                  List of vectors. Each vector contains the minimum and
#'                                maximum values of one of the covariates in \code{covnames}.
#' @param yrange                  Vector containing the minimum and maximum values of the
#'                                outcome variable in the observed dataset.
#' @param compevent_range         Vector containing the minimum and maximum values of the
#'                                competing event variable in the observed dataset.
#' @param outcome_type            Character string specifying the "type" of the outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param subseed                 Integer specifying the seed for this simulation.
#' @param obs_data                Data table containing the observed data.
#' @param time_points             Number of time points to simulate.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param covnames                Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, and \code{"absorbing"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}.
#' @param max_visits              A vector of one or more values denoting the maximum number of times
#'                                a binary covariate representing a visit process may be missed before
#'                                the individual is censored from the data (in the observed data) or
#'                                a visit is forced (in the simulated data). Multiple values exist in the
#'                                vector when the modeling of more than covariate is attached to a visit
#'                                process.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#' @return                        A data table containing simulated data under the specified intervention.
#' @import data.table
simulate <- function(o, fitcov, fitY, fitD,
                     yrestrictions, compevent_restrictions, restrictions,
                     outcome_name, compevent_name, time_name,
                     intvars, interventions, histvars, histories,
                     comprisk, ranges, yrange, compevent_range,
                     outcome_type, subseed, obs_data, time_points, parallel,
                     covnames, covtypes, covparams, covpredict_custom,
                     basecovs, max_visits, ...){
  set.seed(subseed)

  # Mechanism of passing intervention variable and intervention is different for parallel
  # and non-parallel versions
  if (parallel){
    intvar <- intvars[[o]]
    intervention <- interventions[[o]]
  } else {
    intvar <- intvars
    intervention <- interventions
  }

  rmses <- lapply(1:length(fitcov), FUN = rmse_calculate, fits = fitcov, covnames = covnames,
                  covtypes = covtypes, obs_data = obs_data, outcome_name = outcome_name,
                  time_name = time_name, restrictions = restrictions,
                  yrestrictions = yrestrictions, compevent_restrictions = compevent_restrictions)

  # Initialize
  ids_unique <- unique(obs_data$newid)
  data_len <- length(ids_unique)
  restrict_ids <- rep(0, data_len)
  restrict_counts <- rep(list(rep(0, data_len)), length(restrictions))
  if (!is.na(restrictions[[1]][[1]])){
    restrict_covs <- lapply(restrictions, FUN = function(restriction){restriction[[1]]})
  }

  for (t in ((1:time_points) - 1)){
    if (t == 0){
      # Set simulated covariate values at time t = 0 equal to observed covariate values
      if (!is.na(basecovs[[1]])){
        newdf <- obs_data[obs_data[[time_name]] == t, ][, .SD, .SDcols = c(covnames, basecovs)]
      } else {
        newdf <- obs_data[obs_data[[time_name]] == t, ][, .SD, .SDcols = covnames]
      }
      set(newdf, j = 'id', value = ids_unique)
      set(newdf, j = time_name, value = rep(t, data_len))
      if (!is.na(basecovs[[1]])){
        setcolorder(newdf, c('id', time_name, covnames, basecovs))
      } else {
        setcolorder(newdf, c('id', time_name, covnames))
      }
      # Update datatable with specified treatment regime / intervention for this
      # simulation
      intfunc(newdf, pool = newdf, intervention, intvar, time_name, t)
      # Update datatable with new covariates that are functions of history of existing
      # covariates
      make_histories(newdf, histvars = histvars,
                     histories = histories, time_name = time_name, t = t, id = 'id',
                     max_visits = max_visits)
      # Generate outcome probabilities
      if (outcome_type == 'survival'){
        set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        # newdf$Py <- predict(fitY, type = 'response', newdata = newdf)
      } else if (outcome_type == 'continuous'){
        # If outcome is continuous, predict mean
        set(newdf, j = 'Ey', value = stats::predict(fitY, type = 'response', newdata = newdf))
        # newdf$Ey <- predict(fitY, type = 'response', newdata = newdf)
      } else if (outcome_type == 'continuous_eof'){
        if (t < (time_points - 1)){
          # newdf$Ey <- NA
          set(newdf, j = 'Ey', value = as.double(NA))
        } else if (t == (time_points - 1)){
          # newdf$Ey <- predict(fitY, type = 'response', newdata = newdf)
          set(newdf, j = 'Ey', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      } else if (outcome_type == 'binary_eof'){
        if (t < (time_points - 1)){
          # newdf$Py <- NA
          set(newdf, j = 'Py', value = as.double(NA))
        } else if (t == (time_points - 1)){
          # newdf$Py <- predict(fitY, type = 'response', newdata = newdf)
          set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      }
      if (!is.na(yrestrictions[[1]][[1]])){ # Check if there are restrictions on outcome
        # variable simulation
        for (yrestriction in yrestrictions){
          # Set non-modeled outcome variable values equal to user-specified value
          if (outcome_type == 'survival'){
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))]$Py <- yrestriction[2]
            set(newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))], j = 'Py',
                value = yrestriction[2])
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), Py := yrestriction[2]]
          } else if (outcome_type == 'continuous'){
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))]$Ey <- yrestriction[2]
            set(newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))], j = 'Ey',
                value = yrestriction[2])
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), Ey := yrestriction[2]]
          }
        }
      }
      # Simulate outcome variable
      if (outcome_type == 'survival'){
        # newdf$Y <- rbinom(data_len, 1, newdf$Py)
        set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
      }
      # } else if (outcome_type == 'binary_eof'){
      #   if (t < (time_points - 1)){
      #     newdf$Y <- NA
      #   } else if (t == (time_points - 1)){
      #     newdf$Y <- rbinom(data_len, 1, newdf$Py)
      #   }
      # }
      # Set simulated outcome values outside the observed range to the observed min / max
      if (length(newdf[newdf$Y < yrange[1]]$Y) != 0){
        # newdf[newdf$Y < yrange[1]]$Y <- yrange[1]
        set(newdf[newdf$Y < yrange[1]], j = 'Y', value = yrange[1])
      }
      if (length(newdf[newdf$Y > yrange[2]]$Y) != 0){
        # newdf[newdf$Y > yrange[2]]$Y <- yrange[2]
        set(newdf[newdf$Y > yrange[2]], j = 'Y', value = yrange[2])
      }
      if (outcome_type == 'survival')
      {
        if (comprisk){
          # Predict competing event probabilities
          # newdf$Pd <- predict(fitD, type = 'response', newdata = newdf)
          set(newdf, j = 'Pd', value = stats::predict(fitD, type = 'response', newdata = newdf))
          if (!is.na(compevent_restrictions[[1]][[1]])){ # Check if there are restrictions
            # on competing event variable simulation
            for (compevent_restriction in compevent_restrictions){
              # Set non-modeled competing event values equal to user-specified value
              # newdf[!eval(parse(text = compevent_restriction[1]))]$Pd <-
              #   compevent_restriction[2]
              set(newdf[!eval(parse(text = compevent_restriction[1]))], j = 'Pd',
                  value = compevent_restriction[2])
            }
          }
          # Simulate competing event variable
          # newdf$D <- rbinom(data_len, 1, newdf$Pd)
          set(newdf, j = 'D', value = stats::rbinom(data_len, 1, newdf$Pd))
          # Set simulated competing event values outside the observed range to the observed
          # min / max
          if (length(newdf[newdf$D < compevent_range[1]]$D) != 0){
            # newdf[newdf$D < compevent_range[1]]$D <- compevent_range[1]
            set(newdf[newdf$D < compevent_range[1]], j = 'D', value = compevent_range[1])
          }
          if (length(newdf[newdf$D > compevent_range[2]]$D) != 0){
            # newdf[newdf$D > compevent_range[2]]$D <- compevent_range[2]
            set(newdf[newdf$D > compevent_range[2]], j = 'D', value = compevent_range[2])
          }
          # Calculate probability of death by main event rather than competing event at
          # time t
          # newdf$prodp1 <- newdf$Py * (1 - newdf$Pd)
          set(newdf, j = 'prodp1', value = newdf$Py * (1 - newdf$Pd))
        } else {
          # newdf$D <- 0 # No competing event
          set(newdf, j = 'D', value = 0)
          # Calculate probability of death by main event without competing event
          if (outcome_type == 'survival'){
            # newdf$prodp1 <- newdf$Py
            set(newdf, j = 'prodp1', value = newdf$Py)
          }
        }
        set(newdf[newdf$D == 1], j = 'Y', value = NA)
        set(newdf, j = 'prodp0', value = 1 - newdf$Py)
      }
      # newdf[D == 1]$Y <- NA
      # If competing event occurs, outcome cannot also occur because
      # both presumably lead to death
      # Calculate probability of survival or death by competing event at time t
      pool <- newdf # Create datatable containing all simulated values up until time t
    } else {
      # Set initial simulated values at time t to simulated values at time t - 1, to be
      # updated later
      newdf <- pool[pool[[time_name]] == t - 1]
      set(newdf, j = time_name, value = rep(t, data_len))
      if ('categorical time' %in% covtypes){
        newdf[, (paste(time_name, "_f", sep = "")) :=
                factor(newdf[[time_name]],
                       levels = unique(obs_data[obs_data[[time_name]] > 0][[time_name]]))]
      }
      pool <- rbind(newdf, pool)
      make_histories(pool = pool, histvars = histvars, histories = histories,
                     time_name = time_name, t = t, id = 'id', max_visits = max_visits)
      newdf <- pool[pool[[time_name]] == t]
      for (i in 1:length(covnames)){
        if (covtypes[i] == 'binary'){
          set(newdf, j = covnames[i],
              value = predict_binomial(data_len, 1, stats::predict(fitcov[[i]], type = 'response',
                                                            newdata = newdf)))
        } else if (covtypes[i] == 'normal'){
          set(newdf, j = covnames[i],
              value = predict_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                       newdata = newdf),
                                     est_sd = rmses[[i]]))
        } else if (covtypes[i] == 'categorical'){
          set(newdf, j = covnames[i],
              value = stats::predict(fitcov[[i]], type = 'class', newdata = newdf))
        } else if (covtypes[i] == 'zero-inflated normal'){
          set(newdf, j = paste("I_", covnames[i], sep = ""),
              value = predict_binomial(data_len, 1, stats::predict(fitcov[[i]][[1]], type = 'response',
                                                            newdata = newdf)))
          # Where indicator is nonzero, simulate from Gaussian model
          set(newdf, j = covnames[i],
              value = exp(predict_normal(data_len,
                                         stats::predict(fitcov[[i]][[2]], type = 'response',
                                                 newdata = newdf),
                                         est_sd =  rmses[[i]])))
          set(newdf, j = covnames[i],
              value = newdf[[paste("I_", covnames[i], sep = "")]] * newdf[[covnames[i]]])
          # Remove indicator
          set(newdf, j = paste("I_", covnames[i], sep = ""), value = NULL)
        } else if (covtypes[i] == 'bounded normal'){
          if (!is.na(restrictions[[1]][[1]])){
            restrictnames <- lapply(1:length(restrictions), FUN = function(r){
              restrictions[[r]][[1]]})
            # Create list of conditions where covariates are modeled
            conditions <- lapply(1:length(restrictions), FUN = function(r){
              restrictions[[r]][[2]]
            })
            if (covnames[i] %in% restrictnames){
              j <- which(restrictnames %in% covnames[i])
              condition <- ""
              if (length(j) > 1){
                for (k in j){
                  condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
                  condition_var <- sub(">.*$", "", condition_var)
                  condition_var <- sub("=.*$", "", condition_var)
                  condition_var <- sub("!.*$", "", condition_var)
                  condition_var <- sub("%in%.*$", "", condition_var)
                  condition_var <- sub(" ", "", condition_var)
                  if (condition_var %in% names(obs_data)){
                    if (condition[1] == ""){
                      condition <- conditions[j]
                    } else {
                      condition <- paste(condition, conditions[j], sep = "||")
                    }
                  }
                }
              } else {
                condition_var <- sub("<.*$", "", restrictions[[j]][[2]])
                condition_var <- sub(">.*$", "", condition_var)
                condition_var <- sub("=.*$", "", condition_var)
                condition_var <- sub("!.*$", "", condition_var)
                condition_var <- sub("%in%.*$", "", condition_var)
                condition_var <- sub(" ", "", condition_var)
                if (condition_var %in% names(obs_data)){
                  condition <- conditions[j]
                }
              }
              if (condition[1] != ""){
                sub_obs_data <- subset(obs_data, eval(parse(text = condition)))
              } else {
                sub_obs_data <- obs_data
              }
            } else {
              sub_obs_data <- obs_data
            }
          } else {
            sub_obs_data <- obs_data
          }
          set(newdf, j = paste("norm_", covnames[i], sep = ""),
              value = predict_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                       newdata = newdf), est_sd = rmses[[i]]))
          set(newdf, j = covnames[i],
              value = (newdf[[paste("norm_", covnames[i], sep = "")]] *
                         (max(sub_obs_data[[covnames[i]]]) - min(sub_obs_data[[covnames[i]]]))) +
                min(sub_obs_data[[covnames[i]]]))
          set(newdf, j = paste("norm_", covnames[i], sep = ""), value = NULL)
        } else if (covtypes[i] == 'truncated normal'){
          if (covparams$direction[i] == 'left'){
            set(newdf, j = covnames[i],
                value = predict_trunc_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                               newdata = newdf),
                                             est_sd = rmses[[i]], a = covparams$point[i], b = Inf))
          } else if (covparams$direction[i] == 'right'){
            set(newdf, j = covnames[i],
                value = predict_trunc_normal(data_len, stats::predict(fitcov[[i]], type = 'response',
                                                               newdata = newdf),
                                             est_sd = rmses[[i]], a = - Inf, b = covparams$point[i]))
          }
        } else if (covtypes[i] == 'other'){
          if (!is.na(restrictions[[1]][[1]])){
            restrictnames <- lapply(1:length(restrictions), FUN = function(r){
              restrictions[[r]][[1]]})
            # Create list of conditions where covariates are modeled
            conditions <- lapply(1:length(restrictions), FUN = function(r){
              restrictions[[r]][[2]]
            })
            if (covnames[i] %in% restrictnames){
              j <- which(restrictnames %in% covnames[i])
              condition <- ""
              if (length(j) > 1){
                for (k in j){
                  condition_var <- sub("<.*$", "", restrictions[[k]][[2]])
                  condition_var <- sub(">.*$", "", condition_var)
                  condition_var <- sub("=.*$", "", condition_var)
                  condition_var <- sub("!.*$", "", condition_var)
                  condition_var <- sub("%in%.*$", "", condition_var)
                  condition_var <- sub(" ", "", condition_var)
                  if (condition_var %in% names(obs_data)){
                    if (condition[1] == ""){
                      condition <- conditions[j]
                    } else {
                      condition <- paste(condition, conditions[j], sep = "||")
                    }
                  }
                }
              } else {
                condition_var <- sub("<.*$", "", restrictions[[j]][[2]])
                condition_var <- sub(">.*$", "", condition_var)
                condition_var <- sub("=.*$", "", condition_var)
                condition_var <- sub("!.*$", "", condition_var)
                condition_var <- sub("%in%.*$", "", condition_var)
                condition_var <- sub(" ", "", condition_var)
                if (condition_var %in% names(obs_data)){
                  condition <- conditions[j]
                }
              }
            }
          }
          set(newdf, j = covnames[i],
              value = covpredict_custom[[i]](obs_data, newdf, fitcov[[i]],time_name, t,
                                             condition, covnames[i], ...))
        }
        if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
            covtypes[i] == 'truncated normal'){
          classtmp <- class(newdf[[covnames[i]]])
          myclass <- paste('as.', classtmp, sep = "")

          if (length(newdf[newdf[[covnames[i]]] < ranges[[i]][1]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] < ranges[[i]][1], (covnames[i]) :=
                    get(myclass)(ranges[[i]][1])]
          }
          if (length(newdf[newdf[[covnames[i]]] > ranges[[i]][2]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] < ranges[[i]][1], (covnames[i]) :=
                    get(myclass)(ranges[[i]][1])]
          }
        } else if (covtypes[i] == 'zero-inflated normal') {
          classtmp <- class(newdf[[covnames[i]]])
          myclass <- paste('as.', classtmp, sep = "")
          if (length(newdf[newdf[[covnames[i]]] != 0][newdf[newdf[[covnames[i]]] != 0][[covnames[i]]] < ranges[[i]][1]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] != 0][newdf[newdf[[covnames[i]]] != 0][[covnames[i]]] < ranges[[i]][1], (covnames[i]) := get(myclass)(ranges[[i]][1])]
          }
          if (length(newdf[newdf[[covnames[i]]] != 0][newdf[newdf[[covnames[i]]] != 0][[covnames[i]]] > ranges[[i]][2]][[covnames[i]]]) != 0){
            newdf[newdf[[covnames[i]]] != 0][newdf[newdf[[covnames[i]]] != 0][[covnames[i]]] > ranges[[i]][2], (covnames[i]) := get(myclass)(ranges[[i]][2])]
          }
        }
        # Check if there are restrictions on covariate simulation
        if (!is.na(restrictions[[1]][[1]])){
          lapply(1:length(restrictions), FUN = function(r){
            if (restrictions[[r]][[1]] == covnames[i]){
              restrict_ids <- newdf[!eval(parse(text = restrictions[[r]][[2]]))]$id
              if (length(restrict_ids) != 0){
                restrictions[[r]][[3]](newdf, pool[pool[[time_name]] < t], restrictions[[r]], time_name, t)
              }
            }
          })
        }
        pool[pool[[time_name]] == t] <- newdf
        if (covnames[i] %in% histvars){
          make_histories(pool = pool, histvars = c(covnames[i]), histories = histories,
                         time_name = time_name, t = t, id = 'id', max_visits = max_visits)
          newdf <- pool[pool[[time_name]] == t]
        }
      }
      # Update datatable with specified treatment regime / intervention for this
      # simulation
      newdf <- pool[pool[[time_name]] == t]
      intfunc(newdf, pool, intervention, intvar, time_name, t)
      # Update datatable with new covariates that are functions of history of existing
      # covariates
      pool[pool[[time_name]] == t] <- newdf
      make_histories(pool = pool, histvars = histvars,
                     histories = histories, time_name = time_name, t = t, id = 'id',
                     max_visits = max_visits)
      newdf <- pool[pool[[time_name]] == t]

      # Predict outcome probabilities
      if (outcome_type == 'survival'){
        if (comprisk){
          # Predict competing event probabilities
          # newdf$Pd <- predict(fitD, type = 'response', newdata = newdf)
          set(newdf, j = 'Pd', value = stats::predict(fitD, type = 'response', newdata = newdf))
          if (!is.na(compevent_restrictions[[1]][[1]])){ # Check if there are restrictions
            # on competing event variable
            # simulation
            for (compevent_restriction in compevent_restrictions){
              # Set non-modeled competing event values equal to user-specified value
              # newdf[!eval(parse(text = compevent_restriction[1]))]$Pd <-
              #   compevent_restriction[2]
              newdf[!eval(parse(text = compevent_restriction[1])), "Pd" := compevent_restriction[2]]
            }
          }
          # Simulate competing event variable
          # newdf$D <- rbinom(data_len, 1, newdf$Pd)
          set(newdf, j = 'D', value = stats::rbinom(data_len, 1, newdf$Pd))
          # Set simulated competing event values outside the observed range to the observed
          # min / max
          if (length(newdf[newdf$D < compevent_range[1]]$D) != 0){
            # newdf[newdf$D < compevent_range[1]] <- compevent_range[1]
            newdf[newdf$D < compevent_range[1], "D" := compevent_range[1]]
          }
          if (length(newdf[newdf$D > compevent_range[2]]$D) != 0){
            # newdf[newdf$D > compevent_range[2]] <- compevent_range[2]
            newdf[newdf$D > compevent_range[2], "D" := compevent_range[2]]
          }
        } else {
          # newdf$D <- 0
          set(newdf, j = 'D', value = 0)
        }
        # newdf$Py <- predict(fitY, type = 'response', newdata = newdf)
        set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
      } else if (outcome_type == 'continuous'){
        # Predict mean outcome if outcome is continuous
        # newdf$Ey <- (predict(fitY, type = 'response', newdata = newdf) *
        #   (max(obs_data[[outcome_name]]) - min(obs_data[[outcome_name]]))) +
        #   min(obs_data[[outcome_name]])
        set(newdf, j = 'Ey',
            value = stats::predict(fitY, type = 'response', newdata = newdf) *
              (max(obs_data[[outcome_name]]) - min(obs_data[[outcome_name]]))) +
          min(obs_data[[outcome_name]])
        # newdf$Ey <- predict(fitY, type = 'response', newdata = newdf)
      } else if (outcome_type == 'continuous_eof'){
        if (t < (time_points - 1)){
          # newdf$Ey <- NA
          set(newdf, j = 'Ey', value = as.double(NA))
        } else if (t == (time_points - 1)){
          # newdf$Ey <- (predict(fitY, type = 'response', newdata = newdf) *
          #   (max(obs_data[[outcome_name]], na.rm = TRUE) - min(obs_data[[outcome_name]], na.rm = TRUE))) +
          #   min(obs_data[[outcome_name]], na.rm = TRUE)
          set(newdf, j = 'Ey',
              value = (stats::predict(fitY, type = 'response', newdata = newdf) *
                         (max(obs_data[[outcome_name]], na.rm = TRUE) -
                            min(obs_data[[outcome_name]], na.rm = TRUE))) +
                min(obs_data[[outcome_name]], na.rm = TRUE))
        }
      } else if (outcome_type == 'binary_eof'){
        if (t < (time_points - 1)){
          # newdf$Py <- NA
          set(newdf, j = 'Py', value = as.double(NA))
        } else if (t == (time_points - 1)){
          # newdf$Py <- predict(fitY, type = 'response', newdata = newdf)
          set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
        }
      }
      if (!is.na(yrestrictions[[1]][[1]])){ # Check if there are restrictions on outcome
        # variable simulation
        for (yrestriction in yrestrictions){
          # Set non-modeled outcome variable values equal to user-specified value
          if (outcome_type == 'survival'){
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))]$Py <- yrestriction[2]
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Py" := yrestriction[2]]
          } else if (outcome_type == 'continuous'){
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))]$Ey <- yrestriction[2]
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Ey" := yrestriction[2]]
          } else if (outcome_type == 'continuous_eof' && t == (time_points - 1)){
            # newdf[!eval(parse(text = paste("newdf$", yrestriction[1])))]$Ey <- yrestriction[2]
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Ey" := yrestriction[2]]
          } else if (outcome_type == 'binary_eof' && t == (time_points - 1)){
            newdf[!eval(parse(text = paste("newdf$", yrestriction[1]))), "Py" := yrestriction[2]]
          }
        }
      }
      # Calculate probability of survival or death from competing event (if any) at time t
      if (outcome_type == 'survival'){
        # newdf$prodp0 <- 1 - newdf$Py
        set(newdf, j = 'prodp0', value = 1 - newdf$Py)
      }
      # Simulate outcome variable
      if (outcome_type == 'survival'){
        # newdf$Y <- rbinom(data_len, 1, newdf$Py)
        set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
      }
      # } else if (outcome_type == 'binary_eof'){
      #   if (t < (time_points - 1)){
      #     newdf$Y <- NA
      #   } else if (t == (time_points - 1)){
      #     newdf$Y <- rbinom(data_len, 1, newdf$Py)
      #   }
      # }
      # Set simulated outcome values outside the observed range to the observed min / max
      if (length(newdf[newdf$Y < yrange[1]]$Y) != 0){
        # newdf[newdf$Y < yrange[1]] <- yrange[1]
        newdf[newdf$Y < yrange[1], 'Y' := yrange[1]]
      }
      if (length(newdf[newdf$Y > yrange[2]]$Y) != 0){
        # newdf[newdf$Y > yrange[2]] <- yrange[2]
        newdf[newdf$Y > yrange[2], 'Y' := yrange[2]]
      }
      if (outcome_type == 'survival')
        newdf[newdf$D == 1, 'Y' := NA]
      # newdf[D==1]$Y <- NA
      # If competing event occurs, outcome cannot also occur because
      # both presumably lead to death
      # Calculate probability of death from main event at time t
      if (comprisk){
        # newdf$prodp1 <- newdf$Py * tapply(pool$prodp0, pool$id, FUN = prod) *
        #   tapply(1 - pool$Pd, pool$id, FUN = prod) * (1 - newdf$Pd)
        set(newdf, j = 'prodp1',
            value = newdf$Py * tapply(pool[pool[[time_name]] < t]$prodp0,
                                      pool[pool[[time_name]] < t]$id, FUN = prod) *
              tapply(1 - pool[pool[[time_name]] < t]$Pd, pool[pool[[time_name]] < t]$id,
                     FUN = prod) * (1 - newdf$Pd))
      } else if (outcome_type == 'survival'){
        # newdf$prodp1 <- newdf$Py * tapply(pool$prodp0, pool$id, FUN = prod)
        set(newdf, j = 'prodp1', value = newdf$Py * tapply(pool[pool[[time_name]] < t]$prodp0,
                                                           pool[pool[[time_name]] < t]$id,
                                                           FUN = prod))
      }
      # Add simulated data for time t to aggregate simulated data over time
      pool[pool[[time_name]] == t] <- newdf
      # pool <- rbind(pool, newdf)
    }
  }
  colnames(pool)[colnames(pool) == time_name] <- 't0'
  setorder(pool, id, t0)
  colnames(pool)[colnames(pool) == 't0'] <- time_name
  # Calculate probabiity of death from main event at or before time t for each individual
  # at each time point
  if (outcome_type == 'survival'){
    # pool$poprisk <- ave(pool$prodp1, by = pool$id, FUN = cumsum)
    pool[, 'poprisk' := stats::ave(pool$prodp1, by = pool$id, FUN = cumsum)]
    # pool$survival <- ave(pool$prodp0, by = pool$id, FUN = cumprod)
    pool[, 'survival' := stats::ave(pool$prodp0, by = pool$id, FUN = cumprod)]
  }
  pool2 <- copy(pool)
  return (pool2)
}
