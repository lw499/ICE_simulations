#' Bootstrap Observed Data and Simulate Under All Interventions
#'
#' This internal function bootstraps the observed data (i.e., resamples the observed data set with replacement to
#' construct bootstrap confidence intervals and standard error). Then, the function simulates data
#' using the resampled dataset to fit models and generate values.
#'
#' @param r                       Integer specifying the index of the current iteration of the bootstrap.
#' @param time_points             Number of time points to simulate.
#' @param obs_data                Data table containing the observed data.
#' @param bootseeds               Vector of integers specifying the seeds. One seed is used to initialize each bootstrap iteration.
#' @param outcome_type            Character string specifying the "type" of the outcome. The possible "types" are: \code{"survival"}, \code{"continuous_eof"}, and \code{"binary_eof"}.
#' @param intvars                 Vector of character strings specifying the names of the variables to be intervened
#'                                on in each round of the simulation.
#' @param interventions           List of vectors. Each vector contains a function
#'                                implementing a particular intervention, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the risk ratio. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, and \code{"absorbing"}.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}.
#' @param histvars                Vector of character strings specifying the names of the variables for which history functions
#'                                are to be applied.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}.
#' @param ymodel                  Model statement for the outcome variable.
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
#' @param comprisk                Logical scalar indicating the presence of a competing event.
#' @param compevent_model         Model statement for the competing event variable.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}.
#' @param compevent_name          Character string specifying the name of the competing event variable in \code{obs_data}.
#' @param ranges                  List of vectors. Each vector contains the minimum and
#'                                maximum values of one of the covariates in \code{covnames}.
#' @param yrange                  Vector containing the minimum and maximum values of the
#'                                outcome variable in the observed dataset.
#' @param compevent_range         Vector containing the minimum and maximum values of the
#'                                competing event variable in the observed dataset.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of cores to use in parallel
#'                                simulation.
#' @param max_visits              A vector of one or more values denoting the maximum number of times
#'                                a binary covariate representing a visit process may be missed before
#'                                the individual is censored from the data (in the observed data) or
#'                                a visit is forced (in the simulated data). Multiple values exist in the
#'                                vector when the modeling of more than covariate is attached to a visit
#'                                process. A value of \code{NA} should be provided when there is no visit process.
#' @param hazardratio             Logical scalar indicating whether the hazard ratio should be computed between two interventions.
#' @param intcomp                 List of two numbers indicating the pair interventions to be used for the
#'                                hazard ratio calculation. The default is \code{NA}, resulting in no hazard ratio calculation.
#' @return                        A list with the following components
#' \item{Result}{Matrix containing risks over time under the natural course and under each user-specific intervention.}
#' \item{ResultRatio}{Matrix containing risk ratios over time under the natural course and under each user-specific intervention.}
#' \item{ResultDiff}{Matrix containing risk differences over time under the natural course and under each user-specific intervention.}
#'
#' @import data.table
#' @importFrom dplyr funs n summarise_all group_by
bootstrap_helper <- function(r, time_points, obs_data, bootseeds, outcome_type,
                             intvars, interventions, ref_int,
                             covparams, covnames, covtypes, covfits_custom, basecovs, histvars, histories,
                             ymodel, yrestrictions, compevent_restrictions, restrictions,
                             comprisk, compevent_model,
                             time_name, outcome_name, compevent_name,
                             ranges, yrange, compevent_range, parallel, ncores, max_visits,
                             hazardratio, intcomp){

  data_len <- length(unique(obs_data$newid))
  ids <- as.data.table(sample(1:data_len, data_len, replace = TRUE))
  ids[, 'bid' := 1:data_len]
  colnames(ids) <- c("newid", "bid")
  resample_data <- copy(obs_data)
  setkey(resample_data, "newid")
  resample_data <- resample_data[J(ids), allow.cartesian = TRUE]  # create the new data set names "sample"
  resample_data[, 'newid' := resample_data$bid]
  resample_data[, 'bid' := NULL]

  # Fit models for covariates, outcome, and competing event (if any)
  fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, restrictions = restrictions,
                         time_name = time_name, obs_data = resample_data)
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, resample_data)
  if (comprisk){
    fitD <- pred_fun_D(compevent_model, compevent_restrictions, resample_data)
  } else {
    fitD <- NA
  }

  # Include natural course in interventions
  if (!is.na(interventions)[[1]][[1]]){
    comb_interventions <- c(list(c(natural)), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    # If no intervention specified, use natural course
    comb_interventions <- list(c(natural))
    comb_intvars <- list('none')
  }

  # Simulate data under different interventions
  pools <- lapply(1:length(comb_interventions), FUN = function(i){
    simulate(fitcov = fitcov, fitY = fitY, fitD = fitD,
             yrestrictions = yrestrictions,
             compevent_restrictions = compevent_restrictions,
             restrictions = restrictions,
             outcome_name = outcome_name, compevent_name = compevent_name,
             time_name = time_name,
             intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
             histvars = histvars, histories = histories,
             covparams = covparams, covnames = covnames, covtypes = covtypes,
             basecovs = basecovs,
             comprisk = comprisk, ranges = ranges,
             yrange = yrange, compevent_range = compevent_range,
             outcome_type = outcome_type,
             subseed = bootseeds[r], time_points = time_points,
             obs_data = resample_data, parallel = FALSE, max_visits = max_visits)
  })

  nat_pool <- pools[[1]] # Simulated data under natural course
  pools <- pools[-1]     # Simulated data under various interventions

  if (outcome_type == 'survival'){
    param <- 'poprisk'
  } else if (outcome_type == 'continuous'|| outcome_type == 'continuous_eof'){
    param <- 'Ey'
  } else if (outcome_type == 'binary_eof'){
    param <- 'Py'
  }

  # Initialize result matrix
  if (grepl('eof', outcome_type)){
    result_ratio <- result_diff <- int_result <- rep(NA, length(pools) + 1)
  } else {
    result_ratio <- result_diff <- int_result <-
      matrix(NA, nrow = length(pools) + 1, ncol = time_points)
  }

  # Calculate mean risk at each time point under natural course
  if (grepl('eof', outcome_type)){
    nat_result <- mean(nat_pool[[param]], na.rm = TRUE)
  } else {
    nat_result <- tapply(nat_pool[[param]], nat_pool[[time_name]], FUN = mean)
  }

  if (ref_int == 0){
    ref_result <- nat_result
  } else {
    if (grepl('eof', outcome_type)){
      ref_result <- mean(pools[[ref_int]][[param]], na.rm = TRUE)
    } else {
      # Calculate mean risk at each time point for specified reference intervention
      ref_result <- tapply(pools[[ref_int]][[param]], pools[[ref_int]][[time_name]], FUN = mean)
    }
  }

  # Calculate risk ratio under natural course
  # Calculate risk ratio for remaining interventions
  if (grepl('eof', outcome_type)){
    int_result[1] <- nat_result
    int_result[-1] <- sapply(pools, FUN = function(pool){
      mean(pool[[param]], na.rm = TRUE)
    })
    result_ratio <- int_result / ref_result
    result_diff <- int_result - ref_result
  } else {
    int_result[1, ] <- nat_result
    result_ratio[1, ] <- int_result[1, ] / ref_result
    result_diff[1, ] <- int_result[1, ] - ref_result
    if (!is.na(interventions)[[1]][[1]]){
      for (i in 2:(length(pools) + 1)){
        int_result[i, ] <- tapply(pools[[i - 1]][[param]], pools[[i - 1]][[time_name]], FUN = mean)
        result_ratio[i, ] <- int_result[i, ] / ref_result
        result_diff[i, ] <- int_result[i, ] - ref_result
      }
    }
  }

  # Calculate hazard ratio for specified interventions
  if (hazardratio){
    # Generate dataset containing failure/censor time information for each subject
    # under each intervention
    pools_hr <- lapply(1:length(intcomp), FUN = hr_helper, intcomp = intcomp,
                       time_name = time_name, pools = pools)
    data_hr <- rbindlist(pools_hr)
    names(data_hr)[names(data_hr) == time_name] <- "t0"
    names(data_hr)[names(data_hr) == outcome_name] <- "Y"
    # Factor event variable
    data_hr$event <- factor(data_hr$Ycomp, 0:2, labels=c("censor", "Y", "D"))

    if (comprisk){
      # Calculate subdistribution hazard ratio
      hr_data <- survival::finegray(survival::Surv(t0, event) ~ ., data = data_hr, etype = "Y")
      hr_res <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ regime, data = hr_data)
      hr_res <- exp(hr_res$coefficients)
      #hr_comp <- crr(data_hr$t0, data_hr$Ycomp, data_hr$regime, failcode = 1, cencode = 0)
    }
    else {
      # Calculate cause-specific hazard ratio
      hr_res <- survival::coxph(formula = survival::Surv(t0, Y == "1") ~ regime, data = data_hr)
      hr_res <- exp(hr_res$coefficients)
    }
  } else {
    hr_res <- NA
  }



  final <- list(Result = int_result, ResultRatio = result_ratio, ResultDiff = result_diff, ResultHR = hr_res)
  # final <- list(result_ratio)

  return (final)
}
