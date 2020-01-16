#' Estimation of Survival Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this function estimates the risk over time under multiple
#' user-specified interventions using the parametric g-formula. See the documentation of (TO DO - link to vignette) for
#' further details concerning the application and implementation of the parametric g-formula.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}. The default is \code{"id"}.
#' @param time_points             Number of time points to simulate.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}. The default is \code{"t0"}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}. The default is \code{"Y"}.
#' @param compevent_name          Character string specifying the name of the competing event variable in \code{obs_data}. The default is \code{"D"}.
#' @param intvars                 Vector of character strings specifying the names of the variables to be intervened
#'                                on in each round of the simulation. The default is \code{NA}.
#' @param interventions           List of vectors. Each vector contains a function
#'                                implementing a particular intervention, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#'                                The default is \code{NA}.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions} and should not include a description for the natural course.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the risk ratio. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, and \code{"absorbing"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                Vector of character strings specifying the names of the variables for which history functions
#'                                are to be applied. The default is \code{NA}.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param compevent_restrictions  List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the competing event variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the competing event variable takes on the value in the
#'                                second entry. The default is \code{NA}.
#' @param restrictions            List of vectors. Each vector contains as its first entry
#'                                the covariate affected by the restriction; its second entry
#'                                the condition that must be \code{TRUE} for the covariate to be
#'                                modeled; its third entry a function that executes other
#'                                specific actions based on the condition; and its fourth
#'                                entry some value used by the function. The default is \code{NA}.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param compevent_model         Model statement for the competing event variable. The default is \code{NA}.
#' @param intcomp                 List of two numbers indicating the pair interventions to be used for the
#'                                hazard ratio calculation. The default is \code{NA}, resulting in no hazard ratio calculation.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of cores to use in parallel
#'                                simulation.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#' @return                        An object of class \code{gformula_survival}. The object is a list with the following components:
#' \item{result}{Results table containing the estimated risk and risk ratio for all interventions (inculding the natural course) at each time point. If bootstrapping was used, the results table includes the bootstrap mean risk ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{hazardratio_val}{Hazard ratio between two interventions (if applicable).}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation under various interventions only for the first and last time points are printed with the \code{\link{print.gformula_survival}} function. To generate graphs comparing the mean estimated covariate values and risks over time and mean observed covariate values and risks over time, use the \code{\link{plot.gformula_survival}} function.
#'
#'
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
#' @examples
#' ## TO ADD AFTER PAPER IS FINALIZED
#' @import data.table
#' @importFrom dplyr funs n summarise_all group_by
#' @export
gformula_survival <- function(obs_data, id = 'id', time_points,
                              time_name = 't0', covnames, covtypes, covparams,
                              covfits_custom = NA, covpredict_custom = NA,
                              histvars = NA, histories = NA, basecovs = NA,
                              outcome_name = 'Y', ymodel,
                              compevent_name = 'D', compevent_model = NA,
                              intvars = NA, interventions = NA,
                              int_descript = NA, ref_int = 0, intcomp = NA,
                              visitprocess = NA, restrictions = NA,
                              yrestrictions = NA, compevent_restrictions = NA,
                              nsimul = NA, sim_data_b = FALSE, seed,
                              nsamples = 0, parallel = FALSE, ncores = NA, ...){

  colnames(obs_data)[colnames(obs_data) == id] <- 'id'
  for (history in histories) {
    for (t in 0:max(obs_data[[time_name]])) {
      obs_data <- history(pool = obs_data, histvars = histvars,
                          time_name = time_name, t = t, id = 'id',
                          max_visits = NA)[[length(histvars)]]
    }
  }
  colnames(obs_data)[colnames(obs_data) == 'id'] <- id

  outcome_type <- 'survival'
  sample_size <- length(unique(obs_data[[id]]))
  hazardratio <- !(length(intcomp) == 1 && is.na(intcomp))
  comprisk <- !(length(compevent_model) == 1 && is.na(compevent_model))

  for (i in 1:length(covnames)){
    if (covtypes[i] == 'absorbing'){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(covnames[i], paste("lag_", covnames[i], "==0", sep = ""),
                               carry_forward, 1)))
      covtypes[i] <- 'binary'
    }
  }


  max_visits <- NA
  if (!is.na(visitprocess[[1]][[1]])){
    for (vp in visitprocess){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(vp[1], paste("visit_sum_", vp[3], "_", vp[1], "!=0", sep = ""),
                               simple_restriction, 1),
                             c(vp[2], paste(vp[1], "==1", sep = ""), carry_forward)))
      if (is.na(max_visits)){
        max_visits <- as.numeric(vp[3])
      } else {
        max_visits <- c(max_visits, as.numeric(vp[3]))
      }
      if (is.na(histories[1])){
        histories <- c(visit_sum)
      } else {
        histories <- c(visit_sum, histories)
      }
    }
  }

  # Create 1-indexed numerical IDs for observed datasets
  ids <- as.data.table(sort(unique(obs_data[[id]])))
  ids[, 'newid' := seq_len(.N)]
  setkeyv(obs_data, id)
  obs_data <- obs_data[J(ids), allow.cartesian = TRUE]

  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(obs_data$newid))
  }

  error_catch(id, nsimul, intvars, interventions, int_descript,
              covnames, covtypes, basecovs,
              histvars, histories, compevent_model, hazardratio, intcomp,
              time_points, outcome_type, time_name, obs_data, parallel, ncores,
              nsamples, sim_data_b)

  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]

  # Determine ranges of observed covariates and outcome
  ranges <- lapply(1:length(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'zero-inflated normal' || covtypes[i] == 'truncated normal') {
      range(obs_data[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data[obs_data[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })
  yrange <- range(obs_data[[outcome_name]])

  # Fit models to covariates and outcome variable
  fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, restrictions = restrictions,
                         time_name = time_name, obs_data = obs_data)
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, obs_data)

  # If competing event exists, fit model for competing event variable
  if (comprisk){
    fitD <- pred_fun_D(compevent_model, compevent_restrictions, obs_data)
    compevent_range <- range(obs_data[[compevent_name]])
  } else {
    fitD <- NA
    compevent_range <- NA
  }

  len <- length(unique(obs_data$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  if (nsimul < len){
    ids <- as.data.table(sort(sample(unique(obs_data$newid), nsimul, replace = TRUE)))
    colnames(ids) <- "newid"
    ids[, 'sid' := seq_len(.N)]
    obs_data <- merge(ids, obs_data, all.x = TRUE, by = "newid")
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  } else if (nsimul > len){
    ids <- as.data.table(sample(unique(obs_data$newid), nsimul, replace = TRUE))
    ids[, 'newid' := 1:nsimul]
    colnames(ids) <- c("newid", "sid")
    setkeyv(obs_data, "newid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE]
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  }

  # Add natural course to list of interventions
  if (!is.na(interventions)[[1]][[1]]){
    comb_interventions <- c(list(c(natural)), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(c(natural))
    comb_intvars <- list('none')
  }

  # Simulate pooled-over-time datasets containing covariates, outcome, and risk for each
  # subject
  if (parallel){
    if (is.na(ncores)){
      ncores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(ncores)
    parallel::clusterCall(cl, library, package = 'data.table', character.only = TRUE)
    parallel::clusterCall(cl, library, package = 'dplyr', character.only = TRUE)
    if ('categorical' %in% covtypes){
      parallel::clusterCall(cl, library, package = 'nnet', character.only = TRUE)
    }
    if ('truncated normal' %in% covtypes){
      parallel::clusterCall(cl, library, package = 'truncreg', character.only = TRUE)
    }

    clusterAssign(cl, "make_histories", make_histories)
    clusterAssign(cl, "rmse_calculate", rmse_calculate)
    clusterAssign(cl, "pred_fun_cov", pred_fun_cov)
    clusterAssign(cl, "fit_glm", fit_glm)
    clusterAssign(cl, "fit_multinomial", fit_multinomial)
    clusterAssign(cl, "fit_trunc_normal", fit_trunc_normal)
    clusterAssign(cl, "pred_fun_Y", pred_fun_Y)
    clusterAssign(cl, "pred_fun_D", pred_fun_D)
    clusterAssign(cl, "natural", natural)
    clusterAssign(cl, "static", static)
    clusterAssign(cl, "threshold", threshold)

    # Note warning suppression when submitting to CRAN, reference 'future' package
    # as doing the same thing for similar issue
    suppressWarnings(parallel::clusterExport(cl, as.vector(utils::lsf.str())))

    pools <- parallel::parLapply(cl, 1:length(comb_interventions), simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = fitD,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 histvars = histvars, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covfits_custom = covfits_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 yrange = yrange, compevent_range = compevent_range,
                                 outcome_type = outcome_type,
                                 subseed = subseed, time_points = time_points,
                                 obs_data = obs_data, parallel = parallel, max_visits = max_visits, ...)

    # parallel::stopCluster(cl)

  } else {
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
               covfits_custom = covfits_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges, yrange = yrange, compevent_range = compevent_range,
               outcome_type = outcome_type,
               subseed = subseed, time_points = time_points,
               obs_data = obs_data, parallel = parallel, max_visits = max_visits, ...)
    })
  }

  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets

  # Initialize results matrices
  result_ratio <- result_diff <- int_result <-
    matrix(NA, nrow = length(pools) + 1, ncol = time_points)

  # Calculate mean risk over all subjects at each time for natural course
  nat_result <- tapply(nat_pool$poprisk, nat_pool[[time_name]], FUN = mean)

  if (ref_int == 0){
    # Set reference intervention to the natural course
    ref_result <- nat_result
  } else {
    # Set reference intervention as specified
    # Calculate mean risk over all subjects at each time for this intervention
    ref_result <- tapply(pools[[ref_int]]$poprisk, pools[[ref_int]][[time_name]], FUN = mean)
  }

  # Compile results
  int_result[1, ] <- nat_result
  result_ratio[1, ] <- int_result[1, ]/ref_result
  result_diff[1, ] <- int_result[1, ] - ref_result
  # Calculate mean risk over all subjects at each time for all interventions other than
  # the natural course
  if (!is.na(interventions)[[1]][[1]]){
    for (i in 2:(length(pools) + 1)){
      int_result[i, ] <- tapply(pools[[i - 1]]$poprisk, pools[[i - 1]][[time_name]],
                                FUN = mean)
      result_ratio[i, ] <- int_result[i, ]/ref_result
      result_diff[i, ] <- int_result[i, ] - ref_result
    }
  }


  # Calculate user specified number of bootstrap risk ratios
  if (nsamples > 0){
    if (parallel){
      parallel::clusterExport(cl, 'simulate')
      final_bs <- parallel::parLapply(cl, 1:nsamples, bootstrap_helper, time_points = time_points,
                                      obs_data = obs_data, bootseeds = bootseeds,
                                      intvars = intvars, interventions = interventions, ref_int = ref_int,
                                      covparams = covparams, covnames = covnames, covtypes = covtypes,
                                      covfits_custom = covfits_custom, basecovs = basecovs, ymodel = ymodel,
                                      histvars = histvars, histories = histories,
                                      comprisk = comprisk, compevent_model = compevent_model,
                                      yrestrictions = yrestrictions,
                                      compevent_restrictions = compevent_restrictions,
                                      restrictions = restrictions, outcome_type = outcome_type,
                                      ranges = ranges, yrange = yrange, compevent_range = compevent_range,
                                      time_name = time_name, outcome_name = outcome_name,
                                      compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                                      max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp, ...)
      parallel::stopCluster(cl)
    } else {
      final_bs <- lapply(1:nsamples, FUN = bootstrap_helper, time_points = time_points,
                         obs_data = obs_data, bootseeds = bootseeds,
                         intvars = intvars, interventions = interventions, ref_int = ref_int,
                         covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, basecovs = basecovs,
                         ymodel = ymodel,
                         histvars = histvars, histories = histories,
                         comprisk = comprisk, compevent_model = compevent_model,
                         yrestrictions = yrestrictions,
                         compevent_restrictions = compevent_restrictions,
                         restrictions = restrictions,
                         outcome_type = outcome_type,
                         ranges = ranges, yrange = yrange, compevent_range = compevent_range,
                         time_name = time_name, outcome_name = outcome_name,
                         compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                         max_visits = max_visits, hazardratio = hazardratio, intcomp = intcomp, ...)
    }
    comb_result <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$Result))
    }))
    comb_RR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultRatio))
    }))
    comb_RD <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultDiff))
    }))
    comb_HR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(m$ResultHR)
    }))
    # Creat dataframe of risk ratios
    #comb <- rbindlist(lapply(final_bs, FUN = function(m){data.table(t(m))}))

    comb_result$t0 <- comb_RR$t0 <- comb_RD$t0 <-
      rep(min(obs_data[[time_name]]):(min(obs_data[[time_name]]) + time_points - 1))
    # Calculate mean and 95% confidence interval of risk ratios
    CI_result <- comb_result %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))
    CI_RR <- comb_RR %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))
    CI_RD <- comb_RD %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))

    hr_mean <- mean(comb_HR$V1)
    hr_sd <- stats::sd(comb_HR$V1)
    hr_res <- c(hr_mean, hr_mean - 1.96*hr_sd, hr_mean + 1.96*hr_sd)
    names(hr_res) <- c('Est. HR', 'HR lower 95% CI', 'HR upper 95% CI')
  } else {
    # If user does not choose to bootstrap, indicate that no confidence interval has
    # been calculated
    CI_result <- CI_RR <- CI_RD <- 'not available'

    # If user does not choose to bootstrap, compute only point estimate of hazard ratio
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
      }
      else {
        # Calculate cause-specific hazard ratio
        hr_res <- survival::coxph(formula = survival::Surv(t0, Y == "1") ~ regime, data = data_hr)
        hr_res <- exp(hr_res$coefficients)
      }
      names(hr_res) <- "Est. HR"
    } else {
      hr_res <- NA
    }
  }

  plot_info <- get_plot_info(outcome_name = outcome_name,
                             compevent_name = compevent_name,
                             time_name = time_name,
                             time_points = time_points,
                             covnames = covnames,
                             covtypes = covtypes,
                             nat_pool = nat_pool,
                             nat_result = nat_result,
                             comprisk = comprisk,
                             outcome_type = outcome_type,
                             obs_data = obs_data)
  obs_results <- plot_info$obs_results

  # Generate results table
  if (!is.na(interventions)[[1]][[1]]){
    resultdf <- lapply(1:time_points, function(i){
      if (nsamples > 0){
        rowdfs <- lapply(1:dim(int_result)[1], function(k){
          data.table(t = min(obs_data[[time_name]]) + i - 1, Intervention = k - 1, Risk = int_result[k, ][i],
                     Risk_bootmean = eval(parse(text = paste("CI_result$V", k, "_mean", sep = "")))[i],
                     Risk_SE = eval(parse(text = paste("CI_result$V", k, "_se", sep = "")))[i],
                     Risk_CI_LL95 = eval(parse(text = paste("CI_result$V", k, "_CI_lower", sep = "")))[i],
                     Risk_CI_UL95 = eval(parse(text = paste("CI_result$V", k, "_CI_upper", sep = "")))[i],
                     RiskRatio = result_ratio[k, ][i],
                     RR_bootmean = eval(parse(text = paste("CI_RR$V", k, "_mean", sep = "")))[i],
                     RR_SE = eval(parse(text = paste("CI_RR$V", k, "_se", sep = "")))[i],
                     RR_CI_LL95 = eval(parse(text = paste("CI_RR$V", k, "_CI_lower", sep = "")))[i],
                     RR_CI_UL95 = eval(parse(text = paste("CI_RR$V", k, "_CI_upper", sep = "")))[i],
                     RiskDiff = result_diff[k, ][i],
                     RD_bootmean = eval(parse(text = paste("CI_RD$V", k, "_mean", sep = "")))[i],
                     RD_SE = eval(parse(text = paste("CI_RD$V", k, "_se", sep = "")))[i],
                     RD_CI_LL95 = eval(parse(text = paste("CI_RD$V", k, "_CI_lower", sep = "")))[i],
                     RD_CI_UL95 = eval(parse(text = paste("CI_RD$V", k, "_CI_upper", sep = "")))[i])

        })
      } else {
        rowdfs <- lapply(1:dim(int_result)[1], function(k){
          data.table(t = i - 1, Intervention = k - 1, Risk = int_result[k, ][i],
                     RiskRatio = result_ratio[k, ][i], RiskDiff = result_diff[k, ][i])
        })
      }
      rowdfs <- rbindlist(rowdfs)
      if (i <= length(unique(obs_data[[time_name]]))){
        rowdfs[, 'obs_risk' := c(obs_results[[2]][i], rep(NA, dim(rowdfs)[1] - 1))]
      } else {
        rowdfs[, 'obs_risk' := rep(NA, dim(rowdfs)[1])]
      }
      return (rowdfs)
    })
    resultdf <- rbindlist(resultdf)
  } else {
    if (nsamples > 0){
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             Risk = t(int_result), Risk_bootmean = CI_result$mean,
                             Risk_SE = CI_result$se, Risk_CI_LL95 = CI_result$CI_lower,
                             Risk_CI_UL95 = CI_result$CI_upper, RiskRatio = t(result_ratio),
                             RR_bootmean = CI_RR$mean, RR_SE = CI_RR$se,
                             RR_CI_LL95 = CI_RR$CI_lower, RR_CI_UL95 = CI_RR$CI_upper,
                             RiskDiff = t(result_diff), RD_bootmean = CI_RD$mean,
                             RD_SE = CI_RD$se,
                             RD_CI_LL95 = CI_RD$CI_lower, RD_CI_UL95 = CI_RD$CI_upper)
    } else {
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             Risk = t(int_result), RiskRatio = t(result_ratio),
                             RiskDiff = t(result_diff))
    }
    resultdf[, 'obs_risk' := obs_results[[2]]]
  }
  if (nsamples > 0){
    colnames(resultdf) <- c("t", "Interv.", "Est. risk", "Bootstrap mean risk",
                            "Risk SE", "Risk lower 95% CI",
                            "Risk upper 95% CI", "Risk ratio",
                            "Bootstrap mean RR", "RR SE",
                            "RR lower 95% CI", "RR upper 95% CI",
                            "Risk difference", "Bootstrap mean RD", "RD SE",
                            "RD lower 95% CI", "RD upper 95% CI",
                            "Obs. risk")
    setcolorder(resultdf, c("t", "Interv.", "Obs. risk", "Est. risk",
                            "Bootstrap mean risk", "Risk SE",
                            "Risk lower 95% CI", "Risk upper 95% CI",
                            "Risk ratio", "Bootstrap mean RR", "RR SE",
                            "RR lower 95% CI", "RR upper 95% CI",
                            "Risk difference", "Bootstrap mean RD", "RD SE",
                            "RD lower 95% CI", "RD upper 95% CI"))
  } else {
    colnames(resultdf) <- c("t", "Interv.", "Est. risk", "Risk ratio",
                            "Risk difference", "Obs. risk")
    setcolorder(resultdf, c("t", "Interv.", "Obs. risk", "Est. risk",
                            "Risk ratio", "Risk difference"))
  }

  fits <- fitcov
  fits[[length(fits) + 1]] <- fitY
  if (!is.na(fitD)[[1]]){
    fits[[length(fits) + 1]] <- fitD
  }

  # Add list of coefficients for covariates, outcome variable, and competing event
  # variable (if any) to results output
  coeffs <- lapply(fits, FUN = function(fit){
    if (length(fit) == 2){
      return (list(stats::coefficients(fit[[1]]), stats::coefficients(fit[[2]])))
    } else {
      return (stats::coefficients(fit))
    }
  })
  if (!is.na(fitD)[[1]]){
    coeffs <-
      stats::setNames(coeffs,
                      c(as.vector(lapply(covnames,
                                         FUN = function(name){paste("fit", name, sep = "")})),
                        "fitY", "fitD"))
  }
  else {
    coeffs <-
      stats::setNames(coeffs,
                      c(as.vector(lapply(covnames,
                                         FUN = function(name){paste("fit", name, sep = "")})),
                        "fitY"))
  }

  rmses <- lapply(1:length(fits), FUN = rmse_calculate, fits = fits, covnames = covnames,
                  covtypes = covtypes, obs_data = obs_data, outcome_name = outcome_name,
                  time_name = time_name, restrictions = restrictions,
                  yrestrictions = yrestrictions, compevent_restrictions = compevent_restrictions)
  if (!is.na(fitD)[[1]]){
    rmses <-
      stats::setNames(rmses,
                      c(as.vector(lapply(covnames,
                                         FUN = function(name){paste("rmse", name, sep = "")})),
                        paste("rmse", outcome_name, sep=""), paste("rmse", compevent_name, sep="")))
  }
  else {
    rmses <-
      stats::setNames(rmses,
                      c(as.vector(lapply(covnames,
                                         FUN = function(name){paste("rmse", name, sep = "")})),
                        paste("rmse", outcome_name, sep="")))
  }

  # Create header
  header <- get_header(int_descript, sample_size, nsimul, nsamples, ref_int)

  if (sim_data_b){
    sim_data <- pools
    if (!is.na(int_descript[1])){
      names(sim_data) <- int_descript
    }
  } else {
    sim_data <- NA
  }

  res <- list(
    result = resultdf,
    coeffs = coeffs,
    rmses = rmses,
    hazardratio_val = hr_res,
    sim_data = sim_data,
    time_name = time_name,
    time_points = time_points,
    covnames = covnames,
    covtypes = covtypes,
    dt_cov_plot = plot_info$dt_cov_plot,
    dt_out_plot = plot_info$dt_out_plot,
    CI_result = CI_result,
    nsamples = nsamples,
    interventions = interventions,
    comprisk = comprisk,
    header = header
  )
  class(res) <- "gformula_survival"
  return (res)
}

#' Estimation of Continuous End-of-Follow-Up Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this function estimates the outcome mean at end-of-follow-up under
#' multiple user-specified interventions using the parametric g-formula. See the documentation of (TO DO - link to vignette) for
#' further details concerning the application and implementation of the parametric g-formula.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}. The default is \code{"id"}.
#' @param time_points             Number of time points to simulate.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}. The default is \code{"t0"}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}. The default is \code{"Y"}.
#' @param intvars                 Vector of character strings specifying the names of the variables to be intervened
#'                                on in each round of the simulation. The default is \code{NA}.
#' @param interventions           List of vectors. Each vector contains a function
#'                                implementing a particular intervention, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#'                                The default is \code{NA}.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions} and should not include a description for the natural course.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the end-of-follow-up mean ratio. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, and \code{"absorbing"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                Vector of character strings specifying the names of the variables for which history functions
#'                                are to be applied. The default is \code{NA}.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param restrictions            List of vectors. Each vector contains as its first entry
#'                                the covariate affected by the restriction; its second entry
#'                                the condition that must be \code{TRUE} for the covariate to be
#'                                modeled; its third entry a function that executes other
#'                                specific actions based on the condition; and its fourth
#'                                entry some value used by the function. The default is \code{NA}.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of cores to use in parallel
#'                                simulation.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#'
#' @return                        An object of class \code{gformula_continuous_eof}. The object is a list with the following components:
#' \item{result}{Results table containing the estimated mean outcome for all interventions (inculding natural course) at the last time point. If bootstrapping was used, the results table includes the bootstrap end-of-follow-up mean ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation under various interventions for the last time point are printed with the \code{\link{print.gformula_continuous_eof}} function. To generate graphs comparing the mean estimated and observed covariate values over time, use the \code{\link{print.gformula_continuous_eof}} function.
#'
#'
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
#' @examples
#' ## TO ADD AFTER PAPER IS FINALIZED
#'
#' @import data.table
#' @importFrom dplyr funs n summarise_all group_by
#' @export
gformula_continuous_eof <- function(obs_data, id = 'id', time_points,
                                    time_name = 't0', covnames, covtypes,
                                    covparams, covfits_custom = NA,
                                    covpredict_custom = NA, histvars = NA,
                                    histories = NA, basecovs = NA,
                                    outcome_name = 'Y', ymodel,
                                    intvars = NA, interventions = NA,
                                    int_descript = NA, ref_int = 0,
                                    visitprocess = NA, restrictions = NA,
                                    yrestrictions = NA, nsimul = NA,
                                    sim_data_b = FALSE,  seed, nsamples = 0,
                                    parallel = FALSE, ncores = NA, ...){

  colnames(obs_data)[colnames(obs_data) == id] <- 'id'
  for (history in histories) {
    for (t in 0:max(obs_data[[time_name]])) {
      obs_data <- history(pool = obs_data, histvars = histvars,
                          time_name = time_name, t = t, id = 'id',
                          max_visits = NA)[[length(histvars)]]
    }
  }
  colnames(obs_data)[colnames(obs_data) == 'id'] <- id

  outcome_type <- 'continuous_eof'
  comprisk = FALSE
  compevent_model = NA
  compevent_name = NA
  compevent_restrictions = NA
  hazardratio <- FALSE
  intcomp <- NA
  sample_size <- length(unique(obs_data[[id]]))

  for (i in 1:length(covnames)){
    if (covtypes[i] == 'absorbing'){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(covnames[i], paste("lag_", covnames[i], "==0", sep = ""),
                               carry_forward, 1)))
      covtypes[i] <- 'binary'
    }
  }

  max_visits <- NA
  if (!is.na(visitprocess[[1]][[1]])){
    for (vp in visitprocess){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(vp[1], paste("visit_sum_", vp[3], "_", vp[1], "!=0", sep = ""),
                               simple_restriction, 1),
                             c(vp[2], paste(vp[1], "==1", sep = ""), carry_forward)))
      if (is.na(max_visits)){
        max_visits <- as.numeric(vp[3])
      } else {
        max_visits <- c(max_visits, as.numeric(vp[3]))
      }
    }
  }

  # Create 1-indexed numerical IDs for observed datasets
  ids <- as.data.table(sort(unique(obs_data[[id]])))
  ids[, 'newid' := seq_len(.N)]
  setkeyv(obs_data, id)
  obs_data <- obs_data[J(ids), allow.cartesian = TRUE]

  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(obs_data$newid))
  }

  error_catch(id, nsimul, intvars, interventions, int_descript,
              covnames, covtypes, basecovs,
              histvars, histories, compevent_model, hazardratio, intcomp,
              time_points, outcome_type, time_name, obs_data, parallel, ncores,
              nsamples, sim_data_b)

  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]

  # Determine ranges of observed covariates and outcome
  ranges <- lapply(1:length(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'zero-inflated normal' || covtypes[i] == 'truncated normal') {
      range(obs_data[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data[obs_data[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })
  yrange <- range(obs_data[[outcome_name]])

  # Fit models to covariates and outcome variable
  fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, restrictions = restrictions,
                         time_name = time_name, obs_data = obs_data)
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, obs_data)

  len <- length(unique(obs_data$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  if (nsimul < len){
    ids <- as.data.table(sort(sample(unique(obs_data$newid), nsimul, replace = TRUE)))
    colnames(ids) <- "newid"
    ids[, 'sid' := seq_len(.N)]
    obs_data <- merge(ids, obs_data, all.x = TRUE, by = "newid")
    obs_data$newid <- obs_data$sid
    obs_data$sid <- NULL
  } else if (nsimul > len){
    ids <- data.frame(sample(unique(obs_data$newid), nsimul, replace = TRUE))
    ids$newid <- 1:nsimul
    colnames(ids) <- c("newid", "sid")
    setkeyv(obs_data, "newid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE]  # create the new data set names "sample"
    obs_data$newid <- obs_data$sid
    obs_data$sid <- NULL
  }

  # Add natural course to list of interventions
  if (!is.na(interventions)[[1]][[1]]){
    comb_interventions <- c(list(c(natural)), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(c(natural))
    comb_intvars <- list('none')
  }

  if (parallel){
    if (is.na(ncores)){
      ncores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(ncores)
    parallel::clusterCall(cl, library, package = 'data.table', character.only = TRUE)
    parallel::clusterCall(cl, library, package = 'dplyr', character.only = TRUE)
    if ('categorical' %in% covtypes){
      parallel::clusterCall(cl, library, package = 'nnet', character.only = TRUE)
    }
    if ('truncated normal' %in% covtypes){
      parallel::clusterCall(cl, library, package = 'truncreg', character.only = TRUE)
    }

    clusterAssign(cl, "make_histories", make_histories)
    clusterAssign(cl, "rmse_calculate", rmse_calculate)
    clusterAssign(cl, "pred_fun_cov", pred_fun_cov)
    clusterAssign(cl, "fit_glm", fit_glm)
    clusterAssign(cl, "fit_multinomial", fit_multinomial)
    clusterAssign(cl, "fit_trunc_normal", fit_trunc_normal)
    clusterAssign(cl, "pred_fun_Y", pred_fun_Y)
    clusterAssign(cl, "pred_fun_D", pred_fun_D)
    clusterAssign(cl, "natural", natural)
    clusterAssign(cl, "static", static)
    clusterAssign(cl, "threshold", threshold)

    # Note warning suppression when submitting to CRAN, reference 'future' package
    # as doing the same thing for similar issue
    suppressWarnings(parallel::clusterExport(cl, as.vector(utils::lsf.str())))

    pools <- parallel::parLapply(cl, 1:length(comb_interventions), simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = NA,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 histvars = histvars, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covfits_custom = covfits_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 yrange = yrange, compevent_range = NA,
                                 outcome_type = outcome_type,
                                 subseed = subseed, time_points = time_points,
                                 obs_data = obs_data, parallel = parallel, ...)

  } else {
    pools <- lapply(1:length(comb_interventions), FUN = function(i){
      simulate(fitcov = fitcov, fitY = fitY, fitD = NA,
               yrestrictions = yrestrictions,
               compevent_restrictions = compevent_restrictions,
               restrictions = restrictions,
               outcome_name = outcome_name, compevent_name = compevent_name,
               time_name = time_name,
               intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
               histvars = histvars, histories = histories,
               covparams = covparams, covnames = covnames, covtypes = covtypes,
               covfits_custom = covfits_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges, yrange = yrange, compevent_range = NA,
               outcome_type = outcome_type,
               subseed = subseed, time_points = time_points,
               obs_data = obs_data, parallel = parallel, ...)
    })
  }

  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets

  # Initialize results matrices
  result_ratio <- result_diff <- int_result <- rep(NA, length(pools) + 1)

  # Calculate mean outcome over all subjects at each time for natural course
  nat_result <- mean(nat_pool$Ey, na.rm = TRUE)

  if (ref_int == 0){
    # Set reference intervention to the natural course
    ref_mean <- nat_result
  } else {
    # Set reference intervention as specified
    # Calculate mean outcome over all subjects at each time for this intervention
    ref_mean <- mean(pools[[ref_int]]$Ey, na.rm = TRUE)
  }

  # Compile results
  int_result[1] <- nat_result
  # Calculate mean risk over all subjects at each time for all interventions other than
  # the natural course
  int_result[-1] <- sapply(pools, FUN = function(pool){mean(pool$Ey, na.rm = TRUE)})
  result_ratio <- int_result / ref_mean
  result_diff <- int_result - ref_mean

  if (nsamples > 0){
    if (parallel){
      parallel::clusterExport(cl, 'simulate')
      final_bs <- parallel::parLapply(cl, 1:nsamples, bootstrap_helper, time_points = time_points,
                                      obs_data = obs_data, bootseeds = bootseeds,
                                      intvars = intvars, interventions = interventions, ref_int = ref_int,
                                      covparams = covparams, covnames = covnames, covtypes = covtypes,
                                      covfits_custom = covfits_custom, basecovs = basecovs, ymodel = ymodel,
                                      histvars = histvars, histories = histories,
                                      comprisk = comprisk, compevent_model = compevent_model,
                                      yrestrictions = yrestrictions,
                                      compevent_restrictions = compevent_restrictions,
                                      restrictions = restrictions, outcome_type = outcome_type,
                                      ranges = ranges, yrange = yrange, compevent_range = NA,
                                      time_name = time_name, outcome_name = outcome_name,
                                      compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                                      max_visits = max_visits, ...)

      parallel::stopCluster(cl)

    } else {
      final_bs <- lapply(1:nsamples, FUN = bootstrap_helper, time_points = time_points,
                         obs_data = obs_data, bootseeds = bootseeds,
                         intvars = intvars, interventions = interventions, ref_int = ref_int,
                         covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, basecovs = basecovs,
                         ymodel = ymodel,
                         histvars = histvars, histories = histories,
                         comprisk = comprisk, compevent_model = compevent_model,
                         yrestrictions = yrestrictions,
                         compevent_restrictions = compevent_restrictions,
                         restrictions = restrictions,
                         outcome_type = outcome_type,
                         ranges = ranges, yrange = yrange, compevent_range = NA,
                         time_name = time_name, outcome_name = outcome_name,
                         compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                         max_visits = max_visits)
    }
    comb_result <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$Result))
    }))
    comb_MR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultRatio))
    }))
    comb_MD <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultDiff))
    }))
    # Creat dataframe of risk ratios
    #comb <- rbindlist(lapply(final_bs, FUN = function(m){data.table(t(m))}))

    comb_result$t0 <- comb_MR$t0 <- comb_MD$t0 <- max(obs_data[[time_name]])
    # Calculate mean and 95% confidence interval of risk ratios
    CI_result <- comb_result %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))
    CI_MR <- comb_MR %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))
    CI_MD <- comb_MD %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))

    # If user does not choose to bootstrap, indicate that no confidence interval has
    # been calculated
  } else {
    CI_result <- CI_MR <- CI_MD <- 'not available'
  }

  plot_info <- get_plot_info(outcome_name = outcome_name,
                             compevent_name = compevent_name,
                             time_name = time_name,
                             time_points = time_points,
                             covnames = covnames,
                             covtypes = covtypes,
                             nat_pool = nat_pool,
                             nat_result = nat_result,
                             comprisk = comprisk,
                             outcome_type = outcome_type,
                             obs_data = obs_data)
  obs_results <- plot_info$obs_results

  # Generate results table
  if (!is.na(interventions)[[1]][[1]]){
    resultdf <- lapply(1:length(int_result), function(k){
      if (nsamples > 0){
        data.table(t = time_points - 1, Intervention = k - 1,
                   EOFMean = int_result[k],
                   EOFMean_bootmean = eval(parse(text = paste("CI_result$V", k, "_mean", sep = ""))),
                   EOFMean_SE = eval(parse(text = paste("CI_result$V", k, "_se", sep = ""))),
                   EOFMean_CI_LL95 = eval(parse(text = paste("CI_result$V", k, "_CI_lower", sep = ""))),
                   EOFMean_CI_UL95 = eval(parse(text = paste("CI_result$V", k, "_CI_upper", sep = ""))),
                   EOFMeanRatio = result_ratio[k],
                   MR_bootmean = eval(parse(text = paste("CI_MR$V", k, "_mean", sep = ""))),
                   MR_SE = eval(parse(text = paste("CI_MR$V", k, "_se", sep = ""))),
                   MR_CI_LL95 = eval(parse(text = paste("CI_MR$V", k, "_CI_lower", sep = ""))),
                   MR_CI_UL95 = eval(parse(text = paste("CI_MR$V", k, "_CI_upper", sep = ""))),
                   EOFMeanDiff = result_diff[k],
                   MD_bootmean = eval(parse(text = paste("CI_MD$V", k, "_mean", sep = ""))),
                   MD_SE = eval(parse(text = paste("CI_MD$V", k, "_se", sep = ""))),
                   MD_CI_LL95 = eval(parse(text = paste("CI_MD$V", k, "_CI_lower", sep = ""))),
                   MD_CI_UL95 = eval(parse(text = paste("CI_MD$V", k, "_CI_upper", sep = ""))))
      } else {
        data.table(t = time_points - 1, Intervention = k - 1, OutcomeEOFMean = int_result[k],
                   EOFMeanRatio = result_ratio[k], EOFMeanDiff = result_diff[k])
      }
    })
    resultdf <- rbindlist(resultdf)
  } else {
    if (nsamples > 0){
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             EOFMean = int_result, EOFMean_bootmean = CI_result$mean,
                             EOFMean_SE = CI_result$se, EOFMean_CI_LL95 = CI_result$CI_lower,
                             EOFMean_CI_UL95 = CI_result$CI_upper, EOFMeanRatio = result_ratio,
                             MR_bootmean = CI_MR$mean, MR_SE = CI_MR$se,
                             MR_CI_LL95 = CI_MR$CI_lower, MR_CI_UL95 = CI_MR$CI_upper,
                             EOFMeanDiff = result_diff, EOFMD_mean = CI_MD$mean,
                             EOFMD_se = CI_MD$se,
                             MD_CI_LL95 = CI_MD$CI_lower, MD_CI_UL95 = CI_MD$CI_upper)
    } else {
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             OutcomeEOFMean = int_result, EOFMeanRatio = result_ratio,
                             EOFMeanDiff = result_diff)
    }
  }

  resultdf$obs_risk <- c(utils::tail(obs_results[[2]], 1), rep(NA, dim(resultdf)[1] - 1))

  if (nsamples > 0){
    colnames(resultdf) <- c("t", "Interv.", "Est. EOF mean",
                            "Bootstrap EOF mean", "EOF mean SE",
                            "EOF mean lower 95% CI", "EOF mean upper 95% CI",
                            "EOF Mean ratio", "Bootstrap EOF MR",
                            "EOF MR SE", "EOF MR lower 95% CI", "EOF MR upper 95% CI",
                            "EOF mean difference", "Bootstrap EOF MD",
                            "EOF MD SE", "EOF MD lower 95% CI", "EOF MD upper 95% CI",
                            "Obs. EOF mean")
    setcolorder(resultdf, c("t", "Interv.", "Obs. EOF mean", "Est. EOF mean",
                            "Bootstrap EOF mean", "EOF mean SE",
                            "EOF mean lower 95% CI", "EOF mean upper 95% CI",
                            "EOF Mean ratio", "Bootstrap EOF MR",
                            "EOF MR SE", "EOF MR lower 95% CI", "EOF MR upper 95% CI",
                            "EOF mean difference", "Bootstrap EOF MD",
                            "EOF MD SE", "EOF MD lower 95% CI", "EOF MD upper 95% CI"))
  } else {
    colnames(resultdf) <- c("t", "Interv.", "Est. EOF mean", "EOF mean ratio",
                            "EOF mean difference", "Obs. EOF mean")
    setcolorder(resultdf, c("t", "Interv.", "Obs. EOF mean", "Est. EOF mean",
                            "EOF mean ratio", "EOF mean difference"))
  }

  fits <- fitcov
  fits[[length(fits) + 1]] <- fitY

  # Add list of coefficients for covariates, outcome variable, and competing event
  # variable (if any) to results output
  coeffs <- lapply(fits, FUN = function(fit){
    if (length(fit) == 2){
      return (list(stats::coefficients(fit[[1]]), stats::coefficients(fit[[2]])))
    } else {
      return (stats::coefficients(fit))
    }
  })
  coeffs <-
    stats::setNames(coeffs, c(as.vector(lapply(covnames,
                                               FUN = function(name){paste("fit", name, sep = "")})), "fitY"))

  rmses <- lapply(1:length(fits), FUN = rmse_calculate, fits = fits, covnames = covnames,
                  covtypes = covtypes, obs_data = obs_data, outcome_name = outcome_name,
                  time_name = time_name, restrictions = restrictions,
                  yrestrictions = yrestrictions, compevent_restrictions = compevent_restrictions)
  rmses <-
    stats::setNames(rmses,
                    c(as.vector(lapply(covnames,
                                       FUN = function(name){paste("rmse", name, sep = "")})),
                      paste("rmse", outcome_name, sep="")))


  # Create header
  header <- get_header(int_descript, sample_size, nsimul, nsamples, ref_int)

  if (sim_data_b){
    sim_data <- pools
    if (!is.na(int_descript[1])){
      names(sim_data) <- int_descript
    }
  } else {
    sim_data <- NA
  }

  res <- list(
    result = resultdf,
    coeffs = coeffs,
    rmses = rmses,
    sim_data = sim_data,
    time_name = time_name,
    time_points = time_points,
    covnames = covnames,
    covtypes = covtypes,
    dt_cov_plot = plot_info$dt_cov_plot,
    dt_out_plot = plot_info$dt_out_plot,
    header = header
  )
  class(res) <- "gformula_continuous_eof"
  return (res)
}

#' Estimation of Binary End-of-Follow-Up Outcome Under the Parametric G-Formula
#'
#' Based on an observed data set, this function estimates the outcome probability at
#' end-of-follow-up under multiple user-specified interventions using the parametric g-formula. See the documentation of (TO DO - link to vignette) for
#' further details concerning the application and implementation of the parametric g-formula.
#'
#' @param id                      Character string specifying the name of the ID variable in \code{obs_data}. The default is \code{"id"}.
#' @param time_points             Number of time points to simulate.
#' @param obs_data                Data table containing the observed data.
#' @param seed                    Starting seed for simulations and bootstrapping.
#' @param nsimul                  Number of subjects for whom to simulate data. By default, this argument is set
#'                                equal to the number of subjects in \code{obs_data}.
#' @param time_name               Character string specifying the name of the time variable in \code{obs_data}. The default is \code{"t0"}.
#' @param outcome_name            Character string specifying the name of the outcome variable in \code{obs_data}. The default is \code{"Y"}.
#' @param intvars                 Vector of character strings specifying the names of the variables to be intervened
#'                                on in each round of the simulation. The default is \code{NA}.
#' @param interventions           List of vectors. Each vector contains a function
#'                                implementing a particular intervention, optionally
#'                                followed by one or more "intervention values" (i.e.,
#'                                integers used to specify the treatment regime).
#'                                The default is \code{NA}.
#' @param int_descript            Vector of character strings, each describing an intervention. It must
#'                                be in same order as the entries in \code{interventions} and should not include a description for the natural course.
#' @param ref_int                 Integer denoting the intervention to be used as the
#'                                reference for calculating the end-of-follow-up mean ratio. 0 denotes the
#'                                natural course, while subsequent integers denote user-specified
#'                                interventions in the order that they are
#'                                named in \code{interventions}. The default is 0.
#' @param covnames                Vector of character strings specifying the names of the time-varying covariates in \code{obs_data}.
#' @param covtypes                Vector of character strings specifying the "type" of each time-varying covariate included in \code{covnames}. The possible "types" are: \code{"binary"}, \code{"normal"}, \code{"categorical"}, \code{"bounded normal"}, \code{"zero-inflated normal"}, \code{"truncated normal"}, and \code{"absorbing"}.
#' @param covparams               List of vectors, where each vector contains information for
#'                                one parameter used in the modeling of the time-varying covariates (e.g.,
#'                                model statement, family, link function, etc.). Each vector
#'                                must be the same length as \code{covnames} and in the same order.
#'                                If a parameter is not required for a certain covariate, it
#'                                should be set to \code{NA} at that index.
#' @param covfits_custom          Vector containing custom fit functions for time-varying covariates that
#'                                do not fall within the pre-defined covariate types. It should be in
#'                                the same order \code{covnames}. If a custom fit function is not
#'                                required for a particular covariate (e.g., if the first
#'                                covariate is of type \code{"binary"} but the second is of type \code{"custom"}), then that
#'                                index should be set to \code{NA}. The default is \code{NA}.
#' @param covpredict_custom       Vector containing custom prediction functions for time-varying
#'                                covariates that do not fall within the pre-defined covariate types.
#'                                It should be in the same order as \code{covnames}. If a custom
#'                                prediction function is not required for a particular
#'                                covariate, then that index should be set to \code{NA}. The default is \code{NA}.
#' @param basecovs                Vector of character strings specifying the names of baseline covariates in \code{obs_data}. These covariates are not simulated using a model but rather carry their value over all time points from the first time point of \code{obs_data}. These covariates should not be included in \code{covnames}. The default is \code{NA}.
#' @param histvars                Vector of character strings specifying the names of the variables for which history functions
#'                                are to be applied. The default is \code{NA}.
#' @param histories               Vector of history functions to apply to the variables specified in \code{histvars}. The default is \code{NA}.
#' @param ymodel                  Model statement for the outcome variable.
#' @param visitprocess            List of vectors. Each vector contains as its first entry
#'                                the covariate name of a visit process; its second entry
#'                                the name of a covariate whose modeling depends on the
#'                                visit process; and its third entry the maximum number
#'                                of consecutive visits that can be missed before an
#'                                individual is censored. The default is \code{NA}.
#' @param yrestrictions           List of vectors. Each vector containins as its first entry
#'                                a condition and its second entry an integer. When the
#'                                condition is \code{TRUE}, the outcome variable is simulated
#'                                according to the fitted model; when the condition is \code{FALSE},
#'                                the outcome variable takes on the value in the second entry.
#'                                The default is \code{NA}.
#' @param restrictions            List of vectors. Each vector contains as its first entry
#'                                the covariate affected by the restriction; its second entry
#'                                the condition that must be \code{TRUE} for the covariate to be
#'                                modeled; its third entry a function that executes other
#'                                specific actions based on the condition; and its fourth
#'                                entry some value used by the function. The default is \code{NA}.
#' @param nsamples                Integer specifying the number of bootstrap samples to generate.
#'                                The default is 0.
#' @param parallel                Logical scalar indicating whether to parallelize simulations of
#'                                different interventions to multiple cores.
#' @param ncores                  Integer specifying the number of cores to use in parallel
#'                                simulation.
#' @param sim_data_b              Logical scalar indicating whether to return the simulated data set. If bootstrap samples are used (i.e., \code{nsamples} is set to a value greater than 0), this argument must be set to \code{FALSE}. The default is \code{FALSE}.
#' @param ...                     Other arguments, which are passed to the functions in \code{covpredict_custom}.
#'
#' @return An object of class \code{gformula_binary_eof}. The object is a list with the following components:
#' \item{result}{Results table containing the estimated outcome probability for all interventions (inculding natural course) at the last time point. If bootstrapping was used, the results table includes the bootstrap end-of-follow-up mean ratio, standard error, and 95\% confidence interval.}
#' \item{coeffs}{A list of the coefficients of the fitted models.}
#' \item{rmses}{A list of root mean square error (RMSE) values of the fitted models.}
#' \item{sim_data}{A list of data tables of the simulated data. Each element in the list corresponds to one of the interventions. If the argument \code{sim_data_b} is set to \code{FALSE}, a value of \code{NA} is given.}
#' \item{...}{Some additional elements.}
#'
#' The results for the g-formula simulation under various interventions for the last time point are printed with the \code{\link{print.gformula_binary_eof}} function. To generate graphs comparing the mean estimated and observed covariate values over time, use the \code{\link{plot.gformula_binary_eof}} function.
#'
#' @references Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period: application to the healthy worker survivor effect. Mathematical Modelling. 1986;7:1393–1512. [Errata (1987) in Computers and Mathematics with Applications 14, 917.-921. Addendum (1987) in Computers and Mathematics with Applications 14, 923-.945. Errata (1987) to addendum in Computers and Mathematics with Applications 18, 477.].
#' @examples
#' ## TO ADD AFTER PAPER IS FINALIZED
#' @import data.table
#' @importFrom dplyr funs n summarise_all group_by %>%
#' @export
gformula_binary_eof <- function(obs_data, id = 'id', time_points,
                                time_name = 't0', covnames, covtypes, covparams,
                                covfits_custom = NA, covpredict_custom = NA,
                                histvars = NA, histories = NA, basecovs = NA,
                                outcome_name = 'Y', ymodel, intvars = NA,
                                interventions = NA, int_descript = NA,
                                ref_int = 0, visitprocess = NA,
                                restrictions = NA, yrestrictions = NA,
                                nsimul = NA, sim_data_b = FALSE, seed,
                                nsamples = 0, parallel = FALSE, ncores = NA,
                                ...){

  colnames(obs_data)[colnames(obs_data) == id] <- 'id'
  for (history in histories) {
    for (t in 0:max(obs_data[[time_name]])) {
      obs_data <- history(pool = obs_data, histvars = histvars,
                          time_name = time_name, t = t, id = 'id',
                          max_visits = NA)[[length(histvars)]]
    }
  }
  colnames(obs_data)[colnames(obs_data) == 'id'] <- id

  outcome_type <- 'binary_eof'
  comprisk = FALSE
  compevent_model = NA
  compevent_name = NA
  compevent_restrictions = NA
  hazardratio <- FALSE
  intcomp <- NA
  sample_size <- length(unique(obs_data[[id]]))

  # for (i in 1:length(covnames)){
  #   if (covtypes[i] == 'absorbing'){
  #     restrictions <- c(restrictions[!is.na(restrictions)],
  #                       list(c(covnames[i], paste("lag_", covnames[i], "==0", sep = ""),
  #                              carry_forward, 1)))
  #     covtypes[i] <- 'binary'
  #   }
  # }

  max_visits <- NA
  if (!is.na(visitprocess[[1]][[1]])){
    for (vp in visitprocess){
      restrictions <- c(restrictions[!is.na(restrictions)],
                        list(c(vp[1], paste("visit_sum_", vp[3], "_", vp[1], "!=0", sep = ""),
                               simple_restriction, 1),
                             c(vp[2], paste(vp[1], "==1", sep = ""), carry_forward)))
      if (is.na(max_visits)){
        max_visits <- as.numeric(vp[3])
      } else {
        max_visits <- c(max_visits, as.numeric(vp[3]))
      }
    }
  }

  # Create 1-indexed numerical IDs for observed datasets
  ids <- as.data.table(sort(unique(obs_data[[id]])))
  ids[, "newid" := seq_len(.N)]
  setkeyv(obs_data, id)
  obs_data <- obs_data[J(ids), allow.cartesian = TRUE]

  # Set default number of simulated individuals to equal number of individuals in
  # observed dataset
  if (is.na(nsimul)){
    nsimul <- length(unique(obs_data$newid))
  }

  error_catch(id, nsimul, intvars, interventions, int_descript,
              covnames, covtypes, basecovs,
              histvars, histories, compevent_model, hazardratio, intcomp,
              time_points, outcome_type, time_name, obs_data, parallel, ncores,
              nsamples, sim_data_b)

  # Generate seeds for simulations and bootstrapping
  set.seed(seed)
  newseeds <- sample.int(2^30, size = nsamples + 1)
  subseed <- newseeds[1]
  bootseeds <- newseeds[2:(nsamples + 1)]

  # Determine ranges of observed covariates and outcome
  ranges <- lapply(1:length(covnames), FUN = function(i){
    if (covtypes[i] == 'normal' || covtypes[i] == 'bounded normal' ||
        covtypes[i] == 'zero-inflated normal' || covtypes[i] == 'truncated normal') {
      range(obs_data[[covnames[i]]])
    } else if (covtypes[i] == 'zero-inflated normal'){
      range(obs_data[obs_data[[covnames[i]]] > 0][[covnames[i]]])
    } else {
      NA
    }
  })
  yrange <- range(obs_data[[outcome_name]])

  # Fit models to covariates and outcome variable
  fitcov <- pred_fun_cov(covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, restrictions = restrictions,
                         time_name = time_name, obs_data = obs_data)
  fitY <- pred_fun_Y(ymodel, yrestrictions, outcome_type, outcome_name, time_name, obs_data)

  len <- length(unique(obs_data$newid))
  # If the number of user desired simulations differs from the number of individuals in
  # the observed dataset, sample the desired number of observed IDs with replacement
  if (nsimul < len){
    ids <- as.data.table(sort(sample(unique(obs_data$newid), nsimul, replace = TRUE)))
    colnames(ids) <- "newid"
    ids[, 'sid' := seq_len(.N)]
    obs_data <- merge(ids, obs_data, all.x = TRUE, by = "newid")
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  } else if (nsimul > len){
    ids <- data.frame(sample(unique(obs_data$newid), nsimul, replace = TRUE))
    ids[, 'newid' := 1:nsimul]
    colnames(ids) <- c("newid", "sid")
    setkeyv(obs_data, "newid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE]  # create the new data set names "sample"
    obs_data[, 'newid' := obs_data$sid]
    obs_data[, 'sid' := NULL]
  }

  # Add natural course to list of interventions
  if (!is.na(interventions)[[1]][[1]]){
    comb_interventions <- c(list(c(natural)), interventions)
    comb_intvars <- c(list('none'), intvars)
  } else {
    comb_interventions <- list(c(natural))
    comb_intvars <- list('none')
  }

  if (parallel){
    if (is.na(ncores)){
      ncores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(ncores)
    parallel::clusterCall(cl, library, package = 'data.table', character.only = TRUE)
    parallel::clusterCall(cl, library, package = 'dplyr', character.only = TRUE)
    if ('categorical' %in% covtypes){
      parallel::clusterCall(cl, library, package = 'nnet', character.only = TRUE)
    }
    if ('truncated normal' %in% covtypes){
      parallel::clusterCall(cl, library, package = 'truncreg', character.only = TRUE)
    }

    clusterAssign(cl, "make_histories", make_histories)
    clusterAssign(cl, "rmse_calculate", rmse_calculate)
    clusterAssign(cl, "pred_fun_cov", pred_fun_cov)
    clusterAssign(cl, "fit_glm", fit_glm)
    clusterAssign(cl, "fit_multinomial", fit_multinomial)
    clusterAssign(cl, "fit_trunc_normal", fit_trunc_normal)
    clusterAssign(cl, "pred_fun_Y", pred_fun_Y)
    clusterAssign(cl, "pred_fun_D", pred_fun_D)
    clusterAssign(cl, "natural", natural)
    clusterAssign(cl, "static", static)
    clusterAssign(cl, "threshold", threshold)

    # Note warning suppression when submitting to CRAN, reference 'future' package
    # as doing the same thing for similar issue
    suppressWarnings(parallel::clusterExport(cl, as.vector(utils::lsf.str())))

    pools <- parallel::parLapply(cl, 1:length(comb_interventions), simulate,
                                 fitcov = fitcov, fitY = fitY, fitD = NA,
                                 yrestrictions = yrestrictions,
                                 compevent_restrictions = compevent_restrictions,
                                 restrictions = restrictions,
                                 outcome_name = outcome_name, compevent_name = compevent_name,
                                 time_name = time_name,
                                 intvars = comb_intvars, interventions = comb_interventions,
                                 histvars = histvars, histories = histories,
                                 covparams = covparams, covnames = covnames, covtypes = covtypes,
                                 covfits_custom = covfits_custom, basecovs = basecovs,
                                 comprisk = comprisk, ranges = ranges,
                                 yrange = yrange, compevent_range = NA,
                                 outcome_type = outcome_type,
                                 subseed = subseed, time_points = time_points,
                                 obs_data = obs_data, parallel = parallel, ...)

  } else {
    pools <- lapply(1:length(comb_interventions), FUN = function(i){
      simulate(fitcov = fitcov, fitY = fitY, fitD = NA,
               yrestrictions = yrestrictions,
               compevent_restrictions = compevent_restrictions,
               restrictions = restrictions,
               outcome_name = outcome_name, compevent_name = compevent_name,
               time_name = time_name,
               intvars = comb_intvars[[i]], interventions = comb_interventions[[i]],
               histvars = histvars, histories = histories,
               covparams = covparams, covnames = covnames, covtypes = covtypes,
               covfits_custom = covfits_custom, basecovs = basecovs, comprisk = comprisk,
               ranges = ranges, yrange = yrange, compevent_range = NA,
               outcome_type = outcome_type,
               subseed = subseed, time_points = time_points,
               obs_data = obs_data, parallel = parallel, ...)
    })
  }

  nat_pool <- pools[[1]] # Natural course data
  pools <- pools[-1] # List of intervention datasets

  # Initialize results matrices
  result_ratio <- result_diff <- int_result <- rep(NA, length(pools) + 1)

  # Calculate mean outcome over all subjects at each time for natural course
  nat_result <- mean(nat_pool$Py, na.rm = TRUE)

  if (ref_int == 0){
    # Set reference intervention to the natural course
    ref_mean <- nat_result
  } else {
    # Set reference intervention as specified
    # Calculate mean outcome over all subjects at each time for this intervention
    ref_mean <- mean(pools[[ref_int]]$Py, na.rm = TRUE)
  }

  # Compile results
  int_result[1] <- nat_result
  # Calculate mean risk over all subjects at each time for all interventions other than
  # the natural course
  int_result[-1] <- sapply(pools, FUN = function(pool){mean(pool$Py, na.rm = TRUE)})
  result_ratio <- int_result / ref_mean
  result_diff <- int_result - ref_mean

  if (nsamples > 0){
    if (parallel){
      parallel::clusterExport(cl, 'simulate')
      final_bs <- parallel::parLapply(cl, 1:nsamples, bootstrap_helper, time_points = time_points,
                                      obs_data = obs_data, bootseeds = bootseeds,
                                      intvars = intvars, interventions = interventions, ref_int = ref_int,
                                      covparams = covparams, covnames = covnames, covtypes = covtypes,
                                      covfits_custom = covfits_custom, basecovs = basecovs, ymodel = ymodel,
                                      histvars = histvars, histories = histories,
                                      comprisk = comprisk, compevent_model = compevent_model,
                                      yrestrictions = yrestrictions,
                                      compevent_restrictions = compevent_restrictions,
                                      restrictions = restrictions, outcome_type = outcome_type,
                                      ranges = ranges, yrange = yrange, compevent_range = NA,
                                      time_name = time_name, outcome_name = outcome_name,
                                      compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                                      max_visits = max_visits, ...)

      parallel::stopCluster(cl)

    } else {
      final_bs <- lapply(1:nsamples, FUN = bootstrap_helper, time_points = time_points,
                         obs_data = obs_data, bootseeds = bootseeds,
                         intvars = intvars, interventions = interventions, ref_int = ref_int,
                         covparams = covparams, covnames = covnames, covtypes = covtypes,
                         covfits_custom = covfits_custom, basecovs = basecovs,
                         ymodel = ymodel,
                         histvars = histvars, histories = histories,
                         comprisk = comprisk, compevent_model = compevent_model,
                         yrestrictions = yrestrictions,
                         compevent_restrictions = compevent_restrictions,
                         restrictions = restrictions,
                         outcome_type = outcome_type,
                         ranges = ranges, yrange = yrange, compevent_range = NA,
                         time_name = time_name, outcome_name = outcome_name,
                         compevent_name = compevent_name, parallel = parallel, ncores = ncores,
                         max_visits = max_visits)
    }
    comb_result <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$Result))
    }))
    comb_MR <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultRatio))
    }))
    comb_MD <- rbindlist(lapply(final_bs, FUN = function(m){
      as.data.table(t(m$ResultDiff))
    }))
    # Creat dataframe of risk ratios
    #comb <- rbindlist(lapply(final_bs, FUN = function(m){data.table(t(m))}))

    comb_result$t0 <- comb_MR$t0 <- comb_MD$t0 <- max(obs_data[[time_name]])
    # Calculate mean and 95% confidence interval of risk ratios
    CI_result <- comb_result %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))
    CI_MR <- comb_MR %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))
    CI_MD <- comb_MD %>% group_by(t0) %>%
      summarise_all(funs(mean, stats::sd, se = stats::sd(.) / sqrt(n()),
                         CI_lower = mean(.) - 1.96*stats::sd(.) / sqrt(n()),
                         CI_upper = mean(.) + 1.96*stats::sd(.) / sqrt(n())))

    # If user does not choose to bootstrap, indicate that no confidence interval has
    # been calculated
  } else {
    CI_result <- CI_MR <- CI_MD <- 'not available'
  }

  plot_info <- get_plot_info(outcome_name = outcome_name,
                             compevent_name = compevent_name,
                             time_name = time_name,
                             time_points = time_points,
                             covnames = covnames,
                             covtypes = covtypes,
                             nat_pool = nat_pool,
                             nat_result = nat_result,
                             comprisk = comprisk,
                             outcome_type = outcome_type,
                             obs_data = obs_data)
  obs_results <- plot_info$obs_results

  # Generate results table
  if (!is.na(interventions)[[1]][[1]]){
    resultdf <- lapply(1:length(int_result), function(k){
      if (nsamples > 0){
        data.table(t = time_points - 1, Intervention = k - 1,
                   EOFMean = int_result[k],
                   EOFMean_bootmean = eval(parse(text = paste("CI_result$V", k, "_mean", sep = ""))),
                   EOFMean_SE = eval(parse(text = paste("CI_result$V", k, "_se", sep = ""))),
                   EOFMean_CI_LL95 = eval(parse(text = paste("CI_result$V", k, "_CI_lower", sep = ""))),
                   EOFMean_CI_UL95 = eval(parse(text = paste("CI_result$V", k, "_CI_upper", sep = ""))),
                   EOFMeanRatio = result_ratio[k],
                   MR_bootmean = eval(parse(text = paste("CI_MR$V", k, "_mean", sep = ""))),
                   MR_SE = eval(parse(text = paste("CI_MR$V", k, "_se", sep = ""))),
                   MR_CI_LL95 = eval(parse(text = paste("CI_MR$V", k, "_CI_lower", sep = ""))),
                   MR_CI_UL95 = eval(parse(text = paste("CI_MR$V", k, "_CI_upper", sep = ""))),
                   EOFMeanDiff = result_diff[k],
                   MD_bootmean = eval(parse(text = paste("CI_MD$V", k, "_mean", sep = ""))),
                   MD_SE = eval(parse(text = paste("CI_MD$V", k, "_se", sep = ""))),
                   MD_CI_LL95 = eval(parse(text = paste("CI_MD$V", k, "_CI_lower", sep = ""))),
                   MD_CI_UL95 = eval(parse(text = paste("CI_MD$V", k, "_CI_upper", sep = ""))))
      } else {
        data.table(t = time_points - 1, Intervention = k - 1, OutcomeEOFMean = int_result[k],
                   EOFMeanRatio = result_ratio[k], EOFMeanDiff = result_diff[k])
      }
    })
    resultdf <- rbindlist(resultdf)
  } else {
    if (nsamples > 0){
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             EOFMean = int_result, EOFMean_bootmean = CI_result$mean,
                             EOFMean_SE = CI_result$se, EOFMean_CI_LL95 = CI_result$CI_lower,
                             EOFMean_CI_UL95 = CI_result$CI_upper, EOFMeanRatio = result_ratio,
                             MR_bootmean = CI_MR$mean, MR_SE = CI_MR$se,
                             MR_CI_LL95 = CI_MR$CI_lower, MR_CI_UL95 = CI_MR$CI_upper,
                             EOFMeanDiff = result_diff, EOFMD_mean = CI_MD$mean,
                             EOFMD_se = CI_MD$se,
                             MD_CI_LL95 = CI_MD$CI_lower, MD_CI_UL95 = CI_MD$CI_upper)
    } else {
      resultdf <- data.table(t = 0:(time_points - 1), Intervention = 0,
                             OutcomeEOFMean = int_result, EOFMeanRatio = result_ratio,
                             EOFMeanDiff = result_diff)
    }
  }

  resultdf[, 'obs_risk' := c(utils::tail(obs_results[[2]], 1), rep(NA, dim(resultdf)[1] - 1))]

  if (nsamples > 0){
    colnames(resultdf) <- c("t", "Interv.", "Est. EOF mean",
                            "Bootstrap EOF mean", "EOF mean SE",
                            "EOF mean lower 95% CI", "EOF mean upper 95% CI",
                            "EOF Mean ratio", "Bootstrap EOF MR",
                            "EOF MR SE", "EOF MR lower 95% CI", "EOF MR upper 95% CI",
                            "EOF mean difference", "Bootstrap EOF MD",
                            "EOF MD SE", "EOF MD lower 95% CI", "EOF MD upper 95% CI",
                            "Obs. EOF mean")
    setcolorder(resultdf, c("t", "Interv.", "Obs. EOF mean", "Est. EOF mean",
                            "Bootstrap EOF mean", "EOF mean SE",
                            "EOF mean lower 95% CI", "EOF mean upper 95% CI",
                            "EOF Mean ratio", "Bootstrap EOF MR",
                            "EOF MR SE", "EOF MR lower 95% CI", "EOF MR upper 95% CI",
                            "EOF mean difference", "Bootstrap EOF MD",
                            "EOF MD SE", "EOF MD lower 95% CI", "EOF MD upper 95% CI"))
  } else {
    colnames(resultdf) <- c("t", "Interv.", "Est. EOF mean", "EOF mean ratio",
                            "EOF mean difference", "Obs. EOF mean")
    setcolorder(resultdf, c("t", "Interv.", "Obs. EOF mean", "Est. EOF mean",
                            "EOF mean ratio", "EOF mean difference"))
  }

  fits <- fitcov
  fits[[length(fits) + 1]] <- fitY

  # Add list of coefficients for covariates, outcome variable, and competing event
  # variable (if any) to results output
  coeffs <- lapply(fits, FUN = function(fit){
    if (length(fit) == 2){
      return (list(stats::coefficients(fit[[1]]), stats::coefficients(fit[[2]])))
    } else {
      return (stats::coefficients(fit))
    }
  })
  coeffs <-
    stats::setNames(coeffs, c(as.vector(lapply(covnames,
                                               FUN = function(name){paste("fit", name, sep = "")})), "fitY"))

  rmses <- lapply(1:length(fits), FUN = rmse_calculate, fits = fits, covnames = covnames,
                  covtypes = covtypes, obs_data = obs_data, outcome_name = outcome_name,
                  time_name = time_name, restrictions = restrictions,
                  yrestrictions = yrestrictions, compevent_restrictions = compevent_restrictions)
  rmses <-
    stats::setNames(rmses,
                    c(as.vector(lapply(covnames,
                                       FUN = function(name){paste("rmse", name, sep = "")})),
                      paste("rmse", outcome_name, sep="")))


  # Create header
  header <- get_header(int_descript, sample_size, nsimul, nsamples, ref_int)

  if (sim_data_b){
    sim_data <- pools
    if (!is.na(int_descript[1])){
      names(sim_data) <- int_descript
    }
  } else {
    sim_data <- NA
  }

  res <- list(
    result = resultdf,
    coeffs = coeffs,
    rmses = rmses,
    sim_data = sim_data,
    time_name = time_name,
    time_points = time_points,
    covnames = covnames,
    covtypes = covtypes,
    dt_cov_plot = plot_info$dt_cov_plot,
    dt_out_plot = plot_info$dt_out_plot,
    header = header
  )
  class(res) <- "gformula_binary_eof"
  return (res)
}
