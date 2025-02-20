###########################################
###########################################
### FUNCTIONS FOR ASSESSING CALIBRATION ###
###########################################
###########################################

###
### Write a function to estimate survival probabilities based on a fitted model (fit) and baseline hazard (bhaz).
### Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE).
###
est_surv <- function(newdata, fit, bhaz = NULL, t){

  ### Get bhaz if not specified
  if (is.null(bhaz)){bhaz <- survival::basehaz(fit, centered = TRUE)}

  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")

  ### Get the linear predictor for ne wdata
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= t))]))

  return(surv)

}

###
### Assessing calibration using proportional hazards regression approach (graphical calibration curves, Austin et al)
###
est_calib_ph <- function(data, fit, bhaz, t, surv = NULL, nk = 4, pred_plot_range = NULL){

  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }

  ### Add complementary log-log of predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))

  ### Fit calibration model
  fit.calib <- survival::coxph(survival::Surv(time, status) ~ rms::rcs(cloglog, nk), data = data)
  bhaz.calib <- survival::basehaz(fit.calib, centered = TRUE)

  ###
  ### Generate predicted observed values,
  ###

  ### Do so over the range of values pred_plot_range, if specified
  if (!is.null(pred_plot_range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred_plot_range, cloglog = log(-log(1 - pred_plot_range)))

    ### Calculate predicted observed values over pred_plot_range
    pred.obs <- 1 - est_surv(newdata = tmp.data, fit = fit.calib, bhaz = bhaz.calib, t = t)

    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)

  } else {
    ### Otherwise, do so for all individuals in data
    pred.obs <- est_surv(newdata = data, fit = fit.calib, bhaz = bhaz.calib, t = t)

    ### Create plot data
    plot.data <- data.frame("pred.obs" = 1 - pred.obs, "pred" = 1 - data$surv)
  }

  ### Calculate ICI, E50 and E90
  ### Only do so if pred_plot_range isn't specified, otherwise these are meaningless
  if (!is.null(pred_plot_range)){
    ICI <- NULL
    E50 <- NULL
    E90 <- NULL
  } else {
    ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
    E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
    E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9, na.rm = TRUE))
  }

  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle("Proportional hazards") +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")

  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "plotdata" = plot.data)

  return(output.object)

}

###
### Function to estimate calibration using IPCW
###

### Assessing calibration in individuals who are uncensored at time t.
### See [@Pate2024] for an example of this approach applied to multistate models,
### although there are many more resources on inverse probability of censoring weights more generally.

### First Write a function to estimate survival probability for a given model and fitted baseline hazard,
### but survival probability is probability at the time they had their event.
### Variable eventtime is a character string, of the name of the variable denoting the event time.
### Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE).
### This is used to calculate weights for individuals who have an event, and therefore we need to calculate
### probability of being unesored at the time they had their event.
est_surv_eventtime <- function(newdata, bhaz, fit, eventtime, type = "coxph"){

  if (type == "coxph"){

    ### Get lp
    lp <- predict(fit, newdata = newdata, type = "lp", reference = "sample")

    ### Get bhaz specific for each patient at the time they had their event
    bhaz.pat <- unlist(lapply(1:nrow(newdata), function(x) {bhaz$hazard[max(which(bhaz$time <= newdata[x, eventtime]))]}))

    ### Get expected
    surv <- as.numeric(exp(-exp(lp)*bhaz.pat))

  } else if (type == "flexsurv"){

    # Calculate survival probabilities for these individuals
    surv <- unlist(lapply(1:nrow(newdata),
                          function(x) {dplyr::pull(predict(fit,
                                                           newdata = newdata[x, ],
                                                           times = as.numeric(newdata[x, eventtime]),
                                                           type = "survival"),
                                                   .pred_survival)}
    ))

  }

  return(surv)

}


###
### Function to estimate calibration using IPCW
###
est_calib_ipcw <- function(data, fit, bhaz, t, surv = NULL, cens_max_follow = NULL, nk = 4, cens_formula, type = "coxph", pred_plot_range = NULL){

  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }

  ### Assign a variable for whether individuals have had an event by time of interest
  data <- dplyr::mutate(data,
                        status_time = dplyr::case_when(time <= t ~ status,
                                                       time > t ~ 0))

  ### Estimate weights
  weights <- est_ipcw(data = data, t = t, cens_max_follow = cens_max_follow, cens_formula = cens_formula, type = type)

  ### Merge with data
  data <- merge(data, weights, by.x = "id", by.y = "id")

  ### Reduce data to individuals uncensored at time t
  data <- subset(data, !is.na(ipcw))

  ### Convert predict risk onto logit scale
  data$pred <- 1 - data$surv
  data$pred.logit <- log(data$pred/(1 - data$pred))

  ### Fit weighted calibration model and generate predicted observed
  ## Fit model
  rcs.model.stab <- suppressWarnings(glm(status_time ~ rms::rcs(pred.logit, nk),
                                         family = binomial(link = "logit"),
                                         data = data,
                                         weights = data[, "ipcw.stab"]))
  ## Supress warnings due to using weights in binomial glm results in
  ## "Warning message: In eval(family$initialize) : non-integer #successes in a binomial glm!"

  ###
  ### Generate predicted observed values,
  ###

  ### Do so over the range of values pred_plot_range, if specified
  if (!is.null(pred_plot_range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred_plot_range, "pred.logit" = log(pred_plot_range/(1-pred_plot_range)))

    ### Calculate predicted observed values over pred_plot_range
    pred.obs <- predict(rcs.model.stab, newdata = tmp.data, type = "response")

    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)

  } else {
    ### Generate predicted observed
    data$pred.obs <- predict(rcs.model.stab, newdata = data, type = "response")

    ### Create plot data
    plot.data <- data.frame("pred.obs" = data$pred.obs, "pred" = data$pred)

  }


  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle("IPCW") +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")

  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(plot.data$pred.obs - plot.data$pred), na.rm = TRUE)
  E50 <- median(abs(plot.data$pred.obs - plot.data$pred), na.rm = TRUE)
  E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9, na.rm = TRUE))

  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "plotdata" = plot.data)

  ### Return
  return(output.object)

}

###
### Write a function to estimate IPCW's that can be used in multiple other functions
###
est_ipcw <- function(data, t, cens_max_follow = NULL, cens_formula, type){

  ### If there is informative censoring (i.e. all individuals censored after 11 years)
  ### Don't want to model the effect of everybody becoming censored at the final follow up time,
  ### so set these individuals to be uncensored at time point of interest.
  ### This will not be used in the simulation, but can occur in real data so leaving in the functionality.
  if (!is.null(cens_max_follow)){
    data <- dplyr::mutate(data,
                          cens_indicator = dplyr::case_when(cens_time < cens_max_follow + 2 ~ cens_indicator,
                                                            cens_time >= cens_max_follow + 2 ~ 0),
                          cens_time = dplyr::case_when(cens_time < cens_max_follow + 2 ~ cens_time,
                                                       cens_time >= cens_max_follow + 2 ~ cens_max_follow + 2))
  }

  ### Create censoring model and create survival predictions from this model, probability of being uncensored at time t
  if (type == "coxph"){

    ### Create censoring model
    cens.model <- survival::coxph(cens_formula, data = data, model = TRUE)
    cens.model.int <- survival::coxph(survival::Surv(cens_time, cens_indicator) ~ 1, data = data)

    ### Create survival predictions from this model, probability of being uncensored at time t
    cens.bhaz <- survival::basehaz(cens.model, centered = TRUE)
    data$surv.cens <- est_surv(newdata = data, fit = cens.model, bhaz = cens.bhaz, t = t)

  } else if (type == "flexsurv"){

    ### Create censoring model
    cens.model <- flexsurv::flexsurvreg(cens_formula, data = data, dist = "weibull")
    cens.model.int <- flexsurv::flexsurvreg(survival::Surv(cens_time, cens_indicator) ~ 1, data = data, dist = "weibull")

    ### Create survival predictions from this model, probability of being uncensored at time t
    data$surv.cens <- dplyr::pull(predict(cens.model, newdata = data, times = t, type = "survival"), .pred_survival)

  }

  ### For individual who have an event prior to or equal to time t (and are therefore uncensored at time t),
  ### we want probability of being uncensored at t_event. For people censored before or on time t,
  ### they will not be included in the analysis, as they are censored.
  obs.event.prior <- which(data$cens_time <= t & data$cens_indicator == 0)
  obs.censored.prior <- which(data$cens_time <= t & data$cens_indicator == 1)

  ### Calculate weights for individuals who have events prior to time t
  surv.event.prior <- est_surv_eventtime(newdata = data[obs.event.prior, ],
                                         bhaz = cens.bhaz,
                                         fit = cens.model,
                                         eventtime = "cens_time",
                                         type = type)

  ### Replace these values in dataset
  data$surv.cens[obs.event.prior] <- surv.event.prior

  ### Create weights
  data$ipcw <- 1/data$surv.cens

  ### Stabilise
  ### Create censoring model and create survival predictions from this model, probability of being uncensored at time t
  if (type == "coxph"){
    # Get bhaz
    cens.bhaz.int <- survival::basehaz(cens.model.int, centered = TRUE)
    # Get survival probabilities from intercept only model
    data$ipcw.numer <- est_surv(newdata = data, fit = cens.model.int, bhaz = cens.bhaz.int, t = t)
    # Apply stabilisation
    data$ipcw.stab <- data$ipcw.numer*data$ipcw
  } else if (type == "flexsurv"){
    # Get survival probabilities from intercept only model
    data$ipcw.numer <- dplyr::pull(predict(cens.model.int, newdata = data, times = t, type = "survival"), .pred_survival)
    # Apply stabilisation
    data$ipcw.stab <- data$ipcw.numer*data$ipcw
  }

  ### Cap the weights at p1 and p99
  data$ipcw <- pmax(data$ipcw, as.numeric(quantile(data$ipcw, p = .01)))
  data$ipcw <- pmin(data$ipcw, as.numeric(quantile(data$ipcw, p = .99)))

  data$ipcw.stab <- pmax(data$ipcw.stab, as.numeric(quantile(data$ipcw.stab, p = .01)))
  data$ipcw.stab <- pmin(data$ipcw.stab, as.numeric(quantile(data$ipcw.stab, p = .99)))

  ### Cap at 100
  data$ipcw <- pmin(data$ipcw, 100)
  data$ipcw.stab <- pmin(data$ipcw.stab, 100)

  ### Assign NA's to individuals censored prior to time t
  data$ipcw[obs.censored.prior] <- NA
  data$ipcw.stab[obs.censored.prior] <- NA

  ### Filter
  data <- dplyr::select(data, id, ipcw, ipcw.stab)

  return(data)

}


###
### Write a function to estimate IPCW's that can be used to group individuals in the pseudo-value approach
### Note here: we do not want to estimate the weights at min(t, t_abs)
### We are just interested in the probability of being censored assuming individuals are followed up for the same amount of time
###
est_ipcw_pv <- function(data, t, cens_max_follow = NULL, cens_formula, type){

  ### If there is informative censoring (i.e. all individuals censored after 11 years)
  ### Don't want to model the effect of everybody becoming censored at the final follow up time,
  ### so set these individuals to be uncensored at time point of interest.
  ### This will not be used in the simulation, but can occur in real data so leaving in the functionality.
  if (!is.null(cens_max_follow)){
    data <- dplyr::mutate(data,
                          cens_indicator = dplyr::case_when(cens_time < cens_max_follow + 2 ~ cens_indicator,
                                                            cens_time >= cens_max_follow + 2 ~ 0),
                          cens_time = dplyr::case_when(cens_time < cens_max_follow + 2 ~ cens_time,
                                                       cens_time >= cens_max_follow + 2 ~ cens_max_follow + 2))
  }

  ### Create censoring model and create survival predictions from this model, probability of being uncensored at time t
  if (type == "coxph"){

    ### Create censoring model
    cens.model <- survival::coxph(cens_formula, data = data, model = TRUE)
    cens.model.int <- survival::coxph(survival::Surv(cens_time, cens_indicator) ~ 1, data = data)

    ### Create survival predictions from this model, probability of being uncensored at time t
    cens.bhaz <- survival::basehaz(cens.model, centered = TRUE)
    data$surv.cens <- est_surv(newdata = data, fit = cens.model, bhaz = cens.bhaz, t = t)

  } else if (type == "flexsurv"){

    ### Create censoring model
    cens.model <- flexsurv::flexsurvreg(cens_formula, data = data, dist = "weibull")
    cens.model.int <- flexsurv::flexsurvreg(survival::Surv(cens_time, cens_indicator) ~ 1, data = data, dist = "weibull")

    ### Create survival predictions from this model, probability of being uncensored at time t
    data$surv.cens <- dplyr::pull(predict(cens.model, newdata = data, times = t, type = "survival"), .pred_survival)

  }

  ### Create weights
  data$ipcw <- 1/data$surv.cens

  ### Filter
  data <- dplyr::select(data, id, ipcw)

  return(data)

}



###
### Function to estimate calibration using pseudo-value appraoch, where pseudo-values estimated using
### the kaplan-meier estimate of survival
###

### First write a function which will estimate pseudo-values for a given dataset at a time t
est_pv <- function(df, t){

  ### Create prodlim object
  prodlim.obj <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = df)

  ### Calculate psuedo-values
  pv <- prodlim::jackknife(prodlim.obj, times = t)

  return(as.numeric(pv))
}

### Now write the main function
est_calib_pv <- function(data, fit, bhaz, t, surv = NULL, nk = 4, split_n_groups = 1,
                         group_by = c("predrisk", "ipcw"),
                         ipcw_cens_max_follow = NULL, ipcw_cens_formula = NULL, ipcw_type = "coxph",
                         group.var = NULL, pred_plot_range = NULL){

  ### Match.arg
  group_by <- match.arg(group_by)

  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }

  ### Add predicted risk
  data$pred <- 1 - data$surv

  ### Split continuously if group.var is null
  if (is.null(group.var)){
    ### Split data by survival probabilities
    if (group_by == "predrisk"){
      ### Split
      data.pv.split <- split(data,
                             cut(data$surv,
                                 c(-Inf, quantile(data$surv, probs = 1:split_n_groups/split_n_groups))))
      ### Split data by ipcws
    } else if (group_by == "ipcw"){
      ### Estimate weights
      weights <- est_ipcw_pv(data = data, t = t, cens_max_follow = ipcw_cens_max_follow, cens_formula = ipcw_cens_formula, type = ipcw_type)

      ### Merge with data
      data <- merge(data, weights, by.x = "id", by.y = "id")

      ### Split
      data.pv.split <- split(data,
                             cut(data$ipcw,
                                 c(-Inf, quantile(data$ipcw, probs = 1:split_n_groups/split_n_groups))))
      ### Not uncommon for the final few groups to have no events
      ### Want to combine these, until the final group has an event
      ## Set counting variables
      split_n_groups.ticker <- split_n_groups
      n.events.final.group <- sum(data.pv.split[[split_n_groups.ticker]]$status)
      ## Create loop
      while(n.events.final.group <= 0){

        ## Reduce ticker by 1
        split_n_groups.ticker <- split_n_groups.ticker - 1
        ## combine final two groups
        data.pv.split[[split_n_groups.ticker]] <- rbind(data.pv.split[[split_n_groups.ticker]], data.pv.split[[split_n_groups.ticker+1]])
        ## Remove final group
        data.pv.split <- data.pv.split[-(split_n_groups.ticker+1)]
        ## Calculate n events
        n.events.final.group <- sum(data.pv.split[[split_n_groups.ticker]]$status)
      }

    }
  }

  ### Calculate pseudo-values for each data split (note est_pv calculates pseudo-value for survival prob, so want to take 1 - pv)
  pv <- lapply(data.pv.split, est_pv, t = t)
  pvs <- 1 - do.call("c", pv)
  # print(str(pv))

  ### Also get the corresponding ids
  ids <- lapply(data.pv.split, function(x) {x$id})
  ids <- do.call("c", ids)

  ### Combine
  pvs <- data.frame("pv" = pvs, "id" = ids)
  pvs <- dplyr::arrange(pvs, id)

  ### Add the pseudo-values to validation dataset
  data$pv <- pvs$pv

  ### If pseudo-value is NA, this means all individual had had an event prior to time point t.eval
  ### By definition pseudo-values for these individuals are equal to 1, so assign these
  data$pv[is.na(data$pv)] <- 1

  ### Transform est.surv onto logit scale
  data$pred.logit <- log(data$pred/(1-data$pred))

  ### Fit the model using logit link function
  calib.model.pv <- stats::glm(pv ~ rms::rcs(pred.logit, nk),
                               data = data,
                               family = stats::gaussian(link = "logit"),
                               start = rep(0, nk))

  ###
  ### Generate predicted observed values,
  ###

  ### Do so over the range of values pred_plot_range, if specified
  if (!is.null(pred_plot_range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred_plot_range, "pred.logit" = log(pred_plot_range/(1-pred_plot_range)))

    ### Calculate predicted observed values over pred_plot_range
    pred.obs <- predict(calib.model.pv, newdata = tmp.data, type = "response")

    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)

  } else {
    ### Otherwise, do so for all individuals in data
    pred.obs <- predict(calib.model.pv, newdata = data, type = "response")

    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = data$pred)

  }

  ### Calculate ICI, E50 and E90
  ### Only do so if pred_plot_range isn't specified, otherwise these are meaningless
  if (!is.null(pred_plot_range)){
    ICI <- NULL
    E50 <- NULL
    E90 <- NULL
  } else {
    ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
    E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
    E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9, na.rm = TRUE))
  }

  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(("Pseudo-value")) +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")

  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
  E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
  E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9, na.rm = TRUE))

  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "plotdata" = plot.data)

  return(output.object)

}
