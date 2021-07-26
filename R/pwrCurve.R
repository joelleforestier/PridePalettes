#' Simulate and visualize simultaneous power curves
#'
#' Simulate and visualize individual and simultaneous power at a range of sample sizes of a set of statistical tests. Simultaneous power simulations are conducted using \link[SimulPower]{pwrMultivars} or \link[SimulPower]{pwrMultimodels}. A detailed walkthrough and set of vignettes for this and other SimulPower functions is available [here](https://doi.org/10.31219/osf.io/w96uk).
#'
#' The pwrCurve function is the first step in the suggested SimulPower workflow. Start here to visualize the approximate simultaneous power space occupied by your set of tests, then use either \link[SimulPower]{pwrMultivars} or \link[SimulPower]{pwrMultimodels} with lrager numbers of iterations for final power calculations with more stable estimates.
#'
#' When you use this function (and we hope you do!), please cite the package:
#'
#' Le Forestier, J. M. (2020). SimulPower: Simultaneous power analysis for a set of statistical tests. https://doi.org/10.31219/osf.io/w96uk
#'
#' and/or cite the accompanying paper:
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. L. (Forthcoming). Statistical power for a set of tests.
#'
#' @usage pwrCurve(method = NULL, min = NULL, max = NULL,
#' increment = 20, thresholds = c(80, 90, 95), es_units = NULL,
#' es1 = NULL, es2 = NULL, es3...es10 = 0,
#' iv1iv2_cov...iv9iv10_cov = 0,
#' null_effect = 0, iterations = 1000, alpha = .05,
#' bonferroni = FALSE, seed = 1, iv1iv2_cov...iv9iv10_cov = 0)
#'
#' @param method Specify which SimulPower function you would like to use to computer simultaneous power. Accepts either "pwrMultivars" or "pwrMultimodels". This argument has no default.
#' @param min Set the minimum sample size for which you would like to estimate simultaneous power. Accepts any positive number. This argument has no default.
#' @param max Set the maximum sample size for which you would like to estimate simultaneous power. Accepts any positive number. This argument has no default.
#' @param increment Set the increments by which you would like the sample size to increase between power analyses. Accepts any positive whole number. Default = 20.
#' @param thresholds Set the power thresholds for which you would like sample size estimates printed. Accepts either any single number between 0 and 100 or a vector of numbers between 0 and 100. Default = c(80, 90, 95).
#' @param es_units Set the units in which you are specifying your effect sizes. Accepts "d" for Cohen's d, "r" for correlation coefficients, and "r2" for percent of variance accounted for. This argument has no default.
#' @param null_effect For which, if any, of your predictors or models are you computing "null power?" If you want to compute "power" to NOT detect an effect, use this argument to specify which effects are predicted nulls by setting this argument equal to the number(s) corresponding to the predictors you hypothesize to be null. Accepts either a single whole number between 1 and the number of predictors you have specified or a vector of numbers between 1 and the the number of predictors you have specified. Default = no null effects.
#' @param iterations How many times you would like to run your model in random samples drawn from your population? One model will be run in each random sample. Accepts any whole number greater than 0. Default = 1000.
#' @param alpha Set your alpha level. This is the threshold below which p-values will be considered significant. Accepts any number greater than 0 and less than 1. Default = 0.05.
#' @param bonferroni Apply a bonferroni correction? This is suggested if you intend on interpreting the results of multiple tests individually, but not if you intend on assessing a single research question by triangulating across multiple tests (Le Forestier, Page-Gould, & Chasteen, Forthcoming). Accepts TRUE or FALSE. Default = FALSE.
#' @param seed Set a seed to make your results reproducible. Accepts any number. Default = 1.
#' @param es1...es10 The effect size, expressed in units specified in the es_units argument, for the relationship between each predictor and the dependent variable.You should always specify these in order, beginning with es1, and not skipping any. Accepts any number. These arguments have no defaults.
#' @param iv1iv2_cov...iv9iv10_cov The relationships between each set of predictors, specified in correlation coefficients. Specifying relationships between predictors is optional unless your predictors, together, account for more than 100% of the variance in your DV, in which case you must specify relationships between your predictors to make that possible. Accepts any number between -1 and 1. Default = 0.
#' @param dv1dv2_cov...dv9dv10_cov The relationships between each set of dependent variables, specified in correlation coefficients. Specifying relationships between DVs is optional. Accepts any number between -1 and 1. Default = 0.
#' @param iv1dv2_cov...iv9dv10_cov The relationships between each set of predictors and dependent variables from separate models, specified in correlation coefficients. Specifying relationships between separate models' predictors and DVs is optional. Accepts any number between -1 and 1. Default = 0.
#'
#' @return A plot of power curves and a set of sample size estimates for default or user-defined power thresholds.
#'
#' A dataframe containing a power estiamte, expressed as a decimal, for each of the effects individually, and for all the effects simultaneously.
#'
#' @author Joel Le Forestier (joel.leforestier@@mail.utoronto.ca)
#'
#' @references Le Forestier, J. M., Page-Gould, E., & Chasteen, A. (Forthcoming). Statistical power for a set of tests.
#'
#' @examples # An example using pwrMultivars
#'
#' pwrCurve(method = "pwrMultivars", min = 200, max = 500,
#' es_units = "d", es1 = .35, es2 = .45, es3 = .50)
#'
#' # An example using pwrMultimodels
#'
#' pwrCurve(method = "pwrMultimodels", min = 100, max = 400, increment = 40,
#' es_units = "r2", es1 = .05, es2 = .06, es3 = .08, es4 = .06, es5 = .08, es6 = .07, es7 = .04)
#'
#' @export

pwrCurve <- function(method, min, max, increment = 20, thresholds = c(80, 90, 95),
                              es_units, es1, es2,
                              es3 = 0, es4 = 0, es5 = 0, es6 = 0, es7 = 0, es8 = 0, es9 = 0, es10 = 0,
                              iv1iv2_cov = 0, iv1iv3_cov = 0, iv1iv4_cov = 0, iv1iv5_cov = 0, iv1iv6_cov = 0, iv1iv7_cov = 0, iv1iv8_cov = 0, iv1iv9_cov = 0, iv1iv10_cov = 0,
                              iv2iv3_cov = 0, iv2iv4_cov = 0, iv2iv5_cov = 0, iv2iv6_cov = 0, iv2iv7_cov = 0, iv2iv8_cov = 0, iv2iv9_cov = 0, iv2iv10_cov = 0,
                              iv3iv4_cov = 0, iv3iv5_cov = 0, iv3iv6_cov = 0, iv3iv7_cov = 0, iv3iv8_cov = 0, iv3iv9_cov = 0, iv3iv10_cov = 0,
                              iv4iv5_cov = 0, iv4iv6_cov = 0, iv4iv7_cov = 0, iv4iv8_cov = 0, iv4iv9_cov = 0, iv4iv10_cov = 0,
                              iv5iv6_cov = 0, iv5iv7_cov = 0, iv5iv8_cov = 0, iv5iv9_cov = 0, iv5iv10_cov = 0,
                              iv6iv7_cov = 0, iv6iv8_cov = 0, iv6iv9_cov = 0, iv6iv10_cov = 0,
                              iv7iv8_cov = 0, iv7iv9_cov = 0, iv7iv10_cov = 0,
                              iv8iv9_cov = 0, iv8iv10_cov = 0,
                              iv9iv10_cov = 0,

                              dv1dv2_cov = 0, dv1dv3_cov = 0, dv1dv4_cov = 0, dv1dv5_cov = 0, dv1dv6_cov = 0, dv1dv7_cov = 0, dv1dv8_cov = 0, dv1dv9_cov = 0, dv1dv10_cov = 0,
                              dv2dv3_cov = 0, dv2dv4_cov = 0, dv2dv5_cov = 0, dv2dv6_cov = 0, dv2dv7_cov = 0, dv2dv8_cov = 0, dv2dv9_cov = 0, dv2dv10_cov = 0,
                              dv3dv4_cov = 0, dv3dv5_cov = 0, dv3dv6_cov = 0, dv3dv7_cov = 0, dv3dv8_cov = 0, dv3dv9_cov = 0, dv3dv10_cov = 0,
                              dv4dv5_cov = 0, dv4dv6_cov = 0, dv4dv7_cov = 0, dv4dv8_cov = 0, dv4dv9_cov = 0, dv4dv10_cov = 0,
                              dv5dv6_cov = 0, dv5dv7_cov = 0, dv5dv8_cov = 0, dv5dv9_cov = 0, dv5dv10_cov = 0,
                              dv6dv7_cov = 0, dv6dv8_cov = 0, dv6dv9_cov = 0, dv6dv10_cov = 0,
                              dv7dv8_cov = 0, dv7dv9_cov = 0, dv7dv10_cov = 0,
                              dv8dv9_cov = 0, dv8dv10_cov = 0,
                              dv9dv10_cov = 0,

                              iv1dv2_cov = 0, iv1dv3_cov = 0, iv1dv4_cov = 0, iv1dv5_cov = 0, iv1dv6_cov = 0, iv1dv7_cov = 0, iv1dv8_cov = 0, iv1dv9_cov = 0, iv1dv10_cov = 0,
                              iv2dv1_cov = 0, iv3dv1_cov = 0, iv4dv1_cov = 0, iv5dv1_cov = 0, iv6dv1_cov = 0, iv7dv1_cov = 0, iv8dv1_cov = 0, iv9dv1_cov = 0, iv10dv1_cov = 0,
                              iv2dv3_cov = 0, iv2dv4_cov = 0, iv2dv5_cov = 0, iv2dv6_cov = 0, iv2dv7_cov = 0, iv2dv8_cov = 0, iv2dv9_cov = 0, iv2dv10_cov = 0,
                              iv3dv2_cov = 0, iv4dv2_cov = 0, iv5dv2_cov = 0, iv6dv2_cov = 0, iv7dv2_cov = 0, iv8dv2_cov = 0, iv9dv2_cov = 0, iv10dv2_cov = 0,
                              iv3dv4_cov = 0, iv3dv5_cov = 0, iv3dv6_cov = 0, iv3dv7_cov = 0, iv3dv8_cov = 0, iv3dv9_cov = 0, iv3dv10_cov = 0,
                              iv4dv3_cov = 0, iv5dv3_cov = 0, iv6dv3_cov = 0, iv7dv3_cov = 0, iv8dv3_cov = 0, iv9dv3_cov = 0, iv10dv3_cov = 0,
                              iv4dv5_cov = 0, iv4dv6_cov = 0, iv4dv7_cov = 0, iv4dv8_cov = 0, iv4dv9_cov = 0, iv4dv10_cov = 0,
                              iv5dv4_cov = 0, iv6dv4_cov = 0, iv7dv4_cov = 0, iv8dv4_cov = 0, iv9dv4_cov = 0, id10dv4_cov = 0,
                              iv5dv6_cov = 0, iv5dv7_cov = 0, iv5dv8_cov = 0, iv5dv9_cov = 0, iv5dv10_cov = 0,
                              iv6dv5_cov = 0, iv7dv5_cov = 0, iv8dv5_cov = 0, iv9dv5_cov = 0, iv10dv5_cov = 0,
                              iv6dv7_cov = 0, iv6dv8_cov = 0, iv6dv9_cov = 0, iv6dv10_cov = 0,
                              iv7dv6_cov = 0, iv8dv6_cov = 0, iv9dv6_cov = 0, iv10dv6_cov = 0,
                              iv7dv8_cov = 0, iv7dv9_cov = 0, iv7dv10_cov = 0,
                              iv8dv7_cov = 0, iv9dv7_cov = 0, iv10dv7_cov = 0,
                              iv8dv9_cov = 0, iv8dv10_cov = 0,
                              iv9dv8_cov = 0, iv10dv8_cov = 0,
                              iv9dv10_cov = 0,
                              iv10dv9_cov = 0,

                              null_effect = 0, iterations = 1000, alpha = .05, bonferroni = FALSE, seed = 1) {

  # Throw a warning if method is wrong or missing #
  if(method == "pwrmultivars" |
     method == "pwrmultiVars" |
     method == "PwrMultivars" |
     method == "PwrmultiVars" |
     method == "pwrMultiVars" |
     method == "PwrMultiVars" |
     method == "multivars" |
     method == "Multivars" |
     method == "multiVars" |
     method == "MultiVars" |
     method == "vars" |
     method == "Vars" |
     method == "pwrvars" |
     method == "Pwrvars" |
     method == "pwrVars" |
     method == "PwrVars") {
    method <- "pwrMultivars"
  }

  if(method == "pwrmultimodels" |
     method == "pwrmultiModels" |
     method == "PwrMultimodels" |
     method == "PwrmultiModels" |
     method == "pwrMultiModels" |
     method == "PwrMultiModels" |
     method == "multimodels" |
     method == "Multimodels" |
     method == "multiModels" |
     method == "MultiModels" |
     method == "models" |
     method == "Models" |
     method == "pwrmodels" |
     method == "Pwrmodels" |
     method == "pwrModels" |
     method == "PwrModels") {
    method <- "pwrMultimodels"
  }

  if(method != "pwrMultivars" & method != "pwrMultimodels") {
    stop("You have not correctly specified which SimulPower function you would like to use to estimate simutaneous power.")
  }

  # Throw a warning if the increment is 0, negative, or not a whole number #
  if(increment < 1 | round(increment, 0) != increment) {
    stop("You have specified an invalid increment size. Please specify a whole number greater than 0.")
  }

  # Throw a warning if min or max sample size is 0, negative, or not a whole number #
  if(min < 1 | round(min, 0) != min | max < 1 | round(max, 0) != max) {
    stop("You have specified an invalid sample size. Please specify a whole number greater than 0.")
  }

  # Throw a warning if any of the thresholds are > 100 or < 0 #
  error_high <- thresholds > 100
  error_low <- thresholds < 0

  if(TRUE %in% error_high | TRUE %in% error_low) {
    stop("You have specified an impossible power threshold in the \"thresholds\" aregument. Please specify a number, or vector of numbers, between 0 and 1.")
  }

  # Select a function and perform the power analyses #
  if(method == "pwrMultivars") {

    data <- data.frame()

    # let them know it's working
    message(paste("Performing", floor((max - min) / increment + 1), "simultaneous power analyses with", iterations, "sets of tests in each. This may take a few minutes."))

    # Create predictors variable #
    dummy_beta <- 0

    if (1 %in% null_effect) {
      tally_es1 <- "null"
    } else {tally_es1 <- es1}

    if (2 %in% null_effect) {
      tally_es2 <- "null"
    } else {tally_es2 <- es2}

    if (3 %in% null_effect) {
      tally_es3 <- "null"
    } else {tally_es3 <- es3}

    if (4 %in% null_effect) {
      tally_es4 <- "null"
    } else {tally_es4 <- es4}

    if (5 %in% null_effect) {
      tally_es5 <- "null"
    } else {tally_es5 <- es5}

    if (6 %in% null_effect) {
      tally_es6 <- "null"
    } else {tally_es6 <- es6}

    if (7 %in% null_effect) {
      tally_es7 <- "null"
    } else {tally_es7 <- es7}

    if (8 %in% null_effect) {
      tally_es8 <- "null"
    } else {tally_es8 <- es8}

    if (9 %in% null_effect) {
      tally_es9 <- "null"
    } else {tally_es9 <- es9}

    if (10 %in% null_effect) {
      tally_es10 <- "null"
    } else {tally_es10 <- es10}

    specified_params <- table(c(tally_es1,
                                tally_es2,
                                tally_es3,
                                tally_es4,
                                tally_es5,
                                tally_es6,
                                tally_es7,
                                tally_es8,
                                tally_es9,
                                tally_es10,
                                dummy_beta))

    predictors <- 11 - specified_params[["0"]]

    # run the power analyses #
    for (n in seq(from=min, to=max, by=increment)) {
      data <- rbind(data,
                    suppressMessages(SimulPower::pwrMultivars(
                      n = n,
                      es_units = es_units,
                      null_effect = null_effect,
                      iterations = iterations,
                      alpha = alpha,
                      seed = seed,
                      es1 = es1, es2 = es2, es3 = es3, es4 = es4, es5 = es5, es6 = es6, es7 = es7, es8 = es8, es9 = es9, es10 = es10,
                      iv1iv2_cov = iv1iv2_cov, iv1iv3_cov = iv1iv3_cov, iv1iv4_cov = iv1iv4_cov, iv1iv5_cov = iv1iv5_cov, iv1iv6_cov = iv1iv6_cov, iv1iv7_cov = iv1iv7_cov, iv1iv8_cov = iv1iv8_cov, iv1iv9_cov = iv1iv9_cov, iv1iv10_cov = iv1iv10_cov,
                      iv2iv3_cov = iv2iv3_cov, iv2iv4_cov = iv2iv4_cov, iv2iv5_cov = iv2iv5_cov, iv2iv6_cov = iv2iv6_cov, iv2iv7_cov = iv2iv7_cov, iv2iv8_cov = iv2iv8_cov, iv2iv9_cov = iv2iv9_cov, iv2iv10_cov = iv2iv10_cov,
                      iv3iv4_cov = iv3iv4_cov, iv3iv5_cov = iv3iv5_cov, iv3iv6_cov = iv3iv6_cov, iv3iv7_cov = iv3iv7_cov, iv3iv8_cov = iv3iv8_cov, iv3iv9_cov = iv3iv9_cov, iv3iv10_cov = iv3iv10_cov,
                      iv4iv5_cov = iv4iv5_cov, iv4iv6_cov = iv4iv6_cov, iv4iv7_cov = iv4iv7_cov, iv4iv8_cov = iv4iv8_cov, iv4iv9_cov = iv4iv9_cov, iv4iv10_cov = iv4iv10_cov,
                      iv5iv6_cov = iv5iv6_cov, iv5iv7_cov = iv5iv7_cov, iv5iv8_cov = iv5iv8_cov, iv5iv9_cov = iv5iv9_cov, iv5iv10_cov = iv5iv10_cov,
                      iv6iv7_cov = iv6iv7_cov, iv6iv8_cov = iv6iv8_cov, iv6iv9_cov = iv6iv9_cov, iv6iv10_cov = iv6iv10_cov,
                      iv7iv8_cov = iv7iv8_cov, iv7iv9_cov = iv7iv9_cov, iv7iv10_cov = iv7iv10_cov,
                      iv8iv9_cov = iv8iv9_cov, iv8iv10_cov = iv8iv10_cov,
                      iv9iv10_cov = iv9iv10_cov,
                      print_result = FALSE)))
    }

    data$n <- seq(from = min, to = max, by = increment)

    # plot the result #
    color_options <- c("#9F018A", #color blind-friendly palette
                       "#009F81",
                       "#FF5AAF",
                       "#00FCCF",
                       "#8400CD",
                       "#008DF9",
                       "#A40122",
                       "#E20134",
                       "#FF6E3A",
                       "#FFC33B")
    palette <- append(color_options[1:predictors], "#000000") #simultaneous power is always black

    longdata <- tidyr::pivot_longer(data = data,
                                    cols = -n,
                                    names_to = "param",
                                    values_to = "estimate")
    names(longdata) <- c("n", "parameter", "power")

    for(t in 1:predictors) {
      longdata$parameter[longdata$parameter == paste0("es", t, "_power")] <- paste("predictor", t)
    }
    longdata$parameter[longdata$parameter == "simultaneous_power"] <- "simultaneous power"

    plot <- ggplot2::ggplot(data = longdata, mapping = ggplot2::aes(x = n, y = power)) +
      ggplot2::scale_color_manual(name = "Parameter",
                                  values = palette) +
      ggplot2::geom_line(ggplot2::aes(color = parameter), size = 1) +
      ggplot2::scale_y_continuous(limits = c(0, 100),
                         breaks=c(0, 20, 40, 60, 80, 100),
                         labels=c("0%", "20%", "40%", "60%", "80%", "100%")) +
      ggplot2::scale_x_continuous(limits = c(min, max),
                                  breaks = seq(min, max, by = increment)) +
      ggplot2::xlab("Sample Size") +
      ggplot2::ylab("Power") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  } else if (method == "pwrMultimodels") {

    data <- data.frame()

    # let them know it's working
    message(paste("Performing", floor((max - min) / increment + 1), "simultaneous power analyses with", iterations, "sets of tests in each. This may take a few minutes."))

    # Create models variable #
    dummy_es <- 0

    if (1 %in% null_effect) {
      tally_es1 <- "null"
    } else {tally_es1 <- es1}

    if (2 %in% null_effect) {
      tally_es2 <- "null"
    } else {tally_es2 <- es2}

    if (3 %in% null_effect) {
      tally_es3 <- "null"
    } else {tally_es3 <- es3}

    if (4 %in% null_effect) {
      tally_es4 <- "null"
    } else {tally_es4 <- es4}

    if (5 %in% null_effect) {
      tally_es5 <- "null"
    } else {tally_es5 <- es5}

    if (6 %in% null_effect) {
      tally_es6 <- "null"
    } else {tally_es6 <- es6}

    if (7 %in% null_effect) {
      tally_es7 <- "null"
    } else {tally_es7 <- es7}

    if (8 %in% null_effect) {
      tally_es8 <- "null"
    } else {tally_es8 <- es8}

    if (9 %in% null_effect) {
      tally_es9 <- "null"
    } else {tally_es9 <- es9}

    if (10 %in% null_effect) {
      tally_es10 <- "null"
    } else {tally_es10 <- es10}

    specified_models <- table(c(tally_es1,
                                tally_es2,
                                tally_es3,
                                tally_es4,
                                tally_es5,
                                tally_es6,
                                tally_es7,
                                tally_es8,
                                tally_es9,
                                tally_es10,
                                dummy_es))

    models <- 11 - specified_models[["0"]]

    # run the power analyses #
    for (n in seq(from=min, to=max, by=increment)) {
      data <- rbind(data,
                    suppressMessages(SimulPower::pwrMultimodels(
                      n = n,
                      es_units = es_units,
                      null_effect = null_effect,
                      iterations = iterations,
                      alpha = alpha,
                      seed = seed,
                      es1 = es1, es2 = es2, es3 = es3, es4 = es4, es5 = es5, es6 = es6, es7 = es7, es8 = es8, es9 = es9, es10 = es10,
                      iv1iv2_cov = iv1iv2_cov, iv1iv3_cov = iv1iv3_cov, iv1iv4_cov = iv1iv4_cov, iv1iv5_cov = iv1iv5_cov, iv1iv6_cov = iv1iv6_cov, iv1iv7_cov = iv1iv7_cov, iv1iv8_cov = iv1iv8_cov, iv1iv9_cov = iv1iv9_cov, iv1iv10_cov = iv1iv10_cov,
                      iv2iv3_cov = iv2iv3_cov, iv2iv4_cov = iv2iv4_cov, iv2iv5_cov = iv2iv5_cov, iv2iv6_cov = iv2iv6_cov, iv2iv7_cov = iv2iv7_cov, iv2iv8_cov = iv2iv8_cov, iv2iv9_cov = iv2iv9_cov, iv2iv10_cov = iv2iv10_cov,
                      iv3iv4_cov = iv3iv4_cov, iv3iv5_cov = iv3iv5_cov, iv3iv6_cov = iv3iv6_cov, iv3iv7_cov = iv3iv7_cov, iv3iv8_cov = iv3iv8_cov, iv3iv9_cov = iv3iv9_cov, iv3iv10_cov = iv3iv10_cov,
                      iv4iv5_cov = iv4iv5_cov, iv4iv6_cov = iv4iv6_cov, iv4iv7_cov = iv4iv7_cov, iv4iv8_cov = iv4iv8_cov, iv4iv9_cov = iv4iv9_cov, iv4iv10_cov = iv4iv10_cov,
                      iv5iv6_cov = iv5iv6_cov, iv5iv7_cov = iv5iv7_cov, iv5iv8_cov = iv5iv8_cov, iv5iv9_cov = iv5iv9_cov, iv5iv10_cov = iv5iv10_cov,
                      iv6iv7_cov = iv6iv7_cov, iv6iv8_cov = iv6iv8_cov, iv6iv9_cov = iv6iv9_cov, iv6iv10_cov = iv6iv10_cov,
                      iv7iv8_cov = iv7iv8_cov, iv7iv9_cov = iv7iv9_cov, iv7iv10_cov = iv7iv10_cov,
                      iv8iv9_cov = iv8iv9_cov, iv8iv10_cov = iv8iv10_cov,
                      iv9iv10_cov = iv9iv10_cov,
                      dv1dv2_cov = dv1dv2_cov, dv1dv3_cov = dv1dv3_cov, dv1dv4_cov = dv1dv4_cov, dv1dv5_cov = dv1dv5_cov, dv1dv6_cov = dv1dv6_cov, dv1dv7_cov = dv1dv7_cov, dv1dv8_cov = dv1dv8_cov, dv1dv9_cov = dv1dv9_cov, dv1dv10_cov = dv1dv10_cov,
                      dv2dv3_cov = dv2dv3_cov, dv2dv4_cov = dv2dv4_cov, dv2dv5_cov = dv2dv5_cov, dv2dv6_cov = dv2dv6_cov, dv2dv7_cov = dv2dv7_cov, dv2dv8_cov = dv2dv8_cov, dv2dv9_cov = dv2dv9_cov, dv2dv10_cov = dv2dv10_cov,
                      dv3dv4_cov = dv3dv4_cov, dv3dv5_cov = dv3dv5_cov, dv3dv6_cov = dv3dv6_cov, dv3dv7_cov = dv3dv7_cov, dv3dv8_cov = dv3dv8_cov, dv3dv9_cov = dv3dv9_cov, dv3dv10_cov = dv3dv10_cov,
                      dv4dv5_cov = dv4dv5_cov, dv4dv6_cov = dv4dv6_cov, dv4dv7_cov = dv4dv7_cov, dv4dv8_cov = dv4dv8_cov, dv4dv9_cov = dv4dv9_cov, dv4dv10_cov = dv4dv10_cov,
                      dv5dv6_cov = dv5dv6_cov, dv5dv7_cov = dv5dv7_cov, dv5dv8_cov = dv5dv8_cov, dv5dv9_cov = dv5dv9_cov, dv5dv10_cov = dv5dv10_cov,
                      dv6dv7_cov = dv6dv7_cov, dv6dv8_cov = dv6dv8_cov, dv6dv9_cov = dv6dv9_cov, dv6dv10_cov = dv6dv10_cov,
                      dv7dv8_cov = dv7dv8_cov, dv7dv9_cov = dv7dv9_cov, dv7dv10_cov = dv7dv10_cov,
                      dv8dv9_cov = dv8dv9_cov, dv8dv10_cov = dv8dv10_cov,
                      dv9dv10_cov = dv9dv10_cov,
                      iv1dv2_cov = iv1dv2_cov, iv1dv3_cov = iv1dv3_cov, iv1dv4_cov = iv1dv4_cov, iv1dv5_cov = iv1dv5_cov, iv1dv6_cov = iv1dv6_cov, iv1dv7_cov = iv1dv7_cov, iv1dv8_cov = iv1dv8_cov, iv1dv9_cov = iv1dv9_cov, iv1dv10_cov = iv1dv10_cov,
                      iv2dv1_cov = iv2dv1_cov, iv3dv1_cov = iv3dv1_cov, iv4dv1_cov = iv4dv1_cov, iv5dv1_cov = iv5dv1_cov, iv6dv1_cov = iv6dv1_cov, iv7dv1_cov = iv7dv1_cov, iv8dv1_cov = iv8dv1_cov, iv9dv1_cov = iv9dv1_cov, iv10dv1_cov = iv10dv1_cov,
                      iv2dv3_cov = iv2dv3_cov, iv2dv4_cov = iv2dv4_cov, iv2dv5_cov = iv2dv5_cov, iv2dv6_cov = iv2dv6_cov, iv2dv7_cov = iv2dv7_cov, iv2dv8_cov = iv2dv8_cov, iv2dv9_cov = iv2dv9_cov, iv2dv10_cov = iv2dv10_cov,
                      iv3dv2_cov = iv3dv2_cov, iv4dv2_cov = iv4dv2_cov, iv5dv2_cov = iv5dv2_cov, iv6dv2_cov = iv6dv2_cov, iv7dv2_cov = iv7dv2_cov, iv8dv2_cov = iv8dv2_cov, iv9dv2_cov = iv9dv2_cov, iv10dv2_cov = iv10dv2_cov,
                      iv3dv4_cov = iv3dv4_cov, iv3dv5_cov = iv3dv5_cov, iv3dv6_cov = iv3dv6_cov, iv3dv7_cov = iv3dv7_cov, iv3dv8_cov = iv3dv8_cov, iv3dv9_cov = iv3dv9_cov, iv3dv10_cov = iv3dv10_cov,
                      iv4dv3_cov = iv4dv3_cov, iv5dv3_cov = iv5dv3_cov, iv6dv3_cov = iv6dv3_cov, iv7dv3_cov = iv7dv3_cov, iv8dv3_cov = iv8dv3_cov, iv9dv3_cov = iv9dv3_cov, iv10dv3_cov = iv10dv3_cov,
                      iv4dv5_cov = iv4dv5_cov, iv4dv6_cov = iv4dv6_cov, iv4dv7_cov = iv4dv7_cov, iv4dv8_cov = iv4dv8_cov, iv4dv9_cov = iv4dv9_cov, iv4dv10_cov = iv4dv10_cov,
                      iv5dv4_cov = iv5dv4_cov, iv6dv4_cov = iv6dv4_cov, iv7dv4_cov = iv7dv4_cov, iv8dv4_cov = iv8dv4_cov, iv9dv4_cov = iv9dv4_cov, id10dv4_cov = id10dv4_cov,
                      iv5dv6_cov = iv5dv6_cov, iv5dv7_cov = iv5dv7_cov, iv5dv8_cov = iv5dv8_cov, iv5dv9_cov = iv5dv9_cov, iv5dv10_cov = iv5dv10_cov,
                      iv6dv5_cov = iv6dv5_cov, iv7dv5_cov = iv7dv5_cov, iv8dv5_cov = iv8dv5_cov, iv9dv5_cov = iv9dv5_cov, iv10dv5_cov = iv10dv5_cov,
                      iv6dv7_cov = iv6dv7_cov, iv6dv8_cov = iv6dv8_cov, iv6dv9_cov = iv6dv9_cov, iv6dv10_cov = iv6dv10_cov,
                      iv7dv6_cov = iv7dv6_cov, iv8dv6_cov = iv8dv6_cov, iv9dv6_cov = iv9dv6_cov, iv10dv6_cov = iv10dv6_cov,
                      iv7dv8_cov = iv7dv8_cov, iv7dv9_cov = iv7dv9_cov, iv7dv10_cov = iv7dv10_cov,
                      iv8dv7_cov = iv8dv7_cov, iv9dv7_cov = iv9dv7_cov, iv10dv7_cov = iv10dv7_cov,
                      iv8dv9_cov = iv8dv9_cov, iv8dv10_cov = iv8dv10_cov,
                      iv9dv8_cov = iv9dv8_cov, iv10dv8_cov = iv10dv8_cov,
                      iv9dv10_cov = iv9dv10_cov,
                      iv10dv9_cov = iv10dv9_cov,
                      print_result = FALSE)))
    }

  data$n <- seq(from = min, to = max, by = increment)

    # plot the result #
    color_options <- c("#9F018A", #color blind-friendly palette
                       "#009F81",
                       "#FF5AAF",
                       "#00FCCF",
                       "#8400CD",
                       "#008DF9",
                       "#A40122",
                       "#E20134",
                       "#FF6E3A",
                       "#FFC33B")
    palette <- append(color_options[1:models], "#000000") #simultaneous power is always black

    longdata <- tidyr::pivot_longer(data = data,
                                    cols = -n,
                                    names_to = "param",
                                    values_to = "estimate")
    names(longdata) <- c("n", "model", "power")

    for(t in 1:models) {
      longdata$model[longdata$model == paste0("model", t, "_power")] <- paste("model", t)
    }
    longdata$model[longdata$model == "simultaneous_power"] <- "simultaneous power"

    plot <- ggplot2::ggplot(data = longdata, mapping = ggplot2::aes(x = n, y = power)) +
      ggplot2::scale_color_manual(name = "Model",
                                  values = palette) +
      ggplot2::geom_line(ggplot2::aes(color = model), size = 1) +
      ggplot2::scale_y_continuous(limits = c(0, 100),
                                  breaks=c(0, 20, 40, 60, 80, 100),
                                  labels=c("0%", "20%", "40%", "60%", "80%", "100%")) +
      ggplot2::scale_x_continuous(limits = c(min, max),
                                   breaks = seq(min, max, by = increment)) +
      ggplot2::xlab("Sample Size") +
      ggplot2::ylab("Power") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  # print the sample size estimates #
  thresh_n <- list()
  for (x in 1:length(thresholds)) {
    if (thresholds[x] < min(data$simultaneous_power)) {
      thresh_n <- append(thresh_n, NA)
      thresh_n <- append(thresh_n, NA)
    } else {
      thresh_n <- append(thresh_n, data$n[which(data$simultaneous_power >= thresholds[x])[1] - 1])
      thresh_n <- append(thresh_n, data$n[which(data$simultaneous_power >= thresholds[x])[1]])
    }
    if (is.na(thresh_n[[x * 2 - 1]]) | is.na(thresh_n[[x * 2]])) {
      cat("The ", thresholds[x], "% simultaneous power threshold was not crossed.", "\n", sep = "")
    } else {
      cat("The ", thresholds[x], "% simultaneous power threshold was crossed between n = ", thresh_n[[x * 2 - 1]], " and n = ", thresh_n[[x * 2]], "\n", sep = "")
    }
  }

  # return the plot #
  return(plot)
}



