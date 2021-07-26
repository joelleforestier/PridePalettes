#' Simulate simultaneous power for multiple tests in separate models
#'
#' Simulate power to simultaneously detect predicted effects for a set of statistical tests in separate bivariate models. This function simulates data based on a correlation matrix imposed using the mvrnorm function from the MASS package (Venables & Ripley, 2002) and can be used to estimate power for up to 10 bivariate models with single independent variables and single dependent variables. A detailed walkthrough and set of vignettes for this and other SimulPower functions is available [here](https://doi.org/10.31219/osf.io/w96uk).
#'
#' When you use this function (and we hope you do!), please cite the package:
#'
#' Le Forestier, J. M. (2020). SimulPower: Simultaneous power analysis for a set of statistical tests. https://doi.org/10.31219/osf.io/w96uk
#'
#' and/or cite the accompanying paper:
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. L. (Forthcoming). Statistical power for a set of tests.
#'
#' @usage pwrMultimodels(n = NULL, es_units = NULL,
#' es1 = NULL, es2 = NULL, es3...es10 = 0,
#' null_effect = 0, iterations = 5000,
#' alpha = .05, bonferroni = FALSE, seed = 1,
#' iv1iv2_cov...iv9iv10_cov = 0)
#'
#' @param n Set the size of each sample to be drawn from the population. This is the sample size for which you are estimating statistical power. In other words, setting n to equal 100 will estimate statistical power at n = 100. Accepts any positive number. This argument has no default.
#' @param es_units Set the units in which you are specifying your effect sizes. Accepts "d" for Cohen's d, "r" for correlation coefficients, and "r2" for percent of variance accounted for. This argument has no default.
#' @param null_effect For which, if any, of your models are you computing "null power?" If you want to compute "power" to NOT detect an effect, use this argument to specify which models are predicted nulls by setting this argument equal to the number(s) corresponding to the models you hypothesize to be null. Accepts either a single whole number between 1 and the number of models you have specified or a vector of numbers between 1 and the the number of models you have specified. Default = no null effects.
#' @param iterations How many times you would like to run your models in random samples drawn from your population? One model will be run in each random sample. Accepts any whole number greater than 0. Default = 5000.
#' @param alpha Set your alpha level. This is the threshold below which p-values will be considered significant. Accepts any number greater than 0 and less than 1. Default = 0.05.
#' @param bonferroni Apply a bonferroni correction? This is suggested if you intend on interpreting the results of multiple tests individually, but not if you intend on assessing a single research question by triangulating across multiple tests (Le Forestier, Page-Gould, & Chasteen, Forthcoming). Accepts TRUE or FALSE. Default = FALSE.
#' @param seed Set a seed to make your results reproducible. Accepts any number. Default = 1.
#' @param es1...es10 The effect size, expressed in units specified in the es_units argument, for each model. You should always specify these in order, beginning with es1, and not skipping any. Accepts any number. These arguments have no defaults.
#' @param iv1iv2_cov...iv9iv10_cov The relationships between each set of predictors, specified in correlation coefficients. Specifying relationships between predictors is optional. Accepts any number between -1 and 1. Default = 0.
#' @param dv1dv2_cov...dv9dv10_cov The relationships between each set of dependent variables, specified in correlation coefficients. Specifying relationships between DVs is optional. Accepts any number between -1 and 1. Default = 0.
#' @param iv1dv2_cov...iv9dv10_cov The relationships between each set of predictors and dependent variables from separate models, specified in correlation coefficients. Specifying relationships between separate models' predictors and DVs is optional. Accepts any number between -1 and 1. Default = 0.
#' @param print_result Should power analysis results be printed to the console? Accepts TRUE or FALSE. Default = TRUE.
#'
#' @return A dataframe containing a power estiamte, expressed as a decimal, for each of the models individually, and for all the models simultaneously.
#'
#' @author Joel Le Forestier (joel.leforestier@@mail.utoronto.ca)
#'
#' @references Le Forestier, J. M., Page-Gould, E., & Chasteen, A. (Forthcoming). Statistical power for a set of tests.
#'
#' Venables, W. N. & Ripley, B. D. (2002). Modern applied statistics with S. Springer.
#'
#'
#' @examples # A basic example, leaving all the defaults in place.
#'
#' pwrMultimodels(n = 150, es_units = "r", es1 = .2, es2 = .45)
#'
#' # Another example, customizing additional parameters.
#'
#' pwrMultimodels(n = 300, es_units = "r2", es1 = .04, es2 = .01, es3 = .00,
#'      null_effect = 3, alpha = .01, seed = 123)
#'
#' @export

pwrMultimodels <- function(n, es_units, es1, es2,
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
                          null_effect = 0, iterations = 5000, alpha = .05, bonferroni = FALSE, seed = 1, print_result = TRUE) {

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

  # Throw a warning if the user has specified the predictors out of order #
  `%notin%` <- Negate(`%in%`)

  if (
    ((es4 != 0 | 4 %in% null_effect) &
     (es3 == 0 & 3 %notin% null_effect)) |

    ((es5 != 0 | 5 %in% null_effect) &
     ((es3 == 0 & 3 %notin% null_effect) |
      (es4 == 0 & 4 %notin% null_effect))) |

    ((es6 != 0 | 6 %in% null_effect) &
     ((es3 == 0 & 3 %notin% null_effect) |
      (es4 == 0 & 4 %notin% null_effect) |
      (es5 == 0 & 5 %notin% null_effect))) |

    ((es7 != 0 | 7 %in% null_effect) &
     ((es3 == 0 & 3 %notin% null_effect) |
      (es4 == 0 & 4 %notin% null_effect) |
      (es5 == 0 & 5 %notin% null_effect) |
      (es6 == 0 & 6 %notin% null_effect))) |

    ((es8 != 0 | 8 %in% null_effect) &
     ((es3 == 0 & 3 %notin% null_effect) |
      (es4 == 0 & 4 %notin% null_effect) |
      (es5 == 0 & 5 %notin% null_effect) |
      (es6 == 0 & 6 %notin% null_effect) |
      (es7 == 0 & 7 %notin% null_effect))) |

    ((es9 != 0 | 9 %in% null_effect) &
     ((es3 == 0 & 3 %notin% null_effect) |
      (es4 == 0 & 4 %notin% null_effect) |
      (es5 == 0 & 5 %notin% null_effect) |
      (es6 == 0 & 6 %notin% null_effect) |
      (es7 == 0 & 7 %notin% null_effect) |
      (es8 == 0 & 8 %notin% null_effect))) |

    ((es10 != 0 | 9 %in% null_effect) &
     ((es3 == 0 & 3 %notin% null_effect) |
      (es4 == 0 & 4 %notin% null_effect) |
      (es5 == 0 & 5 %notin% null_effect) |
      (es6 == 0 & 6 %notin% null_effect) |
      (es7 == 0 & 7 %notin% null_effect) |
      (es8 == 0 & 8 %notin% null_effect) |
      (es9 == 0 & 9 %notin% null_effect)))
  ) {
    stop("Please be sure to specify the predictors in order. For example, if you aim to specify four effect sizes, you should assign them to es1, es2, es3, and es4.")
  }

  # Throw a warning if the number of iterations is 0, negative, or not a whole number #
  if(iterations < 1 | round(iterations, 0) != iterations) {
    stop("You have specified an invalid number of iterations. Please specify a whole number greater than 0.")
  }

  # Throw a warning if the sample size is 0, negative, or not a whole number #
  if(n < 1 | round(n, 0) != n) {
    stop("You have specified an invalid sample size. Please specify a whole number greater than 0.")
  }

  # Correct effect sizes and throw a warning if es_units is not r, r2, or d #
  if(es_units == "D") {
    es_units <- "d"
  } else if (es_units == "R") {
    es_units <- "r"
  } else if (es_units == "R2") {
    es_units <- "r2"
  }

  if(es_units != "d" & es_units != "r" & es_units != "r2") {
    stop("You have specified an invalid type effect size. Please select from \"r\", \"d\", or \"r2\".")
  }

  # Throw a warning if the alpha level is not a number > 0 and < 1 #
  if(alpha <=0 | alpha >= 1) {
    stop("You have specified an incorrect alpha level. Please specify an alpha level greater than 0 and less than 1.")
  }

  # Throw a warning if Bonferroni input isn't T or F #
  if(bonferroni != TRUE & bonferroni != FALSE) {
    stop("You have specified an invalid value for the bonferroni argument. Please specify whether you wish to bonferroni correct your p-values (bonferroni = TRUE) or not (bonferroni = FALSE). ")
  }

  # Throw a warning if print_result input isn't T or F #
  if(print_result != TRUE & print_result != FALSE) {
    stop("You have specified an invalid value for the print_result argument. Please specify whether you wish to print results to the console (print_result = TRUE) or not (print_result = FALSE). ")
  }

  # Generate the correlation matrix for the population #
  cortable <- matrix(c(1, es1, iv1iv2_cov, iv1dv2_cov, iv1iv3_cov, iv1dv3_cov, iv1iv4_cov, iv1dv4_cov, iv1iv5_cov, iv1dv5_cov, iv1iv6_cov, iv1dv6_cov, iv1iv7_cov, iv1dv7_cov, iv1iv8_cov, iv1dv8_cov, iv1iv9_cov, iv1dv9_cov, iv1iv10_cov, iv1dv10_cov,
                       es1, 1, iv2dv1_cov, dv1dv2_cov, iv3dv1_cov, dv1dv3_cov, iv4dv1_cov, dv1dv4_cov, iv5dv1_cov, dv1dv5_cov, iv6dv1_cov, dv1dv6_cov, iv7dv1_cov, dv1dv7_cov, iv8dv1_cov, dv1dv8_cov, iv9dv1_cov, dv1dv9_cov, iv10dv1_cov, dv1dv10_cov,
                       iv1iv2_cov, iv2dv1_cov, 1, es2, iv2iv3_cov, iv2dv3_cov, iv2iv4_cov, iv2dv4_cov, iv2iv5_cov, iv2dv5_cov, iv2iv6_cov, iv2dv6_cov, iv2iv7_cov, iv2dv7_cov, iv2iv8_cov, iv2dv8_cov, iv2iv9_cov, iv2dv9_cov, iv2iv10_cov, iv2dv10_cov,
                       iv1dv2_cov, dv1dv2_cov, es2, 1, iv3dv2_cov, dv2dv3_cov, iv4dv2_cov, dv2dv4_cov, iv5dv2_cov, dv2dv5_cov, iv6dv2_cov, dv2dv6_cov, iv7dv2_cov, dv2dv7_cov, iv8dv2_cov, dv2dv8_cov, iv9dv2_cov, dv2dv9_cov, iv10dv2_cov, dv2dv10_cov,
                       iv1iv3_cov, iv3dv1_cov, iv2iv3_cov, iv3dv2_cov, 1, es3, iv3iv4_cov, iv3dv4_cov, iv3iv5_cov, iv3dv5_cov, iv3iv6_cov, iv3dv6_cov, iv3iv7_cov, iv3dv7_cov, iv3iv8_cov, iv3dv8_cov, iv3iv9_cov, iv3dv9_cov, iv3iv10_cov, iv3dv10_cov,
                       iv1dv3_cov, dv1dv3_cov, iv2dv3_cov, dv2dv3_cov, es3, 1, iv4dv3_cov, dv3dv4_cov, iv5dv3_cov, dv3dv5_cov, iv6dv3_cov, dv3dv6_cov, iv7dv3_cov, dv3dv7_cov, iv8dv3_cov, dv3dv8_cov, iv9dv3_cov, dv3dv9_cov, iv10dv3_cov, dv3dv10_cov,
                       iv1iv4_cov, iv4dv1_cov, iv2iv4_cov, iv4dv2_cov, iv3iv4_cov, iv4dv3_cov, 1, es4, iv4iv5_cov, iv4dv5_cov, iv4iv6_cov, iv4dv6_cov, iv4iv7_cov, iv4dv7_cov, iv4iv8_cov, iv4dv8_cov, iv4iv9_cov, iv4dv9_cov, iv4iv10_cov, iv4dv10_cov,
                       iv1dv4_cov, dv1dv4_cov, iv2dv4_cov, dv2dv4_cov, iv3dv4_cov, dv3dv4_cov, es4, 1, iv5dv4_cov, dv4dv5_cov, iv6dv4_cov, dv4dv6_cov, iv7dv4_cov, dv4dv7_cov, iv8dv4_cov, dv4dv8_cov, iv9dv4_cov, dv4dv9_cov, id10dv4_cov, dv4dv10_cov,
                       iv1iv5_cov, iv5dv1_cov, iv2iv5_cov, iv5dv2_cov, iv3iv5_cov, iv5dv3_cov, iv4iv5_cov, iv5dv4_cov, 1, es5, iv5iv6_cov, iv5dv6_cov, iv5iv7_cov, iv5dv7_cov, iv5iv8_cov, iv5dv8_cov, iv5iv9_cov, iv5dv9_cov, iv5iv10_cov, iv5dv10_cov,
                       iv1dv5_cov, dv1dv5_cov, iv2dv5_cov, dv2dv5_cov, iv3dv5_cov, dv3dv5_cov, iv4dv5_cov, dv4dv5_cov, es5, 1, iv6dv5_cov, dv5dv6_cov, iv7dv5_cov, dv5dv7_cov, iv8dv5_cov, dv5dv8_cov, iv9dv5_cov, dv5dv9_cov, iv10dv5_cov, dv5dv10_cov,
                       iv1iv6_cov, iv6dv1_cov, iv2iv6_cov, iv6dv2_cov, iv3iv6_cov, iv6dv3_cov, iv4iv6_cov, iv6dv4_cov, iv5iv6_cov, iv6dv5_cov, 1, es6, iv6iv7_cov, iv6dv7_cov, iv6iv8_cov, iv6dv8_cov, iv6iv9_cov, iv6dv9_cov, iv6iv10_cov, iv6dv10_cov,
                       iv1dv6_cov, dv1dv6_cov, iv2dv6_cov, dv2dv6_cov, iv3dv6_cov, dv3dv6_cov, iv4dv6_cov, dv4dv6_cov, iv5dv6_cov, dv5dv6_cov, es6, 1, iv7dv6_cov, dv6dv7_cov, iv8dv6_cov, dv6dv8_cov, iv9dv6_cov, dv6dv9_cov, iv10dv6_cov, dv6dv10_cov,
                       iv1iv7_cov, iv7dv1_cov, iv2iv7_cov, iv7dv2_cov, iv3iv7_cov, iv7dv3_cov, iv4iv7_cov, iv7dv4_cov, iv5iv7_cov, iv7dv5_cov, iv6iv7_cov, iv7dv6_cov, 1, es7, iv7iv8_cov, iv7dv8_cov, iv7iv9_cov, iv7dv9_cov, iv7iv10_cov, iv7dv10_cov,
                       iv1dv7_cov, dv1dv7_cov, iv2dv7_cov, dv2dv7_cov, iv3dv7_cov, dv3dv7_cov, iv4dv7_cov, dv4dv7_cov, iv5dv7_cov, dv5dv7_cov, iv6dv7_cov, dv6dv7_cov, es7, 1, iv8dv7_cov, dv7dv8_cov, iv9dv7_cov, dv7dv9_cov, iv10dv7_cov, dv7dv10_cov,
                       iv1iv8_cov, iv8dv1_cov, iv2iv8_cov, iv8dv2_cov, iv3iv8_cov, iv8dv3_cov, iv4iv8_cov, iv8dv4_cov, iv5iv8_cov, iv8dv5_cov, iv6iv8_cov, iv8dv6_cov, iv7iv8_cov, iv8dv7_cov, 1, es8, iv8iv9_cov, iv8dv9_cov, iv8iv10_cov, iv8dv10_cov,
                       iv1dv8_cov, dv1dv8_cov, iv2dv8_cov, dv2dv8_cov, iv3dv8_cov, dv3dv8_cov, iv4dv8_cov, dv4dv8_cov, iv5dv8_cov, dv5dv8_cov, iv6dv8_cov, dv6dv8_cov, iv7dv8_cov, dv7dv8_cov, es8, 1, iv9dv8_cov, dv8dv9_cov, iv10dv8_cov, dv8dv10_cov,
                       iv1iv9_cov, iv9dv1_cov, iv2iv9_cov, iv9dv2_cov, iv3iv9_cov, iv9dv3_cov, iv4iv9_cov, iv9dv4_cov, iv5iv9_cov, iv9dv5_cov, iv6iv9_cov, iv9dv6_cov, iv7iv9_cov, iv9dv7_cov, iv8iv9_cov, iv9dv8_cov, 1, es9, iv9iv10_cov, iv9dv10_cov,
                       iv1dv9_cov, dv1dv9_cov, iv2dv9_cov, dv2dv9_cov, iv3dv9_cov, dv3dv9_cov, iv4dv9_cov, dv4dv9_cov, iv5dv9_cov, dv5dv9_cov, iv6dv9_cov, dv6dv9_cov, iv7dv9_cov, dv7dv9_cov, iv8dv9_cov, dv8dv9_cov, es9, 1, iv10dv9_cov, dv9dv10_cov,
                       iv1iv10_cov, iv10dv1_cov, iv2iv10_cov, iv10dv2_cov, iv3iv10_cov, iv10dv3_cov, iv4iv10_cov, id10dv4_cov, iv5iv10_cov, iv10dv5_cov, iv6iv10_cov, iv10dv6_cov, iv7iv10_cov, iv10dv7_cov, iv8iv10_cov, iv10dv8_cov, iv9iv10_cov, iv10dv9_cov, 1, es10,
                       iv1dv10_cov, dv1dv10_cov, iv2dv10_cov, dv2dv10_cov, iv3dv10_cov, dv3dv10_cov, iv4dv10_cov, dv4dv10_cov, iv5dv10_cov, dv5dv10_cov, iv6dv10_cov, dv6dv10_cov, iv7dv10_cov, dv7dv10_cov, iv8dv10_cov, dv8dv10_cov, iv9dv10_cov, dv9dv10_cov, es10, 1), 20, 20)

  # Convert effect size into correlation coefficients #
  if (es_units == "d") {
    cortable[1,2] <- cortable[1,2] / sqrt(cortable[1,2]**2 + 4)
    cortable[2,1] <- cortable[2,1] / sqrt(cortable[2,1]**2 + 4)

    cortable[3,4] <- cortable[3,4] / sqrt(cortable[3,4]**2 + 4)
    cortable[4,3] <- cortable[4,3] / sqrt(cortable[4,3]**2 + 4)

    cortable[5,6] <- cortable[5,6] / sqrt(cortable[5,6]**2 + 4)
    cortable[6,5] <- cortable[6,5] / sqrt(cortable[6,5]**2 + 4)

    cortable[7,8] <- cortable[7,8] / sqrt(cortable[7,8]**2 + 4)
    cortable[8,7] <- cortable[8,7] / sqrt(cortable[8,7]**2 + 4)

    cortable[9,10] <- cortable[9,10] / sqrt(cortable[9,10]**2 + 4)
    cortable[10,9] <- cortable[10,9] / sqrt(cortable[10,9]**2 + 4)

    cortable[11,12] <- cortable[11,12] / sqrt(cortable[11,12]**2 + 4)
    cortable[12,11] <- cortable[12,11] / sqrt(cortable[12,11]**2 + 4)

    cortable[13,14] <- cortable[13,14] / sqrt(cortable[13,14]**2 + 4)
    cortable[14,13] <- cortable[14,13] / sqrt(cortable[14,13]**2 + 4)

    cortable[15,16] <- cortable[15,16] / sqrt(cortable[15,16]**2 + 4)
    cortable[16,15] <- cortable[16,15] / sqrt(cortable[16,15]**2 + 4)

    cortable[17,18] <- cortable[17,18] / sqrt(cortable[17,18]**2 + 4)
    cortable[18,17] <- cortable[18,17] / sqrt(cortable[18,17]**2 + 4)

    cortable[19,20] <- cortable[19,20] / sqrt(cortable[19,20]**2 + 4)
    cortable[20,19] <- cortable[20,19] / sqrt(cortable[20,19]**2 + 4)

  } else if (es_units == "r2") {
    cortable[1,2] <- sqrt(cortable[1,2])
    cortable[2,1] <- sqrt(cortable[2,1])

    cortable[3,4] <- sqrt(cortable[3,4])
    cortable[4,3] <- sqrt(cortable[4,3])

    cortable[5,6] <- sqrt(cortable[5,6])
    cortable[6,5] <- sqrt(cortable[6,5])

    cortable[7,8] <- sqrt(cortable[7,8])
    cortable[8,7] <- sqrt(cortable[8,7])

    cortable[9,10] <- sqrt(cortable[9,10])
    cortable[10,9] <- sqrt(cortable[10,9])

    cortable[11,12] <- sqrt(cortable[11,12])
    cortable[12,11] <- sqrt(cortable[12,11])

    cortable[13,14] <- sqrt(cortable[13,14])
    cortable[14,13] <- sqrt(cortable[14,13])

    cortable[15,16] <- sqrt(cortable[15,16])
    cortable[16,15] <- sqrt(cortable[16,15])

    cortable[17,18] <- sqrt(cortable[17,18])
    cortable[18,17] <- sqrt(cortable[18,17])

    cortable[19,20] <- sqrt(cortable[19,20])
    cortable[20,19] <- sqrt(cortable[20,19])
  }

  # Let the user know it's working #
  message(paste("Running", iterations, "sets of tests. This may take a minute.", sep = " "))

  # Set up for parallelization #
  doParallel::registerDoParallel(parallel::detectCores())
  `%dopar%` <- foreach::`%dopar%`

  # Apply Bonferroni correction #
  if (bonferroni == TRUE) {
    alpha <- alpha / models
  }

  # Run the models i times #
  result <- vector()

  result <- foreach::foreach (i=1:iterations, .combine=rbind) %dopar% {
    set.seed(seed + i)

    # Simulate the sample #
    sample <- MASS::mvrnorm(n = n,
                            mu = rep(0, times = 20),
                            Sigma = cortable)

    sample <- data.frame(sample)
    names(sample) <- c("iv1", "dv1",
                     "iv2", "dv2",
                     "iv3", "dv3",
                     "iv4", "dv4",
                     "iv5", "dv5",
                     "iv6", "dv6",
                     "iv7", "dv7",
                     "iv8", "dv8",
                     "iv9", "dv9",
                     "iv10", "dv10")

    if (models == 2) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result)
      result <- round_results
    }
    else if (models == 3) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result)
      result <- round_results
    }
    else if (models == 4) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result)
      result <- round_results
    }
    else if (models == 5) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      es5_result <- model5$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result, es5_result)
      result <- round_results
    }
    else if (models == 6) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      es5_result <- model5$Pr...t..[2] < alpha
      es6_result <- model6$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result, es5_result, es6_result)
      result <- round_results
    }
    else if (models == 7) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      es5_result <- model5$Pr...t..[2] < alpha
      es6_result <- model6$Pr...t..[2] < alpha
      es7_result <- model7$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result, es5_result, es6_result,
                                  es7_result)
      result <- round_results
    }
    else if (models == 8) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      model8 <- data.frame(summary(lm(dv8 ~ iv8, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      es5_result <- model5$Pr...t..[2] < alpha
      es6_result <- model6$Pr...t..[2] < alpha
      es7_result <- model7$Pr...t..[2] < alpha
      es8_result <- model8$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result, es5_result, es6_result,
                                  es7_result, es8_result)
      result <- round_results
    }
    else if (models == 9) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      model8 <- data.frame(summary(lm(dv8 ~ iv8, data = sample))$coefficients)
      model9 <- data.frame(summary(lm(dv9 ~ iv9, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      es5_result <- model5$Pr...t..[2] < alpha
      es6_result <- model6$Pr...t..[2] < alpha
      es7_result <- model7$Pr...t..[2] < alpha
      es8_result <- model8$Pr...t..[2] < alpha
      es9_result <- model9$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result, es5_result, es6_result,
                                  es7_result, es8_result, es9_result)
      result <- round_results
    }
    else if (models == 10) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      model8 <- data.frame(summary(lm(dv8 ~ iv8, data = sample))$coefficients)
      model9 <- data.frame(summary(lm(dv9 ~ iv9, data = sample))$coefficients)
      model10 <- data.frame(summary(lm(dv10 ~ iv10, data = sample))$coefficients)
      es1_result <- model1$Pr...t..[2] < alpha
      es2_result <- model2$Pr...t..[2] < alpha
      es3_result <- model3$Pr...t..[2] < alpha
      es4_result <- model4$Pr...t..[2] < alpha
      es5_result <- model5$Pr...t..[2] < alpha
      es6_result <- model6$Pr...t..[2] < alpha
      es7_result <- model7$Pr...t..[2] < alpha
      es8_result <- model8$Pr...t..[2] < alpha
      es9_result <- model9$Pr...t..[2] < alpha
      es10_result <- model10$Pr...t..[2] < alpha
      round_results <- data.frame(es1_result, es2_result,
                                  es3_result, es4_result, es5_result, es6_result,
                                  es7_result, es8_result, es9_result, es10_result)
      result <- round_results
    }
  }

  # Adjust for "null power" #
  if (1 %in% null_effect) {
    result$es1_result <- ifelse(result$es1_result == TRUE, FALSE,
                                ifelse(result$es1_result == FALSE, TRUE, NA))
  }

  if (2 %in% null_effect) {
    result$es2_result <- ifelse(result$es2_result == TRUE, FALSE,
                                ifelse(result$es2_result == FALSE, TRUE, NA))
  }

  if (3 %in% null_effect) {
    result$es3_result <- ifelse(result$es3_result == TRUE, FALSE,
                                ifelse(result$es3_result == FALSE, TRUE, NA))
  }

  if (4 %in% null_effect) {
    result$es4_result <- ifelse(result$es4_result == TRUE, FALSE,
                                ifelse(result$es4_result == FALSE, TRUE, NA))
  }

  if (5 %in% null_effect) {
    result$es5_result <- ifelse(result$es5_result == TRUE, FALSE,
                                ifelse(result$es5_result == FALSE, TRUE, NA))
  }

  if (6 %in% null_effect) {
    result$es6_result <- ifelse(result$es6_result == TRUE, FALSE,
                                ifelse(result$es6_result == FALSE, TRUE, NA))
  }

  if (7 %in% null_effect) {
    result$es7_result <- ifelse(result$es7_result == TRUE, FALSE,
                                ifelse(result$es7_result == FALSE, TRUE, NA))
  }

  if (8 %in% null_effect) {
    result$es8_result <- ifelse(result$es8_result == TRUE, FALSE,
                                ifelse(result$es8_result == FALSE, TRUE, NA))
  }

  if (9 %in% null_effect) {
    result$es9_result <- ifelse(result$es9_result == TRUE, FALSE,
                                ifelse(result$es9_result == FALSE, TRUE, NA))
  }

  if (10 %in% null_effect) {
    result$es10_result <- ifelse(result$es10_result == TRUE, FALSE,
                                 ifelse(result$es10_result == FALSE, TRUE, NA))
  }

  # Create variable for calculating simultaneous power #
  result$simultaneous <- FALSE

  if (models == 2) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE] <- TRUE
  } else if (models == 3) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE] <- TRUE
  } else if (models == 4) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE] <- TRUE
  } else if (models == 5) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE] <- TRUE
  } else if (models == 6) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                          result$es6_result == TRUE] <- TRUE
  } else if (models == 7) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                          result$es6_result == TRUE & result$es7_result == TRUE] <- TRUE
  } else if (models == 8) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                          result$es6_result == TRUE & result$es7_result == TRUE & result$es8_result == TRUE] <- TRUE
  } else if (models == 9) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                          result$es6_result == TRUE & result$es7_result == TRUE & result$es8_result == TRUE & result$es9_result == TRUE] <- TRUE
  } else if (models == 10) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                          result$es6_result == TRUE & result$es7_result == TRUE & result$es8_result == TRUE & result$es9_result == TRUE & result$es10_result == TRUE] <- TRUE
  }

  # Convert T/F results to numeric for ease of calculation #
  result[result == "TRUE"] <- 1

  # Calculate power #
  simultaneous_power <- mean(result$simultaneous)*100
  model1_power <- mean(result$es1_result)*100
  model2_power <- mean(result$es2_result)*100
  if (models > 2) {
    model3_power <- mean(result$es3_result)*100
  }
  if (models > 3) {
    model4_power <- mean(result$es4_result)*100
  }
  if (models > 4) {
    model5_power <- mean(result$es5_result)*100
  }
  if (models > 5) {
    model6_power <- mean(result$es6_result)*100
  }
  if (models > 6) {
    model7_power <- mean(result$es7_result)*100
  }
  if (models > 7) {
    model8_power <- mean(result$es8_result)*100
  }
  if (models > 8) {
    model9_power <- mean(result$es9_result)*100
  }
  if (models > 9) {
    model10_power <- mean(result$es10_result)*100
  }

  # Summarize results

  if (models == 2) {
    power <- data.frame(model1_power, model2_power, simultaneous_power)
  } else if ( models == 3) {
    power <- data.frame(model1_power, model2_power, model3_power, simultaneous_power)
  } else if ( models == 4) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, simultaneous_power)
  } else if ( models == 5) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, simultaneous_power)
  } else if ( models == 6) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, simultaneous_power)
  } else if ( models == 7) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, simultaneous_power)
  } else if ( models == 8) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, model8_power, simultaneous_power)
  } else if ( models == 9) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, model8_power, model9_power, simultaneous_power)
  } else if ( models == 10) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, model8_power, model9_power, model10_power, simultaneous_power)
  }

  if(print_result == TRUE) {
    for (p in 1:models) {
      cat(paste0("Model ", p, " Power: ", power[p], "%"), sep = "\n")
      }
    cat(paste0("Simultaneous Power: ", power[length(power)], "%"))
    }

  invisible(power)

}
