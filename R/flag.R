#' Generate a Pride Flag
#'
#' Using the pride_palette function and ggplot2, create a Pride flag. The function's only argument (flag) is used to select the flag you want to generate.
#'
#' @usage flag(flag = NULL)
#'
#' @param palette Which Pride flag would you like to generate? Options are: the six-color Pride Flag ("pride"), the Philadelphia People Of Color Pride Flag ("philly_poc_pride"), Gilbert Baker's 1978 Pride Flag ("gilbert_baker_pride"), the Bisexual Pride Flag ("bisexual_pride"), the Pansexual Pride Flag ("pansexual_pride"), the Asexual Pride Flag ("asexual_pride"), the Trans Pride Flag ("trans_pride"), the Genderfluid Pride Flag ("genderfluid_pride"), the Genderqueer Pride Flag ("genderqueer_pride"), the Polysexual Pride Flag ("polysexual_pride"), the Agender Pride Flag ("agender_pride"), the Aromantic Pride Flag ("aromantic_pride"), and the Nonbinary Pride Flag ("nonbinary_pride").
#'
#' @author Joel Le Forestier (Twitter: @@JoelLeForestier; Website: joelleforestier.com; Email: joel.leforestier@@mail.utoronto.ca)
#'
#' @examples # Generate a Trans Pride Flag
#'
#' flag(flag = "trans_pride)
#'
#' @export

flag <- function(flag) {

  if (flag == "pride") {
    x <- rep(1:2, times = 6)
    y <- rep(1:6, each = 2)
    pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "philly_poc_pride") {
    x <- rep(1:2, times = 8)
    y <- rep(1:8, each = 2)
    philly_poc_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = philly_poc_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "philly_poc_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "gilbert_baker_pride") {
    x <- rep(1:2, times = 8)
    y <- rep(1:8, each = 2)
    gilbert_baker_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = gilbert_baker_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "gilbert_baker_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "bisexual_pride") {
    x <- rep(1:2, times = 5)
    y <- rep(1:5, each = 2)
    bisexual_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = bisexual_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = c(pride_palette(palette = "bisexual_pride")[1],
                                   pride_palette(palette = "bisexual_pride")[1],
                                   pride_palette(palette = "bisexual_pride")[2],
                                   pride_palette(palette = "bisexual_pride")[3],
                                   pride_palette(palette = "bisexual_pride")[3])) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "pansexual_pride") {
    x <- rep(1:2, times = 3)
    y <- rep(1:3, each = 2)
    pansexual_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = pansexual_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "pansexual_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "asexual_pride") {
    x <- rep(1:2, times = 4)
    y <- rep(1:4, each = 2)
    asexual_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = asexual_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "asexual_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "trans_pride") {
    x <- rep(1:2, times = 5)
    y <- rep(1:5, each = 2)
    trans_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = trans_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = c(pride_palette(palette = "trans_pride")[1],
                                   pride_palette(palette = "trans_pride")[2],
                                   pride_palette(palette = "trans_pride")[3],
                                   pride_palette(palette = "trans_pride")[2],
                                   pride_palette(palette = "trans_pride")[1])) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "genderfluid_pride") {
    x <- rep(1:2, times = 5)
    y <- rep(1:5, each = 2)
    genderfluid_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = genderfluid_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "genderfluid_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "genderqueer_pride") {
    x <- rep(1:2, times = 3)
    y <- rep(1:3, each = 2)
    genderqueer_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = genderqueer_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "genderqueer_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "polysexual_pride") {
    x <- rep(1:2, times = 3)
    y <- rep(1:3, each = 2)
    polysexual_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = polysexual_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "polysexual_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "agender_pride") {
    x <- rep(1:2, times = 7)
    y <- rep(1:7, each = 2)
    agender_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = agender_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = c(pride_palette(palette = "agender_pride")[1],
                                   pride_palette(palette = "agender_pride")[2],
                                   pride_palette(palette = "agender_pride")[3],
                                   pride_palette(palette = "agender_pride")[4],
                                   pride_palette(palette = "agender_pride")[3],
                                   pride_palette(palette = "agender_pride")[2],
                                   pride_palette(palette = "agender_pride")[1])) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "aromantic_pride") {
    x <- rep(1:2, times = 5)
    y <- rep(1:5, each = 2)
    aromantic_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = aromantic_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "aromantic_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")

  } else if (flag == "nonbinary_pride") {
    x <- rep(1:2, times = 4)
    y <- rep(1:4, each = 2)
    nonbinary_pride_data <- data.frame(x, y)

    ggplot2::ggplot(data = nonbinary_pride_data, aes(x, ymin = rev(y), ymax = rev(y+1), fill=factor(y))) +
      scale_fill_manual(values = pride_palette(palette = "nonbinary_pride")) +
      geom_ribbon() +
      theme_void() +
      theme(legend.position = "none")
  }
}

