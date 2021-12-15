#' Pride Flag Color Schemes
#'
#' This function contains a set of palettes based on the colors of various Pride
#' flags. The function's only argument (\code{palette}) is used to select the flag
#' you want to use as the basis of your color scheme. This function works well
#' wrapped inside the \code{scale_color_manual} and \code{scale_fill_manual}
#' functions in \pkg{ggplot2}.
#'
#' @param palette \[\code{character(1)}\]\cr
#'  Which Pride flag would you like to use as the basis of your color scheme?
#' 
#' @details
#' Available color schemes are:
#' \itemize{
#' \item \code{"pride"} -- the six-color Pride Flag
#' \item \code{"philly_poc_pride"} -- the Philadelphia People Of Color Pride Flag
#' \item \code{"gilbert_baker_pride"} -- Gilbert Baker's 1978 Pride Flag
#' \item \code{"bisexual_pride"} -- the Bisexual Pride Flag
#' \item \code{"pansexual_pride"} -- the Pansexual Pride Flag
#' \item \code{"asexual_pride"} -- the Asexual Pride Flag
#' \item \code{"trans_pride"} -- the Trans Pride Flag
#' \item \code{"genderfluid_pride"} -- the Genderfluid Pride Flag
#' \item \code{"genderqueer_pride"} -- the Genderqueer Pride Flag
#' \item \code{"polysexual_pride"} -- the Polysexual Pride Flag
#' \item \code{"agender_pride"} -- the Agender Pride Flag
#' \item \code{"aromantic_pride"} -- the Aromantic Pride Flag
#' \item \code{"nonbinary_pride"} -- the Nonbinary Pride Flag
#' }
#'
#' @author Joel Le Forestier (Twitter: @@JoelLeForestier; Website: joelleforestier.com; Email: joel.leforestier@@mail.utoronto.ca)
#'
#' @examples
#' # A scatterplot using the Philadelphia People Of Color Pride Flag
#'
#' var1 <- rnorm(n = 80)
#' var2 <- rnorm(n = 80)
#' var3 <- as.factor(rep(1:8, times = 10))
#' data <- data.frame(var1, var2, var3)
#'
#' ggplot(data = data, mapping = aes(x = var1, y = var2)) +
#'   scale_color_manual(values = pride_palette("philly_poc_pride")) +
#'   geom_point(aes(color = var3)) +
#'   theme_minimal()
#'
#' @export
pride_palette <- function(palette) {
  if (!is.character(palette) | length(palette) != 1)
    stop("Palette must be a single string")
  
  if (palette %in% names(PRIDE_PALETTES))
    PRIDE_PALETTES[[palette]]
  else
    stop("No such pride flag! Available pride palettes are described under Details section.")
}

PRIDE_PALETTES <- list(
  # Pride Flag
  pride = c("#e40303", # red
            "#ff8c00", # orange
            "#ffed00", # yellow
            "#008026", # green
            "#004dff", # blue
            "#750787"), # purple
  
  # Philly People of Color Pride Flag
  philly_poc_pride = c("#000000", # black
                       "#785017", # brown
                       "#e40303", # red
                       "#ff8c00", # orange
                       "#ffed00", # yellow
                       "#008026", # green
                       "#004dff", # blue
                       "#750787"), # purple
  
  # Gilbert Baker Pride Flag
  gilbert_baker_pride = c("#f564ae", # pink
                          "#e40303", # red
                          "#ff8c00", # orange
                          "#ffed00", # yellow
                          "#008026", # green
                          "#52ced9", # light blue
                          "#391294", # dark purple
                          "#750787"), # red-purple
  
  # Bisexual Pride Flag
  bisexual_pride = c("#D70270", # pink
                     "#734F96", # purple
                     "#0038A8"), # blue
  
  # Pansexual Pride Flag
  pansexual_pride = c("#ff218c", # pink
                      "#ffd800", # yellow
                      "#21b1ff"), # blue
  
  # Asexual Pride Flag
  asexual_pride = c("#000000", # black
                    "#a3a3a3", # grey
                    "#ffffff", # white
                    "#800080"), # purple
  
  # Trans Pride Flag
  trans_pride = c("#5bcdfa", # blue
                  "#f5a9b8", # pink
                  "#ffffff"), # white
  
  # Gender Fluid Pride Flag
  genderfluid_pride = c("#ff75a3", # pink
                        "#ffffff", # white
                        "#bd18d6", # purple
                        "#000000", # black
                        "#333fbd"), # blue
  
  # Genderqueer Pride Flag
  genderqueer_pride = c("#B77FDD", # purple
                        "#ffffff", # white
                        "#48821E"), # green
  
  # Polysexual Pride Flag
  polysexual_pride = c("#F61CB9", # pink
                       "#07D569", # green
                       "#1C92F6"), # blue
  
  # Agender Pride Flag
  agender_pride = c("#000000", # black
                    "#a3a3a3", # grey
                    "#ffffff", # white
                    "#b8f483"), # green
  
  # Aromantic Pride Flag
  aromantic_pride = c("#3da542", # green
                      "#a7d379", # light green
                      "#ffffff", # white
                      "#a3a3a3", # grey
                      "#000000"), # black
  
  # Nonbinary Pride Flag
  nonbinary_pride = c("#fff530", # yellow
                      "#ffffff", # white
                      "#9d59d1", # purple
                      "#000000") # black
)
