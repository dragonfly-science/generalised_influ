default_palette <- c(
  'current'  = 'dodgerblue1', 
  'previous' = '#E41A1CCC', 
  'reduced'  = 'grey85',
  'main'     = 'black', 
  'grid'     = 'grey90', 
  'band'     = 'grey80',
  'dharm'    = 'grey30', 
  'helper'   = 'black',
  'extra1'   = '#F5B915FF',
  'extra2'   = '#08235FFF',
  'extra3'   = '#4D9221',
  'extra4'   = "purple4" ,  
  'extra5'   = "violetred",
  'ref'      = 'blue',
  'above'    = '#E41A1CCC',
  'gradient1' = "darkred", 
  'gradient2' = "tomato", 
  'gradient3' = "grey90",
  'gradient4' = "cornflowerblue", 
  'gradient5'  = "darkblue"
)

my_darken <- function(color, factor = 0.7) {
  # Convert Hex to RGB matrix
  rgb_vals <- col2rgb(color) 
  
  # Do the math and convert back to Hex
  rgb(
    rgb_vals["red", ] * factor, 
    rgb_vals["green", ] * factor, 
    rgb_vals["blue", ] * factor, 
    maxColorValue = 255
  )
}