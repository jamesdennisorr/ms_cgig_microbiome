# Define the shared colour, shape and linetype scales used by every figure
#  sourced by 02 through 07, not run on its own
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# User set variables
# Treatment strings must match the Treatment column of the phyloseq sample data
TREATMENT_ORDER <- c("Control", "Antibiotics", "High temperature", "Antibiotics + HT")

# Cap on classes drawn individually in Figures 3 and 4. Wong (2011) supplies
# seven hues plus a grey, so asking for more forces colours the palette cannot
# give and the figure stops being colourblind safe.
N_CLASSES_SHOWN <- 7


#### 01. Treatment scales ####
# Wong (2011) colourblind-safe palette
TREATMENT_COLORS <- c("Control"          = "#999999"
                      , "Antibiotics"      = "#E69F00"
                      , "High temperature" = "#56B4E9"
                      , "Antibiotics + HT" = "#009E73"
                      )

TREATMENT_SHAPES <- c("Control"          = 16   # filled circle
                      , "Antibiotics"      = 17   # filled triangle
                      , "High temperature" = 15   # filled square
                      , "Antibiotics + HT" = 18   # filled diamond
                      )

TREATMENT_LINETYPES <- c("Control"          = "solid"
                         , "Antibiotics"      = "dashed"
                         , "High temperature" = "dotted"
                         , "Antibiotics + HT" = "twodash"
                         )


#### 02. Class scales ####
# Seven Wong hues for the classes shown, grey reserved for the pooled remainder.
# Order is fixed so the same class takes the same colour in Figures 3 and 4.
CLASS_PALETTE <- c("#E69F00"   # orange
                   , "#56B4E9"   # sky blue
                   , "#009E73"   # bluish green
                   , "#F0E442"   # yellow
                   , "#0072B2"   # blue
                   , "#D55E00"   # vermillion
                   , "#CC79A7"   # reddish purple
                   )

CLASS_OTHER_COLOR <- "#999999"

stopifnot(length(CLASS_PALETTE) >= N_CLASSES_SHOWN)

# End
