##############################################################################
##############################    SETTING UP    ##############################
##############################################################################

# Getting Package Statsecol from github which contains the data
remotes::install_github("https://github.com/chrissuthy/statsecol")
# Loading in the data package
library(statsecol)
# Set up libraries 
library(Distance)
library(tidyverse)
library(remotes)
library(MuMIn)

# First fit the most basic model to explore it.
# Distance model with half-normal key function, no adjustment.
hn <- ds(data = bowhead_LT, 
         # Half-normal
         key = "hn", 
         # No adjustment argument
         adjustment = NULL) 

# Setting plot size
par(mfrow = c(1, 2), mar = c(5.5, 4, 2, 2) + 0.1, oma = c(1, 1, 1, 1))

# Plotting the detection function over the distances distribution
plot(hn, 
     which = 2, 
     pl.col = adjustcolor("seagreen", 0.5), 
     border = NULL,
     ylab = "Detection probability (g(x))", 
     xlab = "Distance", 
     main = "A",
     lwd = 1.4,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21)

# Plotting and running the Cramer-von Mises test
gof_ds(hn, 
       main = "B",
       lwd = 1.4,
       pch = 21,
       col = "darkgoldenrod4",
       bg = alpha("darkgoldenrod2", 0.3))
# Create the labels for figures
subtitle = ("Figure 1: Half-Normal No Adjustment Model Fit. 
            A is the model fit over the distribution \nof distances. 
            B is the QQ-plot of the expected vs observed CDF")
# Print out the labels for figures
mtext(side = 1, adj = 0, cex = 0.7, subtitle, outer = TRUE)


##############################################################################
########################    EXPLORATORY ANALYSIS    ##########################
##############################################################################

# Inspect data
str(bowhead_LT)

# Looking at the distribution of the distances
# Setting up the plot size
par(mfrow = c(1,1), 
    mar = c(5, 4, 4, 2) + 0.1) 
# Plot histogram of distances
hist(bowhead_LT$distance,
     # Naming x-axis 
     xlab = "Distance (km)", 
     # Naming y axis
     ylab = "Number of Bowhead Whales",  
     # Naming plot 
     main = "Distances of Bowhead Whales from Line Transects",
     # Changing colour and intensity of bars
     col = adjustcolor("seagreen", 0.5), 
     # Changing orientation of x-axis tick labels
     las = 1) 

##############################################################################
##############################    MODELLING    ###############################
##############################################################################

# Due to the long shoulder seen in the histogram, first fit a hazard rate,
#   then fit a half-normal

# Fitting hazard-rate detection functions with all adjustment terms 

# Hazard-rate with no adjustment 
hr <- ds(data = bowhead_LT, 
         key = "hr", 
         adjustment = NULL ) 

# Hazard-rate with a cosine adjustment 
hr_cos <- ds(data = bowhead_LT, 
             key = "hr",
             adjustment = "cos")

# Hazard-rate with a hermite polynomial adjustment 
hr_herm <- ds(data = bowhead_LT, 
              key = "hr", 
              adjustment = "herm" )

# Hazard-rate with a simple polynomial adjustment 
hr_poly <- ds(data = bowhead_LT, 
              key = "hr", 
              adjustment = "poly" ) 



# Fitting half-normal detection functions with all adjustment terms 

# Half-normal with no adjustment 
hn <- ds(data = bowhead_LT, 
         key = "hn", 
         adjustment = NULL ) 

# Half-normal with a cosine adjustment
hn_cos <- ds(data = bowhead_LT, 
             key = "hn", 
             adjustment = "cos" ) 

# Half-normal with a hermite polynomial adjustment
hn_herm <- ds(data = bowhead_LT, 
              key = "hn", 
              adjustment = "herm" )

# Half-normal with a simple polynomial adjustment
hn_poly <- ds(data = bowhead_LT, 
              key = "hn", 
              adjustment = "poly" )

##############################################################################
###########################    MODEL SELECTION    ############################
##############################################################################

# Comparing Models graphically by fit of detection probability 
#   over distance distribution

# Half-normal detection function over histogram plot
# Setting up plot size and margins
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1) 

# Plotting all 4 half-normal models
plot(hn, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "No Adjustment") 
plot(hn_cos,
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Cosine")
plot(hn_poly, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Polynomial")
plot(hn_herm, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Hermite Polynomial")
title("Half-Normal Models", 
      line = -1, 
      # Naming title and setting location
      outer = TRUE) 


# Hazard-rate Detection function over histogram plot
par(mfrow = c(2, 2), 
    # Setting up plot size and margins
    mar = c(5, 4, 4, 2) + 0.1) 

# Plotting all 4 hazard-rate models
plot(hr, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "No Adjustment")
plot(hr_cos, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Cosine")
plot(hr_poly, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Polynomial")
plot(hr_herm, 
     pl.col = adjustcolor("seagreen", 0.5),
     lwd = 1.2,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Hermite Polynomial")
title("Hazard-Rate Models", 
      line = -1, 
      # Naming title and setting location
      outer = TRUE) 

# Comparing AIC of all the models 
summarize_ds_models(hn, hn_cos, hn_herm, hn_poly, 
                    hr, hr_cos, hr_herm, hr_poly,
                    # Setting Plain instead of equations
                    output = "plain") 

# As the dataset is small, checking whether AICc is a more appropriate measure
# To do this, Takezawa (2014) has said AICc should be used when the ratio of 
# your parameters to number of data points is less than 1:40

# Print summary of hazard-rate distance model
summary(hr)
# Parameters = 2 : Observations = 58
# Therefore ratio is 1:29
#   This ratio will be even smaller for the models with adjustments 

# Print summary of half-normal distance model
summary(hn)
# This is the only model that meets the assumptions of AIC however as you 
#   Can"t compare across model selection parameters we will use AICc

# Using AICc to select models 
AICc(hn, hn_cos, hn_herm, hn_poly, 
     hr, hr_cos, hr_herm, hr_poly)

# Despite this there is no change in the best model

# The half-normal will be chosen as it has the smallest AIC and AICc
# No adjustment will be chosen as the adjustments don"t model extra variability 
#   in the data 

##############################################################################
###############################    MODEL FIT    ##############################
##############################################################################

# Comparing the detection function to the cramer-von mises test

# Setting plot size
par(mfrow = c(1,2), mar = c(5, 4, 4, 2) + 0.1)

# Plotting the detetion function over the distances distribution
plot(hn, 
     which = 2, 
     pl.col = adjustcolor("seagreen", 0.5), 
     border = NULL,
     ylab = "Detection probability (g(x))", 
     xlab = "Distance", 
     las = 1,
     lwd = 1.4,
     bg = alpha("darkgoldenrod2", 0.3),
     pch = 21,
     main = "Half-normal Model No Adjustments")

# Plotting and running the Cramer-von Mises test and bootstrap Kolmogorov-Smirnov      
#   test for goodness-of-fit

gof_ds(hn, 
       main = "Expected vs Observed CDF", 
       lwd = 1.4,
       pch = 21,
       bg = alpha("darkgoldenrod2", 0.3),
       ks = TRUE,
       col = "darkgoldenrod4")

# The Cramer-von Mises test gives a test statistic of 0.0732325 and a p-value 
#   of 0.731882
# The Kolomogorov-Smirnov test gives a test statistic stat of 0.0725551 and a             
#   p-value of 1

# Therefore the model has a good fit as the p values are much more than 0.05

##############################################################################
############################    MODEL INFERENCE    ###########################
##############################################################################
# Print summary of the model
summary(hn)

# Creating subsets of the column combinations needed for dht() density and 
#   abundance estimate and variances function.
# Region labels and areas. 
region_table <- unique(bowhead_LT[, c("Region.Label", 
                                      "Area")])

# Region and sample labels, effort
sample_table <- unique(bowhead_LT[, c("Region.Label", 
                                      "Sample.Label", 
                                      "Effort")])

# Object, region and sample labels
observation_table <- unique(bowhead_LT[, c("object", 
                                           "Region.Label", 
                                           "Sample.Label")])

# Finding the abundance and CI
abund_bio_hn <- dht(model = hn$ddf,
                    region_table, 
                    sample_table, 
                    observation_table)

# Viewing the abundance and CI and getting the tables for plotting

# Estimated Abundance of individuals
N_ind <- abund_bio_hn$individuals$N 

# Estimated Abundance of clusters
N_clu <- abund_bio_hn$clusters$N 

# Estimated density of individuals 
d_ind <- abund_bio_hn$individuals$D

# Estimated density of clusters
d_clu <- abund_bio_hn$clusters$D

# Average probability of detection for individuals
abund_bio_hn$individuals$average.p 

# Average probability of detection for clusters
abund_bio_hn$clusters$average.p 

# Changing the Region label to a factor
# This allows it to be plotted in a desired order
fct_N_ind <- N_ind %>%
  mutate(Label = as.factor(Label),
         Label = fct_relevel(Label, 
                             level = c("2", "3", "9", 
                                       "11", "12", "15",
                                       "Total"))) 


ggplot(fct_N_ind)+
  geom_point(aes(x = Label, y = Estimate), stat = "identity", colour = "red") +
  geom_linerange(aes(x = Label, ymin = lcl, ymax = ucl)) +
  xlab("Regions") +
  ylab("")
  