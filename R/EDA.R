# MATH 688: Data Analytics Capstone II
# Capstone Project: Exploratory Data Analyses
# Author: Isaiah Steinke
# Last Modified: June 17, 2021
# Written, Debugged, and Tested in R v. 4.1.0

# Load required libraries
# -----------------------
library(ggplot2) # v. 3.3.4
library(ggpubr) # v. 0.4.0

# Import dataset
# --------------
heat <- read.csv("Full_Dataset.csv", header = TRUE)

# Look at the number of components
# --------------------------------
hist(heat$NComp, main = NULL, xlab = "Number of Components")
sum(heat$NComp == 1) # 64
sum(heat$NComp == 2) # 2743
sum(heat$NComp == 3) # 6783

# Summary statistics & plots: original data excluding Magpie features
# -------------------------------------------------------------------
summary(subset(heat, select = c(Volume, Mass, NoAtoms, HeatForm)))
cor(subset(heat, select = c(Volume, Mass, NoAtoms, HeatForm)))
pairs(subset(heat, select = c(Volume, Mass, NoAtoms, HeatForm)))

# Notably, none of the original predictors have significant
# correlations with HeatForm (the largest is 0.16 in magnitude).

# Violin plots
v1 <- ggplot(heat, aes(factor(NComp), Volume)) + geom_violin() +
      xlab("Number of Components") +
      ylab(expression("Volume of Unit Cell [Å"^{~3}*"]")) +
      theme_bw()
v2 <- ggplot(heat, aes(factor(NComp), Mass)) + geom_violin() +
      xlab("Number of Components") +
      ylab("Atomic Mass [g/mol]") +
      theme_bw()
v3 <- ggplot(heat, aes(factor(NComp), NoAtoms)) + geom_violin() +
      xlab("Number of Components") +
      ylab("Number of Atoms") +
      theme_bw()
v4 <- ggplot(heat, aes(factor(NComp), HeatForm)) + geom_violin() +
      xlab("Number of Components") +
      ylab("Heat of Formation [eV/atom]") +
      theme_bw()
ggarrange(v1, v2, v3, v4, ncol = 2, nrow = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"))

# Scatterplots
ggplot(heat, aes(Volume, Mass, color = factor(NComp))) +
  geom_point(shape = 1) + theme_bw()
ggplot(heat, aes(Volume, NoAtoms, color = factor(NComp))) +
  geom_point(shape = 1) + theme_bw()
ggplot(heat, aes(Volume, HeatForm, color = factor(NComp))) +
  geom_point(shape = 1) + theme_bw()
ggplot(heat, aes(Mass, NoAtoms, color = factor(NComp))) +
  geom_point(shape = 1) + theme_bw()
ggplot(heat, aes(Mass, HeatForm, color = factor(NComp))) +
  geom_point(shape = 1) + theme_bw()
ggplot(heat, aes(NoAtoms, HeatForm, color = factor(NComp))) +
  geom_point(shape = 1) + theme_bw()

# It's difficult to distinguish points in these scatterplots since
# there are a lot of data points. Even using "alpha" argument (instead
# of "shape") doesn't improve things much.

# Set global plotting parameters
par(las = 1, mar = c(5.1, 6.1, 2.1, 2.1))

# Boxplots
par(mfrow = c(2, 2))
boxplot(Volume ~ factor(NComp), data = heat,
        xlab = "Number of Components",
        ylab = expression("Volume of Unit Cell [Å"^{~~3}*"]"))
text(0.6, 2150, "(a)")
boxplot(Mass ~ factor(NComp), data = heat,
        xlab = "Number of Components",
        ylab = "Atomic Mass [g/mol]")
text(0.6, 7250, "(b)")
boxplot(NoAtoms ~ factor(NComp), data = heat,
        xlab = "Number of Components",
        ylab = "Number of Atoms")
text(0.6, 175, "(c)")
boxplot(HeatForm ~ factor(NComp), data = heat,
        xlab = "Number of Components",
        ylab = "Heat of Formation [eV/atom]")
text(0.6, -3.9, "(d)")
par(mfrow = c(1, 1))

# Scatterplots; separate by number of components since the number of
# data points to plot is high
plot(heat$Mass[heat$NComp == 1], heat$Volume[heat$NComp == 1],
     main = NULL)
plot(heat$Mass[heat$NComp == 2], heat$Volume[heat$NComp == 2],
     main = NULL)
plot(heat$Mass[heat$NComp == 3], heat$Volume[heat$NComp == 3],
     main = NULL)

# It's nice to see the how the data separate by NComp. However, this
# would result in a large number of plots (6 × 3 = 18), and there
# doesn't appear to be any deep insights.

# Reset plot margins
par(mar = c(5.1, 4.1, 4.1, 2.1))

# Instead, let's just look at HeatForm vs. Mass, Volume, and NoAtoms
# for two and three components.
par(mfrow = c(2, 3))
plot(heat$Volume[heat$NComp == 2], heat$HeatForm[heat$NComp == 2],
     main = NULL, xlab = expression("Volume of Unit Cell [Å"^{3}*"]"),
     ylab = "Heat of Formation [eV/atom]")
text(2125, -4, "(a)")
plot(heat$Mass[heat$NComp == 2], heat$HeatForm[heat$NComp == 2],
     main = NULL, xlab = expression("Atomic Mass [g/mol]"),
     ylab = "Heat of Formation [eV/atom]")
text(7200, -4, "(b)")
plot(heat$NoAtoms[heat$NComp == 2], heat$HeatForm[heat$NComp == 2],
     main = NULL, xlab = expression("Number of Atoms"),
     ylab = "Heat of Formation [eV/atom]")
text(175, -4, "(c)")
plot(heat$Volume[heat$NComp == 3], heat$HeatForm[heat$NComp == 3],
     main = NULL, xlab = expression("Volume of Unit Cell [Å"^{3}*"]"),
     ylab = "Heat of Formation [eV/atom]")
text(2175, -4, "(d)")
plot(heat$Mass[heat$NComp == 3], heat$HeatForm[heat$NComp == 3],
     main = NULL, xlab = expression("Atomic Mass [g/mol]"),
     ylab = "Heat of Formation [eV/atom]")
text(7250, -4, "(e)")
plot(heat$NoAtoms[heat$NComp == 3], heat$HeatForm[heat$NComp == 3],
     main = NULL, xlab = expression("Number of Atoms"),
     ylab = "Heat of Formation [eV/atom]")
text(168, -4, "(f)")
par(mfrow = c(1, 1))

# Reset plot area
par(las = 0)

# Summary statistics: all features
# ---------------------------------
summary(heat[, 2:150])
cor(heat[, 2:149], heat$HeatForm)

# Most features have correlations with HeatForm that are <0.5.
# Some of the features related to CovalentRadius, Electronegativity,
# NpValence, SpaceGroupNumber, frac_pValence, frac_dValence, MaxIonicChar,
# and MeanIonicChar have correlations with HeatForm that are >0.5. The
# largest is for MeanIonicChar (-0.788977848). In addition, there are
# features that have NA for the correlations; this is because these
# features have the same value of zero for every compound. These features
# are the six features related to NfUnfilled. Thus, these features should
# be removed before building models since they will be useless for
# prediction.