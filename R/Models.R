# MATH 688: Data Analytics Capstone II
# Capstone Project: Model Building and Assessment
# Author: Isaiah Steinke
# Last Modified: July 2, 2021
# Written, Debugged, and Tested in R v. 4.1.0

# ===================================================================
# Load required libraries
# ===================================================================
library(ggplot2) # v. 3.3.5
library(ggpubr)  # v. 0.4.0
library(caret)   # v. 6.0-88
library(leaps)   # v. 3.1
library(glmnet)  # v. 4.1-2
library(rpart)   # v. 4.1-15
library(ranger)  # v. 0.12.1
library(gbm)     # v. 2.1.8

# ===================================================================
# Import dataset
# ===================================================================
heat <- read.csv("Full_Dataset.csv", header = TRUE)

# ===================================================================
# Process dataset
# ===================================================================
# From our exploratory data analyses, we found that six features
# had the same value of zero for all compounds. Thus, these features
# can be removed since they will have no explanatory use in the
# models. In addition, we will remove all one-component compounds
# since the heat of formation is zero for these compounds. Thus, 
# our final model process will first check for the number of 
# components before applying the ML model to two- and three-component
# materials.

# Remove useless features
del.features <- c("mean_NfUnfilled", "maxdiff_NfUnfilled",
                  "dev_NfUnfilled", "max_NfUnfilled",
                  "min_NfUnfilled", "most_NfUnfilled")
heat <- within(heat, rm(list = del.features))

# Alternative method
del.features.indexes <- which(names(heat) %in% del.features)
heat <- heat[, -del.features.indexes]

# Remove one-component materials
del.onecomp <- which(heat$NComp == 1)
heat <- heat[-del.onecomp, ]

# Clean up environment
rm(del.features, del.onecomp)

# Create training and test sets
# -----------------------------
# Since we have a good amount of data, we'll use a 70/30 split
# for the training/test sets.
set.seed(26448)
indexes <- createDataPartition(heat$HeatForm, times = 1,
                               p = 0.7, list = FALSE)
heat.train <- heat[indexes, ]
heat.test <- heat[-indexes, ]

# The first column, "Formula," is simply the name of the compound
# and only useful as a label. To make things easier in the
# following analyses, we'll store the formulas in separate vectors
# for the training/test sets and just have the predictors and 
# response variables in our dataset.

heat.train.labels <- heat.train$Formula
heat.train <- heat.train[, -1]
heat.test.labels <- heat.test$Formula
heat.test <- heat.test[, -1]

# Create a model matrix for the test data that will be useful
# when computing predictions
testdata <- model.matrix(HeatForm ~ ., heat.test)

# Clean up environment
rm(indexes)

# ===================================================================
# Set the cross-validation parameters for caret
# ===================================================================
cv.params <- trainControl(method = "repeatedcv", number = 10,
                          repeats = 5)

# Note: In the following, we build models with the same random
# number seed to ensure that the cross-validation procedure is
# performed in the same way for each method.

# ===================================================================
# Regression models: Subset selection methods
# ===================================================================
# We'll be using regsubsets in caret via the train function.
# Since the output of regsubsets doesn't provide easy access to
# predictions, we'll have to do some more work with the output.

# Forward selection
# -----------------
# Rather fast, so we'll have it fit up to 100 variables.
set.seed(723)
fwd.model <- train(HeatForm ~ ., data = heat.train,
                   method = "leapForward", trControl = cv.params,
                   tuneGrid = data.frame(nvmax = 1:100))

# Analyze results; compute test-set metrics
fwd.model$results
write.csv(fwd.model$results, file = "fwd.csv", row.names = FALSE)
plot(fwd.model, metric = "RMSE")
plot(fwd.model, metric = "Rsquared")
plot(fwd.model, metric = "MAE")

# Select model with 50 predictors since models with more
# predictors suffer from diminishing returns.

coef.values <- coef(fwd.model$finalModel, 50)
coef.names <- names(coef.values)
fwd.preds <- testdata[, coef.names] %*% coef.values
fwd.err <- heat.test$HeatForm - fwd.preds
sqrt(mean(fwd.err^2)) # RMSE
mean(abs(fwd.err)) # MAE

# Backward selection
# ------------------
# Again, this is fairly fast with only one thread.
set.seed(723)
bwd.model <- train(HeatForm ~ ., data = heat.train,
                   method = "leapBackward", trControl = cv.params,
                   tuneGrid = data.frame(nvmax = 1:100))

# Analyze results; compute test-set metrics
bwd.model$results
write.csv(bwd.model$results, file = "bwd.csv", row.names = FALSE)
plot(bwd.model, metric = "RMSE")
plot(bwd.model, metric = "Rsquared")
plot(bwd.model, metric = "MAE")

# Select model with 50 predictors since models with more
# predictors suffer from diminishing returns.

coef.values <- coef(bwd.model$finalModel, 50)
coef.names <- names(coef.values)
bwd.preds <- testdata[, coef.names] %*% coef.values
bwd.err <- heat.test$HeatForm - bwd.preds
sqrt(mean(bwd.err^2)) # RMSE
mean(abs(bwd.err)) # MAE

# Stepwise selection
# ------------------
# We have to limit nvmax to 20; higher values require a long
# time to finish, even with multiple threads. (In fact, I had
# to force-close R after an hour with nvmax = 35.)
set.seed(723)
step.model <- train(HeatForm ~ ., data = heat.train,
                    method = "leapSeq", trControl = cv.params,
                    tuneGrid = data.frame(nvmax = 1:20))

# Analyze results; compute test-set metrics
step.model$results
write.csv(step.model$results, file = "step.csv", row.names = FALSE)
plot(step.model, metric = "RMSE")
plot(step.model, metric = "Rsquared")
plot(step.model, metric = "MAE")

# Select model with 20 predictors. It does appear that performance
# can be improved with more predictors, but it takes too long if
# I increase nvmax (even to 25).

coef.values <- coef(step.model$finalModel, 20)
coef.names <- names(coef.values)
step.preds <- testdata[, coef.names] %*% coef.values
step.err <- heat.test$HeatForm - step.preds
sqrt(mean(step.err^2)) # RMSE
mean(abs(step.err)) # MAE

# Clean up environment
rm(coef.values, coef.names)

# ===================================================================
# Regression models: Ridge regression & lasso
# ===================================================================
# Use glmnet; alpha = 0 for ridge regression and alpha = 1 for
# lasso. Tune lambda values.

# Ridge regression
# ----------------
# Set tuning parameters
rr.params <- expand.grid(alpha = 0,
                         lambda = 10^seq(3, -6, length = 10))

# Build model
set.seed(723)
rr.model <- train(HeatForm ~ ., data = heat.train,
                  method = "glmnet", trControl = cv.params,
                  tuneGrid = rr.params)

# Analyze results; compute test-set metrics
rr.model$results
write.csv(rr.model$results, file = "rr.csv", row.names = FALSE)

# Select lambda = 1e-02 as the best parameter value.

rr.preds <- predict(rr.model$finalModel, s = 1e-02,
                    newx = testdata[, -1])
rr.err <- heat.test$HeatForm - rr.preds
sqrt(mean(rr.err^2)) # RMSE
mean(abs(rr.err)) # MAE

# Lasso
# -----
# Set tuning parameters
lasso.params <- expand.grid(alpha = 1,
                            lambda = 10^seq(3, -6, length = 10))

# Build model
set.seed(723)
lasso.model <- train(HeatForm ~ ., data = heat.train,
                     method = "glmnet", trControl = cv.params,
                     tuneGrid = lasso.params)

# Analyze results; compute test-set metrics
lasso.model$results
write.csv(lasso.model$results, file = "lasso.csv", row.names = FALSE)

# Select lambda = 1e-05 as the best parameter value.

lasso.preds <- predict(lasso.model$finalModel, s = 1e-05,
                       newx = testdata[, -1])
lasso.err <- heat.test$HeatForm - lasso.preds
sqrt(mean(lasso.err^2)) # RMSE
mean(abs(lasso.err)) # MAE
lasso.coef <- predict(lasso.model$finalModel, type = "coefficients",
                      s = 1e-05)

# Clean up environment
rm(rr.params, lasso.params)

# ===================================================================
# Tree-based methods
# ===================================================================
# For the basic decision tree, we'll use rpart. For random forests, 
# we'll use ranger since we can leverage the built-in parallel
# processing functionality. Finally, we'll use gbm for boosted
# decision trees.

# Decision tree
# -------------
set.seed(723)
tree.model <- train(HeatForm ~ ., data = heat.train,
                    method = "rpart", trControl = cv.params)

# Analyze results; compute test-set metrics
tree.model$results
tree.preds <- predict(tree.model$finalModel, newdata = heat.test)
tree.err <- heat.test$HeatForm - tree.preds
sqrt(mean(tree.err^2)) # RMSE
mean(abs(tree.err)) # MAE

# Variable importance
tree.vi <- tree.model$finalModel[[12]]
barchart(tree.vi[12:1], xlab = "Variable Importance")

# Random forests
# --------------
# Since caret doesn't allow the number of trees to be tuned,
# we'll set it by using the rule of thumb of 10 × the number of
# predictors (~1500). We'll tune mtry and min.node.size and 
# leave splitrule set to "variance." For mtry, the default for
# regression trees is 1/3 × the number of predictors (~47); so,
# we will tune around this value.

# Set tuning parameters
rf.params <- expand.grid(mtry = c(10, 20, 30, 40, 50, 60, 70,
                                  80, 90, 100),
                         splitrule = "variance",
                         min.node.size = c(5, 10))

# Build model
set.seed(723)
rf.model <- train(HeatForm ~ ., data = heat.train,
                  method = "ranger", trControl = cv.params,
                  tuneGrid = rf.params, num.trees = 1500,
                  importance = "impurity", num.threads = 12)

# Analyze results; compute test-set metrics
rf.model$results
write.csv(rf.model$results, file = "rf.csv", row.names = FALSE)

# Best hyperparameter values: mtry = 30, min.node.size = 5

rf.preds <- predict(rf.model$finalModel, data = heat.test)
rf.err <- heat.test$HeatForm - rf.preds$predictions
sqrt(mean(rf.err^2)) # RMSE
mean(abs(rf.err)) # MAE

# Variable importance
rf.vi <- rf.model$finalModel[[6]] 
rf.vi <- sort(rf.vi, decreasing = TRUE) 
barchart(rf.vi[25:1], xlab = "Variable Importance") # plot top 25

# Boosted decision trees
# ----------------------
# For boosted decision trees, we'll use gbm in caret. We'll tune
# the number of trees, the interaction depth, and the shrinkage and
# leave n.minobsinnode to the default value.

# Set tuning parameters
boost.params <- expand.grid(n.trees = 1500,
                            interaction.depth = c(3, 5, 7),
                            shrinkage = c(0.1, 0.5, 1),
                            n.minobsinnode = 10)

# Build model
set.seed(723)
boost.model <- train(HeatForm ~ ., data = heat.train,
                     method = "gbm", trControl = cv.params,
                     tuneGrid = boost.params,
                     distribution = "gaussian")

# Analyze results; compute test-set metrics
boost.model$results
write.csv(boost.model$results, file = "boost.csv",
          row.names = FALSE)

# Best hyperparameter values: interaction.depth = 7, shrinkage = 0.1

boost.preds <- predict(boost.model$finalModel, newdata = heat.test,
                       n.trees = 1500)
boost.err <- heat.test$HeatForm - boost.preds
sqrt(mean(boost.err^2)) # RMSE
mean(abs(boost.err)) # MAE

# Variable importance
summary(boost.model$finalModel, plotit = FALSE)

boost.vi <- summary(boost.model$finalModel, plotit = FALSE)
boost.vi.vec <- boost.vi$rel.inf
names(boost.vi.vec) <- boost.vi$var
barchart(boost.vi.vec[25:1], xlab = "Variable Importance")

# Clean up environment
rm(rf.params, boost.params)

# ===================================================================
# Other figures
# ===================================================================
# Histograms of the test-set errors
h1 <- ggplot(NULL, aes(x = fwd.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(a) Forward",
               size = 4, hjust = 0)
h2 <- ggplot(NULL, aes(x = bwd.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(b) Backward",
               size = 4, hjust = 0)
h3 <- ggplot(NULL, aes(x = step.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(c) Stepwise",
               size = 4, hjust = 0)
h4 <- ggplot(NULL, aes(x = rr.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(d) Ridge Regression",
               size = 4, hjust = 0)
h5 <- ggplot(NULL, aes(x = lasso.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(e) Lasso",
               size = 4, hjust = 0)
h6 <- ggplot(NULL, aes(x = tree.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(f) Decision Tree",
               size = 4, hjust = 0)
h7 <- ggplot(NULL, aes(x = rf.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(g) Random Forests",
               size = 4, hjust = 0)
h8 <- ggplot(NULL, aes(x = boost.err)) +
      geom_histogram(binwidth = 0.1,
                     fill = "white",
                     color = "black") +
      xlim(-1.5, 1.5) + ylim(0, 1500) +
      xlab("Test-Set Error") + ylab("Count") + theme_bw() +
      annotate("text", x = -1.5, y = 1500,
               label = "(h) Boosting",
               size = 4, hjust = 0)
ggarrange(h1, h2, h3, h4, h5, h6, h7, h8,
          ncol = 3, nrow = 3)