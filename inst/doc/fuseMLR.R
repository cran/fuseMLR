## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----libraries, warning = FALSE-----------------------------------------------
library(fuseMLR)

## ----data_exam, include=TRUE, eval=TRUE---------------------------------------
data("multi_omics")
# This list contains two lists of data: training and test.
# Each sublist contains three omics data as described above.
str(object = multi_omics, max.level = 2L)

## ----training, include=TRUE, eval=TRUE----------------------------------------
training <- createTraining(id = "training",
                           ind_col = "IDS",
                           target = "disease",
                           target_df = multi_omics$training$target,
                           verbose = TRUE)
print(training)

## ----geneexpr, include=TRUE, eval=TRUE----------------------------------------
# Create gene expression layer.
createTrainLayer(training = training,
                 train_layer_id = "geneexpr",
                 train_data = multi_omics$training$geneexpr,
                 varsel_package = "Boruta",
                 varsel_fct = "Boruta",
                 varsel_param = list(num.trees = 1000L,
                                     mtry = 3L,
                                     probability = TRUE),
                 lrner_package = "ranger",
                 lrn_fct = "ranger",
                 param_train_list = list(probability = TRUE,
                                         mtry = 1L),
                 param_pred_list = list(),
                 na_action = "na.keep")

## ----proteinexpr, include=TRUE, eval=TRUE-------------------------------------
# Create gene protein abundance layer
createTrainLayer(training = training,
                 train_layer_id = "proteinexpr",
                 train_data = multi_omics$training$proteinexpr,
                 varsel_package = "Boruta",
                 varsel_fct = "Boruta",
                 varsel_param = list(num.trees = 1000L,
                                     mtry = 3L,
                                     probability = TRUE),
                 lrner_package = "ranger",
                 lrn_fct = "ranger",
                 param_train_list = list(probability = TRUE,
                                         mtry = 1L),
                 param_pred_list = list(type = "response"),
                 na_action = "na.keep")

## ----methylation, include=TRUE, eval=TRUE-------------------------------------
# Create methylation layer
createTrainLayer(training = training,
                 train_layer_id = "methylation",
                 train_data = multi_omics$training$methylation,
                 varsel_package = "Boruta",
                 varsel_fct = "Boruta",
                 varsel_param = list(num.trees = 1000L,
                                     mtry = 3L,
                                     probability = TRUE),
                 lrner_package = "ranger",
                 lrn_fct = "ranger",
                 param_train_list = list(probability = TRUE,
                                         mtry = 1L),
                 param_pred_list = list(),
                 na_action = "na.keep")

## ----training_meta_layers, include=TRUE, eval=TRUE, class.output="scroll-100"----
# Create meta layer with imputation of missing values.
createTrainMetaLayer(training = training,
                     meta_layer_id = "meta_layer",
                     lrner_package = NULL,
                     lrn_fct = "weightedMeanLearner",
                     param_train_list = list(),
                     param_pred_list = list(na_rm = TRUE),
                     na_action = "na.rm")
print(training)

## ----upsetplot, include=TRUE, eval=TRUE, fig.width=5, fig.align='center'------
upsetplot(object = training, order.by = "freq")

## ----varsel, include=TRUE, eval=TRUE, warning=FALSE---------------------------
# Variable selection
set.seed(5467)
var_sel_res <- varSelection(training = training)
print(var_sel_res)

## ----training_varsel_disp, include=TRUE, eval=TRUE, warning=FALSE-------------
print(training)

## ----lrner_train, include=TRUE, eval=TRUE, message=TRUE-----------------------
set.seed(5462)
fusemlr(training = training,
        use_var_sel = TRUE)

## ----display_lrner_trained, include=TRUE, eval=TRUE, message=TRUE-------------
print(training)

## ----training_summary, include=TRUE, eval=TRUE, message=TRUE------------------
summary(training)

## ----basic_lrnr, include=TRUE, eval=TRUE--------------------------------------
models_list <- extractModel(training = training)
str(object = models_list, max.level = 1L)

## ----basic_data, include=TRUE, eval=TRUE--------------------------------------
data_list <- extractData(object = training)
str(object = data_list, max.level = 1)

## ----testing, include=TRUE, eval=TRUE-----------------------------------------
# Create testing for predictions
testing <- createTesting(id = "testing",
                         ind_col = "IDS")

## ----testing_ge, include=TRUE, eval=TRUE--------------------------------------
# Create gene expression layer
createTestLayer(testing = testing,
                test_layer_id = "geneexpr",
                test_data = multi_omics$testing$geneexpr)

## ----testing_pr, include=TRUE, eval=TRUE--------------------------------------
# Create gene protein abundance layer
createTestLayer(testing = testing,
                test_layer_id = "proteinexpr",
                test_data = multi_omics$testing$proteinexpr)

## ----testing_me, include=TRUE, eval=TRUE--------------------------------------
# Create methylation layer
createTestLayer(testing = testing,
                test_layer_id = "methylation",
                test_data = multi_omics$testing$methylation)

## ----testing_summary, include=TRUE, eval=TRUE, message=TRUE-------------------
summary(testing)

## ----basic_test_data, include=TRUE, eval=TRUE---------------------------------
data_list <- extractData(object = testing)
str(object = data_list, max.level = 1)

## ----upsetplot_new, include=TRUE, eval=TRUE, fig.width=5, fig.align='center'----
upsetplot(object = testing, order.by = "freq")

## ----new_pred, include=TRUE, eval=TRUE----------------------------------------
predictions <- predict(object = training, testing = testing)
print(predictions)

## ----brier, include=TRUE, eval=TRUE, message = FALSE--------------------------
pred_values <- predictions$predicted_values
actual_pred <- merge(x = pred_values,
                     y = multi_omics$testing$target,
                     by = "IDS",
                     all.y = TRUE)
y <- as.numeric(actual_pred$disease == "1")

# On all patients
perf_bs <- sapply(X = actual_pred[ , 2L:5L], FUN = function (my_pred) {
  bs <- mean((y[complete.cases(my_pred)] - my_pred[complete.cases(my_pred)])^2)
  roc_obj <- pROC::roc(y[complete.cases(my_pred)], my_pred[complete.cases(my_pred)])
  auc <- pROC::auc(roc_obj)
  performances = rbind(bs, auc)
  return(performances)
})
rownames(perf_bs) <- c("BS", "AUC")
print(perf_bs)

## ----interface, include=TRUE, eval=FALSE--------------------------------------
#  # Re-create the gene expression layer with support vector machine as learner.
#  createTrainLayer(training = training,
#                   train_layer_id = "geneexpr",
#                   train_data = multi_omics$training$geneexpr,
#                   varsel_package = "Boruta",
#                   varsel_fct = "Boruta",
#                   varsel_param = list(num.trees = 1000L,
#                                       mtry = 3L,
#                                       probability = TRUE),
#                   lrner_package = "e1071",
#                   lrn_fct = "svm",
#                   param_train_list = list(type = 'C-classification',
#                                           kernel = 'radial',
#                                           probability = TRUE),
#                   param_pred_list = list(probability = TRUE),
#                   na_action = "na.rm",
#                   x_lrn = "x",
#                   y_lrn = "y",
#                   object = "object",
#                   data = "newdata", # Name discrepancy resolved.
#                   extract_pred_fct = function (pred) {
#                     pred <- attr(pred, "probabilities")
#                     return(pred[ , 1L])
#                   }
#  )
#  # Variable selection
#  set.seed(5467)
#  var_sel_res <- varSelection(training = training)
#  set.seed(5462)
#  training <- fusemlr(training = training,
#                      use_var_sel = TRUE)
#  
#  print(training)

## ----wrap_lasso, include=TRUE, eval=FALSE-------------------------------------
#  # We wrap the original functions
#  mylasso <- function (x, y,
#                       nlambda = 25,
#                       nfolds = 5) {
#    # Perform cross-validation to find the optimal lambda
#    cv_lasso <- glmnet::cv.glmnet(x = as.matrix(x), y = y,
#                          family = "binomial",
#                          type.measure = "deviance",
#                          nfolds = nfolds)
#    best_lambda <- cv_lasso$lambda.min
#    lasso_best <- glmnet::glmnet(x = as.matrix(x), y = y,
#                         family = "binomial",
#                         alpha = 1,
#                         lambda = best_lambda
#    )
#    lasso_model <- list(model = lasso_best)
#    class(lasso_model) <- "mylasso"
#    return(lasso_model)
#  }

## ----wrap_predict, include=TRUE, eval=FALSE-----------------------------------
#  # We extend the generic predict function mylasso.
#  predict.mylasso <- function (object, data) {
#    glmnet_pred <- predict(object = object$model,
#                           newx = as.matrix(data),
#                           type = "response",
#                           s = object$model$lambda)
#    return(as.vector(glmnet_pred))
#  }
#  
#  # Re-create the gene expression layer with support vector machine as learner.
#  createTrainMetaLayer(training = training,
#                       meta_layer_id = "meta_layer",
#                       lrner_package = NULL,
#                       lrn_fct = "mylasso",
#                       param_train_list = list(nlambda = 100L),
#                       na_action = "na.impute")
#  set.seed(5462)
#  training <- fusemlr(training = training,
#                      use_var_sel = TRUE)
#  print(training)

## ----add_wrap, include=TRUE, eval=FALSE---------------------------------------
#  # Re-create the gene expression layer with support vector machine as learner.
#  createTrainMetaLayer(training = training,
#                       meta_layer_id = "meta_layer",
#                       lrner_package = NULL,
#                       lrn_fct = "mylasso",
#                       param_train_list = list(nlambda = 100L),
#                       na_action = "na.impute")
#  set.seed(5462)
#  training <- fusemlr(training = training,
#                      use_var_sel = TRUE)
#  print(training)

## ----implemented_learners, include=TRUE, echo=FALSE---------------------------
# Load knitr package
library(knitr)

# Create a data frame
data <- data.frame(
  Leaner = c("weightedMeanLearner", "bestLayerLearner", "cobra"),
  Description = c("The weighted mean meta learner. It uses modality-specific predictions to estimate the weights of the modality-specific models", "The best layer-specific model is used as meta model.", "cobra implements the COBRA (COmBined Regression Alternative), an aggregation method for combining predictions from multiple individual learners")
)

# Generate the table
kable(data, caption = "")

