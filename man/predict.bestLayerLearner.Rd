% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.bestLayerLearner.R
\name{predict.bestLayerLearner}
\alias{predict.bestLayerLearner}
\title{Best specific Learner prediction.}
\usage{
\method{predict}{bestLayerLearner}(object, data, ...)
}
\arguments{
\item{object}{\code{bestLayerLearner} \cr
An object from class \link{bestLayerLearner}}

\item{data}{\code{data.frame} \cr
New data to predicted.}

\item{...}{\code{any}\cr
Further arguments passed to or from other methods.}
}
\value{
Predicted target values are returned.
}
\description{
Predict function for models from class \code{bestLayerLearner}.
}
\examples{
set.seed(20240625)
x = data.frame(x1 = runif(n = 50L, min = 0, max = 1))
y <- sample(x = 0:1, size = 50L, replace = TRUE)
my_model <- bestLayerLearner(x = x, y = y)
x_new <- data.frame(x1 = rnorm(10L))
my_predictions <- predict(object = my_model, data = x_new)


}
