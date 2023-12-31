% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{PredictiveParameters}
\alias{PredictiveParameters}
\title{Calculation of the predictive parameters}
\usage{
PredictiveParameters(prediction_dt)
}
\arguments{
\item{prediction_dt}{Data frame with the prediction made, the actual 
classification, and the score used for discrimination (the difference 
between the distance of the coordinates of each individual with each 
of the maxima of the two classes). It must be the output of 
the function "ClassAssignment".}
}
\value{
List of calculated predictive parameters (relative confussion
matrix, negative and positive predictive values, F1 score, Cohen's
kappa and Matthew's correlation coefficients)
}
\description{
It calculates from the data frame constructed by the function 
"ClassAssignment "the values of the confusion matrix in relative 
frequencies and other predictive parameters (predictive values, accuracy, 
F1 score, Cohen's Kappa and Matthews Correlation coefficients).
}
