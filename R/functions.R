#' @title Data preparation
#' @description 
#' Transforms the data.frame with predictor and response 
#' variables into a data frame with dummy variables. 
#' @param db Database predictor and response variables. 
#' The two last columns must contain the response variable.
#' The variables must be dummy encoded. You can use the makeNominalData function 
#' of the package ExPosition.
#' @param n_model Number of individuals to build the map. The number 
#' of individuals in the test set is calculated by difference.
#' @return The database with the dummy variables. The antepenultimate 
#' and penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test) 
#' to which the individual has been assigned.
#' @import ExPosition
#' @import dendextend
#' @export
DataPreparation = function(db, n_model) {
  # divide the database into training and validation groups
  n_test <- nrow(db) - n_model
  index_model <- sample(nrow(db), size=n_model)
  index_test <- setdiff(1:nrow(db),index_model)
  group <- as.vector(rep(NA, nrow(db)))
  group[index_model] <- rep("model", length(index_model))
  group[index_test] <- rep("test", length(index_test))
  db_post <- data.frame(db, group)
  
  colnames(db_post)[c(length(db_post)-2, length(db_post)-1)] <- c("class0",
                                                                  "class1")
  
  return(db_post) 
}

#' @title Construction of the base
#' @description 
#' It gets the new base obtained after the singular value decomposition
#' @param db_post The database with the dummy variables. The antepenultimate 
#' and penultimate column must correspond to the response variable. 
#' The last column must correspond to the group (training or test) 
#' to which the individual has been assigned. It must be the output 
#' of the function "DataPreparation".
#' @return The base obtained after the singular value decomposition 
#' of the density matrix. Only vectors whose singular values are 
#' greater than 0 are retained.
#' @export
BaseConstruction= function(db_post) {
  predictors <- as.matrix(db_post[db_post$group == "model", 
                                  c(1:(length(db_post)-3))])
  response <- as.matrix(db_post[db_post$group == "model", 
                                c("class0", "class1")])
  D <- t(response)%*%predictors
  X <- sqrt(D)
  
  F <- t(X)%*%X
  F <- (1/sum(diag(F)))*F  # density matrix
  svd_result <- svd(F)
  
  base <-  svd(F)$u[, which(svd_result$d > .Machine$double.eps)]
}

#' @title Calculation of the new coordinates
#' @description 
#' It calculates the coordinates in the new basis of the vectors 
#' of the categorical co-variables
#' @param db_post The database with the dummy variables. the antepenultimate 
#' and penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test) 
#' to which the individual has been assigned. It must be the output 
#' of the function "DataPreparation".
#' @param base The base obtained after the singular value decomposition 
#' of the density matrix
#' @return The data frame with the coordinates in the new basis. The 
#' penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test).
#' @export
CalculateNewCoordinates <- function(db_post, base) {
  predictors <- as.matrix(db_post[, c(1:(length(db_post)-3))])
  new_coordinates <- t(apply(predictors, 1, 
                             function(y) t(base)%*%(y/sqrt(sum(y)))))
  new_coordinates <- data.frame(new_coordinates, db_post$class1, 
                                db_post$group) 
  colnames(new_coordinates) <- c("coor1", "coor2", "response", "group")
  
  return(new_coordinates)
}

#' @title Calculation of the polar coordinates
#' @description 
#' It calculates the polar coordinates from the data frame containing 
#' the coordinate information in the new base.
#' @param new_coordinates The data frame with the coordinates in the new basis. 
#' The penultimate column must correspond to the response variable. 
#' The last column must correspond to the group (training or test). 
#' It must be the output of the function "CalculateNewCoordinates".
#' @param db_post The database with the dummy variables. The antepenultimate 
#' and penultimate column must correspond to the response variable. 
#' The last column must correspond to the group (training or test) 
#' to which the individual has been assigned. It must be the output 
#' of the function "DataPreparation".
#' @return The data frame with the polar coordinates. The 
#' penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test).
#' @export
CalculatePolarCoordinates <- function(new_coordinates, db_post) {
  theta <- atan(new_coordinates[,"coor2"]/new_coordinates[,"coor1"])
  rho <- sqrt((new_coordinates[,"coor1"]^2) + (new_coordinates[,"coor2"]^2))
  
  polar_coor <- data.frame(theta, rho, db_post$class1,  db_post$group)
  colnames(polar_coor) <- c("theta", "rho", "response", "group")
  
  return(polar_coor)
} 

#' @title Calculation of the maximum of the density functions of each class
#' @description 
#' It calculates the maximum of the density functions of the coordinate(s)
#' of the "model" group (training set) of each class. If two-coordinate maxima 
#' are calculated, the Gaussian kernel is used, if one-coordinate maxima 
#' are calculated, the Epanechnikov kernel is used.
#' @param coor The data frame  with the coordinates for which you want to 
#' calculate the density maxima. The penultimate column must correspond to 
#' the response variable. The last column must correspond to the 
#' group (training or test). It must be the output of the function 
#' "CalculateNewCoordinates" or "CalculatePolarCoordinates".
#' @param coordinates Coordinates for which the maxima of density functions
#' are to be calculated. By default, it is calculated for two coordinates, 
#' corresponding to the first two columns of the "coor" data frame. 
#' If you only want to do it for one, you must indicate the name of 
#' the column containing that coordinate.
#' @return Matrix with the maxima of the density functions of each class (rows) 
#' @import ks
#' @import stats
#' @export
MaximumDensityFunctions <- function(coor, coordinates = "two") {
  if (coordinates == "two") {
    den_class0 <- kde(coor[coor$group == "model" & coor$response == "0", 
                           c(1,2)])
    # Position of the maximum density value
    max_pos_class0 <- which(den_class0$estimate == max(den_class0$estimate), 
                            arr.ind = TRUE)[1,]
    # Value of the coordinates at the point of maximum density
    class0 <- c(den_class0$eval.points[[1]][max_pos_class0[1]], 
                den_class0$eval.points[[2]][max_pos_class0[2]])
    
    den_class1 <- kde(coor[coor$group == "model" & coor$response == "1", 
                           c(1,2)])
    max_pos_class1 <- which(den_class1$estimate == max(den_class1$estimate), 
                            arr.ind = TRUE)[1,]
    class1 <- c(den_class1$eval.points[[1]][max_pos_class1[1]], 
                den_class1$eval.points[[2]][max_pos_class1[2]])
    
    maxima <- rbind(class0, class1)
    colnames(maxima) <- colnames(coor[,c(1,2)]) 
    
    return(maxima)
  }
  else { 
    stopifnot(coordinates %in% colnames(coor[,c(1,2)]))
    den_class0 <- stats::density(coor[coor$group == "model" & coor$response == "0", 
                               c(coordinates)], kernel = "epanechnikov")
    # Position of the maximum density value
    max_pos_class0 <- which.max(den_class0$y)
    # Value of the coordinates at the point of maximum density
    class0 <- den_class0$x[max_pos_class0] 
    
    den_class1 <- stats::density(coor[coor$group == "model" & coor$response == "1", 
                               c(coordinates)], kernel = "epanechnikov")
    max_pos_class1 <- which.max(den_class1$y)
    class1 <- den_class1$x[max_pos_class1] 
    
    maxima <- rbind(class0, class1)
    colnames(maxima) <- coordinates
    
    return(maxima)}
}

#' @title Assignment of the map-predicted class
#' @description 
#' It calculates for each individual in the validation group ("test") the 
#' class predicted by the map. It builds a data frame with the prediction 
#' made, the actual classification, whether these are matched and the 
#' score used for discrimination.
#' @param coor The data frame with coordinates of individuals. The penultimate 
#' column must correspond to the response variable. The last column must 
#' correspond to the group (training or test).It must be the output of 
#' the function "CalculateNewCoordinates" or "CalculatePolarCoordinates".
#' @param maxima_points A list from the function "MaximumDensityFunctions" 
#' with the maxima of the density functions of each class (class0 and
#' class1). 
#' @param coordinates Coordinates used in the discrimination. 
#' By default, they are two coordinates, corresponding to the first 
#' two columns of the "coor" data frame. If you only want to do it for one, 
#' you must indicate the name of the column containing that coordinate.
#' @return Data frame with the prediction made, the actual classification, 
#' and the score used for discrimination (the difference between the distance 
#' of the coordinates of each individual with each of the maxima of the 
#' two classes).
#' @export
ClassAssignment <- function(coor, maxima_points, coordinates = "two") {
  coor_test_set <- coor[coor$group == "test",]
  
  if (coordinates == "two") {
    # Differences with the maximums of each class
    distance_to_class0 <- sqrt((coor_test_set[,1] - maxima_points["class0",1])^2 
                               + (coor_test_set[,2] 
                                  - maxima_points["class0", 2])^2)  
    distance_to_class1 <- sqrt((coor_test_set[,1] - maxima_points["class1",1])^2 
                               + (coor_test_set[,2] 
                                  - maxima_points["class1", 2])^2)  
  }
  else {
    stopifnot(coordinates %in% colnames(coor[,c(1,2)]))
    distance_to_class0 <- abs(coor_test_set[,coordinates] 
                              - maxima_points["class0", coordinates]) 
    distance_to_class1 <- abs(coor_test_set[,coordinates] 
                              - maxima_points["class1", coordinates]) }
  
  difference <- distance_to_class0 - distance_to_class1
  # each individual is assigned to the class with the 
  # maximum density closest to its coordinates.
  prediction <- ifelse(difference < 0, 0, 1)
  prediction_dt <- data.frame(coor_test_set$response, prediction, difference)
  colnames(prediction_dt) <- c("real", "predicted", "difference")
  
  return(prediction_dt)
}

#' @title Calculation of the predictive parameters
#' @description 
#' It calculates from the data frame constructed by the function 
#' "ClassAssignment "the values of the confusion matrix in relative 
#' frequencies and other predictive parameters (predictive values, accuracy, 
#' F1 score, Cohen's Kappa and Matthews Correlation coefficients).
#' @param prediction_dt Data frame with the prediction made, the actual 
#' classification, and the score used for discrimination (the difference 
#' between the distance of the coordinates of each individual with each 
#' of the maxima of the two classes). It must be the output of 
#' the function "ClassAssignment".
#' @return List of calculated predictive parameters (relative confussion
#' matrix, negative and positive predictive values, F1 score, Cohen's
#' kappa and Matthew's correlation coefficients)
#' @export
PredictiveParameters <- function(prediction_dt) {
  TN <- nrow(prediction_dt[prediction_dt$real == 0 
                           & prediction_dt$predicted == 0, ])
  TP <- nrow(prediction_dt[prediction_dt$real == 1 
                           & prediction_dt$predicted == 1, ])
  FN <- nrow(prediction_dt[prediction_dt$real == 1 
                           & prediction_dt$predicted == 0, ])
  FP <- nrow(prediction_dt[prediction_dt$real == 0 
                           & prediction_dt$predicted == 1, ])
  TN_rel <- TN/(TN + FP)
  TP_rel <- TP/(TP + FN)
  FN_rel <- FN/(TP + FN)
  FP_rel <- FP/(TN + FP)
  
  NPV <- TN/(TN + FN)
  PPV <- TP/(TP + FP)
  
  accuracy <- (TN + TP)/(TN + TP + FN + FP)
  f1_score <- (2*TP)/((2*TP) + FN + FP)
  cohens_kappa <- 2*((TP*TN)-(FN*FP))/((TP+FP)*(FP+TN)+(TP+FN)*(FN+TN))
  mcc <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  return(list(confusion_matrix =(c(TN = TN_rel, FP = FP_rel, 
                                   FN = FN_rel, TP = TP_rel)),
              predictive_parameters = c(NPV = NPV, PPV = PPV, 
                                        accuracy = accuracy, 
                                        f1_score = f1_score, 
                                        cohens_kappa = cohens_kappa, 
                                        mcc = mcc)))
}

#' @title Categorical Data Classification Analysis
#' @description 
#' It performs the Categorical Data Classification Analysis from a 
#' data.frame of categorical data
#' @param db Database predictor and response variables. 
#' It must contain either character data or factors. The last 
#' column must contain the response variable
#' @param coordinates Coordinates used in the discrimination. 
#' By default, the two polar coordinates will be used ("two").
#' If you only want to do it for one, you must indicate the name of the 
#' coordinate ("theta" or "rho").
#' @param n_model Number of individuals to build the map. The number 
#' of individuals in the test set is calculated by difference.
#' @return CDCA_object. List of the results of the analysis: 
#' - preprocessed_db The database with the dummy variables. The antepenultimate 
#' and penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test) 
#' to which the individual has been assigned.
#' - base The base obtained after the singular value decomposition 
#' of the density matrix. Only vectors whose singular values are 
#' greater than 0 are retained.
#' - new_coordinates The data frame with the coordinates in the new basis. The 
#' penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test).
#' - polar_coor The data frame with the polar coordinates. The 
#' penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test).
#' - maxima Matrix with the maxima of the density functions of each 
#' class (rows) 
#' - prediction Data frame with the prediction made, the actual classification, 
#' and the score used for discrimination (the difference between the distance 
#' of the coordinates of each individual with each of the maxima of the 
#' two classes).
#' - results List of calculated predictive parameters (relative confussion
#' matrix, negative and positive predictive values, F1 score, Cohen's
#' kappa and Matthew's correlation coefficients). 
#' @import ks
#' @import ExPosition
#' @import dendextend
#' @export
CategoricalDataClassificationAnalysis <- function(db, n_model, 
                                                  coordinates = "two") {
  stopifnot(is.data.frame(db) & all((apply(db, 2, is.character)) 
                                    | all(apply(db, 2, is.factor))) 
            & coordinates %in% c("rho", "theta", "two"))
  stopifnot(dendextend::is.natural.number(n_model) & n_model <= nrow(db))
  
  db_dummy <- makeNominalData(db)
  preprocessed_db <- DataPreparation(db_dummy, n_model)
  base <- BaseConstruction(preprocessed_db)
  new_coordinates <- CalculateNewCoordinates(preprocessed_db, base)
  polar_coor <- CalculatePolarCoordinates(new_coordinates, preprocessed_db) 
  maxima  <- MaximumDensityFunctions(polar_coor, coordinates)
  prediction <- ClassAssignment(polar_coor, maxima, coordinates)
  results <- PredictiveParameters(prediction)
  
  CDCA_object <- list("preprocessed_db" = preprocessed_db,
                      "base" = base,
                      "new_coordinates" = new_coordinates,
                      "polar_coor" = polar_coor,
                      "maxima" = maxima,
                      "prediction" = prediction,
                      "results" = results)
  
  class(CDCA_object) <- "CDCA_object"
  
  return(CDCA_object)
}