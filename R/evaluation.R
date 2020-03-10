#' Error performance measures
#' 
#' Various error measures evaluating the quality of imputations
#'  
#' @name evaluation
#' @rdname evaluation
#' @aliases evaluation nrmse pfc msecov msecor
#' @param x matrix or data frame
#' @param y matrix or data frame of the same size as x 
#' @param m the indicator matrix for missing cells
#' @param vartypes a vector of length ncol(x) specifying the variables types, like factor or numeric
#' @return the error measures value
#' @author Matthias Templ
#' @references M. Templ, A. Kowarik, P. Filzmoser (2011) Iterative stepwise
#' regression imputation using standard and robust methods.  \emph{Journal of
#' Computational Statistics and Data Analysis}, Vol. 55, pp. 2793-2806.
#' 
# seealso \code{\link{robCompositions::rdcm}}
#' @details This function has been mainly written for procudures 
#' that evaluate imputation or replacement of rounded zeros. The ni parameter can thus, e.g. be
#' used for expressing the number of rounded zeros.
#' @keywords manip
#' @export
#' @examples
#' data(iris)
#' iris_orig <- iris_imp <- iris
#' iris_imp$Sepal.Length[sample(1:nrow(iris), 10)] <- NA
#' iris_imp$Sepal.Width[sample(1:nrow(iris), 10)] <- NA
#' iris_imp$Species[sample(1:nrow(iris), 10)] <- NA
#' m <- is.na(iris_imp)
#' iris_imp <- VIM::kNN(iris_imp, imp_var = FALSE)
#' \dontrun{
#' evaluation(iris_orig, iris_imp, m = m, vartypes = c(rep("numeric", 4), "factor"))
#' msecov(iris_orig[, 1:4], iris_imp[, 1:4])
# nrmse <- function(x, y, m){
#   return(sqrt( (sum((x[m] - y[m])^2) / sum(m)) / var(x[m])) )
# }
#' }
evaluation <- function(x, y, m, vartypes = "guess"){
  guessType <- function(x){
    type <- as.character(unlist(lapply(x, function(x) class(x)[1])))
    indexNumeric <- type == "numeric" | type =="integer"
    indexNominal <- type == "factor"
    indexOrdered <- type == "ordered"
    indexMixed <- rep(FALSE, length(type))
    indexBinary <- rep(FALSE, length(type))
    
    if(any(type == "character")){
      indexNominal <- type == "character" | indexNominal
    }
    
    checkMixed <- function(x, constant = 0.2){
      if(any(table(x) > constant * length(x)) & length(table(x)) > length(x) * 1/10){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    
    checkNominal <- function(x, constant = 0.2){
      if(length(table(x)) < length(x) * 1/10 & !any(x[!is.na(x)] %% 1 > .Machine$double.eps)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    
    checkBinary <- function(x){
      if(length(table(x)) == 2){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    
    checkMixedConst <- function(x, constant = 0.2){
      return(as.numeric(names(which(table(x) > constant * length(x)))))
    }
    
    # if numeric
    #   is it binary?
    #   is it nominal?
    #   is it mixed?
    # if nominal
    #   is it binary?
    
    if(any(indexNumeric)){ 
      for(i in which(indexNumeric)){
        if(checkBinary(x[, i])){
          indexBinary[i] <- TRUE
          indexNumeric[i] <- FALSE    
        }
        if(checkNominal(x[, i]) & !checkBinary(x[, i])){
          indexNominal[i] <- TRUE  
          indexNumeric[i] <- FALSE   
        }
        if(checkMixed(x[, i] & !checkBinary(x[, i]) & !checkNominal(x[, i]))){
          indexMixed[i] <- TRUE  
          indexNumeric[i] <- FALSE           
        }
      }
    }
    
    if(any(indexNominal)){
      for(i in which(indexNominal)){   
        if(checkBinary(x[, i])){
          indexBinary[i] <- TRUE
          indexNominal[i] <- FALSE    
        }      
      }
    }
    
    m <- rbind(indexBinary, indexNominal, indexOrdered, indexNumeric, indexMixed)
    colnames(m) <- colnames(x)
    type <- numeric(ncol(x))
    for(i in 1:ncol(x)){
      type[i] <- names(which(m[, i]))
    }
    type <- tolower(substr(type, 6, nchar(type)))
    res <- list("indices" = m,
                "type" = type)
    return(res)
  }
  if(vartypes == "guess"){
    vartypes <- guessType(x)$type
  }
  err_num <- err_cat <- err_mixed <- 0
  if(any(vartypes == "numeric")){
    # err_num <- sum((x[, vartypes == "numeric"] - y[, vartypes == "numeric"])^2) / sum(m[, vartypes == "numeric"])
    xx <- x[, vartypes == "numeric"]
    yy <- y[, vartypes == "numeric"]
    mm <- m[, vartypes == "numeric"]
    err_num <- median(abs(xx[mm] - yy[mm]))
  }
  if(any(vartypes == "factor")){
    err_cat <- sum(x[, vartypes == "factor"] != y[, vartypes == "factor"]) / sum(m[, vartypes == "factor"])
  } 
  results <- list("err_num" = err_num,
                  "err_cat" = err_cat,
                  "error" = err_num + err_cat + err_mixed)
  return(results)
}

#' @rdname evaluation
#' @export

nrmse <- function(x, y, m){
  return(sqrt( mean((x[m] - y[m])^2)  / var(x[m])) )
}
# nrmse <- function(x, y, m){
#   bias <- x[m] - y[m]
#   variance <- var(x[m] - y[m]) / var(x[m])
#   variance
# }

#' @rdname evaluation
#' @export

pfc <- function(x, y, m){
  return(sum(x != y) / sum(m))
}

#' @rdname evaluation
#' @export

mape <- function(x, y, m){
  return(median(abs(x[m] - y[m])^2) * 100)
}

#' @rdname evaluation
#' @export

lr <- function(x, y, m){
  return(median(abs(x[m] / y[m])))
}

#' @rdname evaluation
#' @export

smape <- function(x, y, m){
  return((abs(x[m] - y[m]))/((abs(x[m]) + abs(y[m]))))
}

# ced <- function(x, y, m){
#   return(robCompositions::aDist(x, y) / sum(m))
# }

#' @rdname evaluation
#' @export

msecov <- function(x, y){
  sum((cov(x) - cov(y))^2) / ncol(x)
}

#' @rdname evaluation
#' @export

msecor <- function(x, y){
  sum((cor(x) - cor(y))^2) / ncol(x)
}

