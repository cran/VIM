## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.align = "center"
)

## ----setup, message = FALSE---------------------------------------------------
library(VIM)
library(magrittr)
dataset <- sleep[, c("Dream", "NonD", "BodyWgt", "Span")]
dataset$BodyWgt <- log(dataset$BodyWgt)
dataset$Span <- log(dataset$Span)
aggr(dataset)
str(dataset)

## -----------------------------------------------------------------------------
imp_regression <- regressionImp(NonD ~ BodyWgt + Span, dataset)
imp_ranger <- rangerImpute(NonD ~ BodyWgt + Span, dataset)
aggr(imp_regression, delimiter = "_imp")

## ---- fig.height=5------------------------------------------------------------
imp_regression[, c("NonD", "BodyWgt", "NonD_imp")] %>% 
  marginplot(delimiter = "_imp")

## ---- fig.height=5------------------------------------------------------------
imp_ranger[, c("NonD", "BodyWgt", "NonD_imp")] %>% 
  marginplot(delimiter = "_imp")
imp_ranger[, c("NonD", "Span", "NonD_imp")] %>% 
  marginplot(delimiter = "_imp")

## -----------------------------------------------------------------------------
imp_regression <- regressionImp(Dream + NonD ~ BodyWgt + Span, dataset)
imp_ranger <- rangerImpute(Dream + NonD ~ BodyWgt + Span, dataset)
aggr(imp_regression, delimiter = "_imp")

