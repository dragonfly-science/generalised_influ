#' Get the standardised indices
#' 
#' Get the standardised indices each year with associated uncertainty and return either as a table or a ggplot.
#' 
#' @param fit An object of class \code{brmsfit}.
#' @param year The year or time label (e.g. year, Year, fishing_year, etc).
#' @param probs The percentiles to be computed by the \code{quantile} function.
#' @param rescale How to re-scale the series. Choose from "raw" to retain the raw series, "unstandardised" to re-scale to the geometric mean of the unstandardised series, or a number to re-scale by. 
#' @param do_plot Return a \code{ggplot} object instead of a \code{data.frame}.
#' @param ... Additional parameters passed to \code{fitted}.
#' @return a \code{data.frame} or a \code{ggplot} object.
#' @importFrom stats fitted
#' @importFrom brms is.brmsfit
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @export
#' 
get_index <- function(fit, year = NULL, do_plot = FALSE, pred_grid = NULL, predictor = NULL, ...) {
  
  if  (!inherits(fit, c("sdmTMB", "glm", "survreg", "brmsfit"))) stop("This model class is not supported.")
  
  is_sdm <- inherits(fit, 'sdmTMB')
  
  if (is_sdm){
    
    #extract formula for this predictor
    newFormula <- (fit$formula)[[predictor]]
    
    #extract tems from this formula
    terms <- stats::terms(newFormula)
    terms_labels <- attr(terms, "term.labels")
    
    for(termCount in 0:length(terms_labels)){
      
      if(termCount>0){
        
        term = terms_labels[termCount]
        
        #Update both formula and model
        
        newFormula = update.formula(newFormula,
                                    formula(paste("~ 0 + ",
                                                  paste(paste(terms_labels[1:termCount],
                                                              collapse='+')))))
        
        model = update(fit, formula = newFormula, evaluate = TRUE)
        
        # STOPPED HERE ON fEB 11
        
        #Get index for this model
        index = rep(NA, nrow(.$indices))
        index_index = 1:length(index) %w/o% .$excl
        index[index_index] = .$effects(model,excl=.$excl)
        
        #Add column to .$indices
        .$indices = cbind(.$indices,index)
        
        #Give index column the right hand side of formula as name
        names(.$indices)[ncol(.$indices)] = if(termCount==1) term else paste('+',term)
      } else {
        term = 'intercept'
        model = update(.$model,.~1)
      }
      
      type = class(model)[1]
      logLike =  switch(type,
                        survreg = model$loglik[2],
                        logLik(model)
      )
      fitted = switch(type,
                      survreg = predict(model,type='response'),
                      fitted(model)
      )
      
      #Sums-of-squared based R2 based on log(observed) and log(fitted)
      if(termCount==0) r2 = 0
      else r2 = cor(log(observed),
                    log(fitted))^2
      
      #Deviance pseudo-R2
      r2Dev = (model$null.deviance-model$deviance)/model$null.deviance
      if(length(r2Dev)==0)
        r2Dev = NA
      
      #Negelkerke pseudo-R2
      if(termCount==0)
        logLikeInterceptOnly = logLike
      
      n = length(observed)
      r2Negel = (1-exp((logLikeInterceptOnly-logLike)*(2/n)))/(1-exp(logLikeInterceptOnly*(2/n)))
      
      .$summary = rbind(.$summary,data.frame(term = term,
                                             k = length(coef(model)),
                                             logLike = logLike,
                                             aic = extractAIC(model)[2],
                                             r2 = r2,
                                             r2Dev = r2Dev,
                                             r2Negel = r2Negel)
      )
    }
    
  }
}