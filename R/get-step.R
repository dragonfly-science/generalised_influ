#' Get the standardised indices with terms consecutively added
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
get_step <- function(fit, indices = indices, year = NULL, do_plot = FALSE, pred_grid = NULL, predictor = NULL, ...) {
  
  if  (!inherits(fit, c("sdmTMB", "glm", "survreg", "brmsfit"))) stop("This model class is not supported.")
  
  year <- get_first_term(fit = fit)
  
  is_sdm <- inherits(fit, 'sdmTMB')
  
  if (is_sdm){
    
    if  (is.null(predictor)) stop("Argument 'predictor' is missing. Please specify 1 for the first part or 2 for the hurdle part.")
    
    #extract formula for this predictor
    newFormula <- (fit$formula)[[predictor]]
    
    # extract terms from this formula
    terms <- stats::terms(newFormula)
    terms_labels <- attr(terms, "term.labels")
    
    # initiate effects list
    effects <- list()
    
    # Create models with terms successively added
    for(termCount in 0:(length(terms_labels)+2)){
      
      if(termCount>0){
        
        term <- terms_labels[termCount]
        
        #Update both formula and model
        
        newFormula <- update.formula(newFormula,
                                     formula(paste("~ 0 + ",
                                                   paste(paste(terms_labels[1:min(termCount,length(terms_labels)) ],
                                                               collapse='+')))))
        
        
        
        # Turn spatial and/or spatio-temporal component on
        spatial <- ifelse(termCount<=length(terms_labels), 'off', fit$spatial[predictor])
        spatiotemp <- ifelse(termCount<=length(terms_labels)+1, 'off', fit$spatiotemporal[predictor])
        
        # refit the model
        fit_reduced <- update(fit, formula = newFormula, evaluate = TRUE, spatial = spatial, spatiotemporal = spatiotemp)
        
        # Get index for this model
        idx_reduced <- get_index (fit_reduced,  pred_grid = pred_grid, predictor = predictor)
        print(summary(fit_reduced))
        # Generate the right hand side of formula as name for index
        idx_name <- case_when(
          termCount == 1                           ~ term,
          termCount == length(terms_labels) + 2    ~ "+ spatiotemporal",
          termCount == length(terms_labels) + 1    ~ "+ spatial",
          TRUE                                     ~ paste("+", term)
        )
        
        # Store column of indices
        effects[[idx_name]] <- idx_reduced$stan
        
      } else {
        term = 'intercept'
        fit_reduced = update(fit,.~1)
      }
      
      # TO DO calculate summary statistics here
    }
    #  Combine all the effect columns into one wide data frame
    all_idx <- do.call(cbind, effects)
    
    #  Final bind indices from last iteration with unstan and CI + all the steps 
    indices <- cbind(idx_reduced, all_idx)    
    rownames(indices) <- NULL       
  }
  return(indices)
}
      
      
      