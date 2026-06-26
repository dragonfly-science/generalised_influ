#' Plot Residual Implied Index
#'
#' @param fit A fitted model object
#' @param grouping_var A character string specifying the column name to group by. 
#' Defaults to 'stat_area'.
#' @importFrom splines ns
#' @importFrom cowplot theme_cowplot
#' @return A ggplot2 object showing the implied vs. scaled indices.
#' @export


plot_RIC <- function(fit, grouping_var = 'stat_area', min_records = 10,  add.rho = TRUE, custom_palette =  default_palette){
  year <- get_first_term(fit)
  grouping_var <- sym(grouping_var)
  
  component <- if(length(fit$family$family)==2) {
    'Combined'
  } else if (!is.null(fit$family$family) && (fit$family$family %in% c("bernoulli", "binomial"))){
    'Binomial'
  } else {
    'Positive'
  }
  
  
  idx <- get_index(fit, format = "wide")
  idx <- idx %>%
    filter(Index==component)
  
  
  raw_data <- fit$data
  
  
  if (is.null(raw_data) && inherits(fit, "survreg")) {
    raw_data <- eval(fit$call$data)
  }

  if(inherits(fit, "survreg")) logobs <- log(fit$model[,1][,1]) else if (grepl("log", formula(fit)[2])) logobs <- fit$model[,1] else logobs <- log(fit$model[,1])
  if (grepl("log", formula(fit)[2])) logpred <- predict(fit) else logpred <- log(predict(fit))

  #browser()
  ric_data <- tibble(level = raw_data[[year]],
                     !!grouping_var := raw_data[[as.character(grouping_var)]],                              
                     resid  = logobs - logpred) %>%
    add_count(level, !!grouping_var) %>% 
    filter(n >= min_records) %>%
    dplyr::left_join(idx %>% select(level, stan, stan_unscaled)) %>%
    mutate(implied = exp(resid + log(stan_unscaled))) %>%
    group_by(!!grouping_var) %>%    
    mutate(base_imp = {
      # Calculate arithmetic mean for each year 
      yearly_means <- tapply(implied, level, mean, na.rm = TRUE)
      # Calculate geometric mean of those yearly means
      exp(mean(log(yearly_means), na.rm = TRUE))
    },
    imp_scaled = implied / base_imp,
    idx_scaled = stan,
    stan_base = gmean(unique(stan))) %>%
    group_by(!!grouping_var, level) %>%
    summarise(n = n(),
              se = sd(imp_scaled)/ sqrt(n()),
              implied = mean(implied),
              idx = mean(stan_unscaled),
              imp_scaled = mean(imp_scaled)*unique(stan_base), # rescale to geo mean of stan base during overlapping perdiod
              idx_scaled = mean(idx_scaled)) %>%
    mutate(level  = as.integer(as.character(level))) %>%
    complete(level = full_seq(level, 1)) %>% 
    ungroup()
    
  imp_count <- ric_data %>%
    group_by(!!grouping_var) %>%
    summarise(Num = sum(n, na.rm = T),
           rho=cor(implied, idx_scaled, use="pairwise.complete.obs"))
  
  
  
  
  p <-   ggplot(ric_data,
                aes(x=level,
                    y=imp_scaled))+
    geom_point(aes(size=n, colour = "Implied index"),  alpha = 0.3)+
    geom_line(aes(colour = "Implied index"))+
    geom_errorbar(aes(ymin=(imp_scaled-1.96*se),
                      ymax=(imp_scaled+1.96*se),
                    colour = "Implied index"),
                  linewidth=0.3,
                  width=0.3)+
    geom_hline(yintercept=1,
               linetype=3,
               color = custom_palette['main'])+
    geom_line(aes(y=idx_scaled, group = 1,  colour = 'Standardised index'))+
    facet_wrap(as.formula(paste("~", as.character(grouping_var))),
               ncol=2,scales='free_y')+
    labs(x='Fishing year',
         y='Index',
         size="Records",
        colour=NULL) +
    guides(size = guide_legend(override.aes = list(colour = custom_palette['current']))) +
    scale_colour_manual(values = c("Implied index" = custom_palette[['current']], 
                                 "Standardised index" = custom_palette[['previous']])) +
     scale_y_continuous(limits = function(x) c(0, max(pretty(x))), 
                        expand = c(0, 0))+
    scale_x_continuous(breaks = function(x) {
      # x is the range of years present in the data
      years <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
      
      if (length(years) > 15) {
        # If more than 15 years, show every 2nd year
        round(seq(min(years), max(years), by = 2))
      } else {
        # Otherwise show all
        round(years)
      }
    }) +
    
    theme_cowplot()+
    theme(axis.text.x = element_text(hjust = 0, angle = 90, size = 12),
          legend.position = "top",
          legend.justification = "right",
          legend.direction = "horizontal", 
          legend.box = "vertical",
          legend.box.just = "right",
          legend.key.spacing.x = unit(1.5, 'lines'),
          legend.spacing.y = unit(0, 'mm'),
          panel.border = element_blank(),
          panel.spacing = unit(30, "pt")) +
    (if(add.rho) {
      list(
          geom_text(data = imp_count, aes(x = Inf, y = Inf, label = paste0("N = ", Num)), 
                    vjust = -3.2, hjust = 1.1, colour = "#03576E"),
          geom_text(data = imp_count, aes(x = Inf, y = Inf, label = paste0('rho == ', round(rho, 2))), 
                    vjust = -0.7, hjust = 1.1, colour = "#03576E", parse = TRUE),
          coord_cartesian(clip = "off"),
          theme( strip.text = element_text(margin = margin(t = 15, b = 20, unit = "pt") ))
  )
       
    } else {
      NULL
    })
  
  return(p)
}

