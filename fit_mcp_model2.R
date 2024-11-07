fit_mcp_model2 <- function(data) {
  library(mcp)
  max_time <- max(data$Row_shift_scale, na.rm = TRUE)
  get_dynamic_priors <- function(data) {
    
    max_time <- max(data$Row_shift_scale, na.rm = TRUE)
    
    # Base priors
    priors <- list(
      int_1 = "dnorm(10, 70)",
      cp_1 = "dunif(0, 40)",
      cp_2 = "dunif(70, 120)",
      cp_3 = "dunif(90, 130)",
      cp_4 = "dunif(150, 250)",
      Row_shift_scale_1 = "dunif(-0.5, 1)", 
      Row_shift_scale_2 = "dunif(1, 5)",
      Row_shift_scale_3 = "dunif(-1, -5)",
      Row_shift_scale_4 = "dunif(-0.5, 1)",
      Row_shift_scale_5 = "dunif(1, 5)",
      Row_shift_scale_6 = "dunif(-0.5, 0.8)"
    )
    
    # Conditionally include cp_5 and Row_shift_scale_7
    if (max_time > 350 && max_time <= 450) {
      priors$cp_5 <- sprintf("dunif(350, %s)", max_time)
      priors$Row_shift_scale_7 <- "dunif(-5, -1)"
    } else if (max_time > 450) {
      priors$cp_5 <- "dunif(350, 450)"
      priors$Row_shift_scale_7 <- "dunif(-5, -1)"
    }
    
    return(priors)
  }
  
  priors <- get_dynamic_priors(data)
  
  get_dynamic_model <- function(max_time) {
    
    base_model <- RescaledIntensity ~ 1 + Row_shift_scale
    
    additional_model_components <- list(
      ~ 0 + Row_shift_scale,
      ~ 0 + Row_shift_scale,
      ~ 0 + Row_shift_scale,
      ~ 0 + Row_shift_scale
    )
    
    if (max_time > 350) {
      additional_model_components <- append(additional_model_components, list(~ 0 + Row_shift_scale))
      additional_model_components <- append(additional_model_components, list(~ 0 + Row_shift_scale))
    }
    
    model <- c(list(base_model), additional_model_components)
    
    return(model)
  }
  
  
  model <- get_dynamic_model(max_time = max_time)
  
  fit = mcp(model, data = data, prior = priors)
  # Consider saving the plots with ggsave() or similar
  
  return(summary(fit))
}
