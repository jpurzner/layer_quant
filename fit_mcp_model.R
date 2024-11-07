fit_mcp_model <- function(data) {
  library(mcp)
  
  max_time <- max(data$Row_shift_scale, na.rm = TRUE)
  
  # Base model and priors
  model_list <- list(RescaledIntensity ~ 1 + Row_shift_scale, # Linear segment spans the iEGL
                     ~ 0 + Row_shift_scale, # Linear without intercept start of oEGL, positive slope, start cp_1 and end cp_2 
                     ~ 0 + Row_shift_scale, # Linear without intercept, peak of oEGL to end, negative slope, start cp_2, end cp_3
                     ~ 0 + Row_shift_scale, # Relatively Flat (no slope), molecular layer, start cp_3 end cp_4 
                     ~ 0 + Row_shift_scale, # Linear without intercept, start of IGL, start cp_4 and cp_5
                     ~ 0 + Row_shift_scale) # Linear without intercept, flat part of the IGL
  
  priors <- list(
    int_1 = "dnorm(10, 70)", # normal prior with mean=0, precision=
    cp_1 = "dunif(0, 30)",  #start of iegl
    cp_2 = "dunif(50, 100)", #peak of iegl
    cp_3 = "dunif(100, 130)", # start of molecular layer
    cp_4 = "dunif(150, 250)", # start of igl raise
    cp_5 = if (max_time < 300) sprintf("dunif(225, %s)", max_time-10) else "dunif(225, 300)", # start of igl flat part
    Row_shift_scale_1 = "dunif(-0.5, 1)", # flat slope of the oEGL   
    Row_shift_scale_2 = "dunif(1, 5)", # start of oEGL, positive slope
    Row_shift_scale_3 = "dunif(-5, -0.5)", # peak of oEGL to end, negative slope, 
    Row_shift_scale_4 = "dunif(-0.5, 1)", # molecular layer, relatively flat  
   Row_shift_scale_5 = "dunif(1, 5)", #start of the IGL, positive slope 
   Row_shift_scale_6 = "dunif(-0.5, 0.8)" # within the IGL, relatively flat 
  )
  
  

  
  # Add additional segments if max_time conditions are met
  if (max_time > 350) {
   model_list <- c(model_list, list(~ 0 + Row_shift_scale))
   priors$Row_shift_scale_7 <- "dunif(-5, -1)"
   priors$cp_6 <- if (max_time <= 450) sprintf("dunif(350, %s)", max_time-10) else "dunif(350, 450)"
 
  }
  
  model <- do.call(c, model_list)
  fit = mcp(model, data = data, prior = priors)
  
  #fit = mcp(model_list, data = data)
  
  
  #return(model)
  return(summary(fit))
}
