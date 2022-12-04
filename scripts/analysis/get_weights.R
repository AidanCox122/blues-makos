
# create a function to calculate AIC weights for different models
get_weight <- 
  # function will accept inputs where base is a vector of established predictors,
  # and test will be a vector of predictors to include in the next iteration of forward selection
  function(base = NULL, test) {
    # create a repository
    models <- vector(
      mode = 'list',
      length = length(test)) %>% 
      set_names(unique(test))
    # find AIC for each test model
    for(x in names(models)) {
      form <- 
        paste('clus2~', paste0(base,
                               collapse = '',
                               sep = '+'),
              x,
              sep = '')
      
      m1 <- 
        mblogit(formula = as.formula(form),
                random = ~1|ptt,
                data = combo_mod)
      
      models[x] <- 
        AIC(m1)}
    
    # now create a reference model with only the base variables
    if(is.null(base)){
      # if base is null, compare to a model with only a random effect of individual
      m2 <- 
        mblogit(formula = clus2 ~ 1,
                random = ~1|ptt,
                data = combo_mod)
      # store the results in the models list
      null <- length(test) + 1
      models[null] <- 
        AIC(m2)
    } else{
      # otherwise, if base is not null, build a new formula
      form2 <- 
        paste('clus2~', paste0(base,
                               collapse = '+'),
              sep = '')
      m2 <- 
        mblogit(formula = as.formula(form2),
                random = ~1|ptt,
                data = combo_mod)
      null <- length(test) + 1
      models[null] <- 
        AIC(m2)
    }
    
    # calculate AIC weights
    vec_AIC <- unlist(models)
    dAIC <- vec_AIC - min(vec_AIC)
    AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
    return(AICw)
  }