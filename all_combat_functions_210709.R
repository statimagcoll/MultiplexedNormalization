### --- Script to adjust raw data with ComBat algorithm

### -------SETUP-------
require(dplyr)
require(tidyr)

## ---- COMBAT functions

## data adjustment function
remove_single_images = function(chan, image_var){
  ## count cells by images
  sub_chan = chan %>% group_by_at(image_var) %>% count()
  sub_chan$bool = sub_chan$n <= 1
  
  ## mark cells that are n-of-1 in an image
  nof1s = sub_chan[sub_chan$bool == TRUE,image_var]
  
  ## return dataset with no N-of-1s
  return(chan[!(chan$Pos %in% nof1s$Pos),])
}

## internal function for delta functions
sqerr = function(x){sum((x - mean(x))^2)}

## update each iteration of the algo
update_gamma = function(batch_chan, gamma_c, tau_c,channel, slide_var){
  ## create numerator value
  #batch_chan$gamma_num = (batch_chan[,channel] - batch_chan$alpha_c)/batch_chan$delta_ijc
  batch_chan$gamma_num = batch_chan[,channel]
  
  countr = batch_chan %>%
    group_by_at(slide_var) %>%
    count()
  
  gamma_num = batch_chan %>%
    group_by_at(slide_var) %>%
    summarise(avg = mean(gamma_num),.groups='drop')
  
  gamma_num$avg = (countr$n * tau_c * gamma_num$avg) + gamma_c * unique(batch_chan$delta_ijc)
  
  ## create denominator value
  # gamma_denom = batch_chan %>%
  #   group_by_at(slide_var) %>%
  #   summarise(avg = mean(delta_ijc_inv),.groups='drop')
  #gamma_denom$avg = gamma_denom$avg + (1/tau_c)
  gamma_denom = countr$n * tau_c + unique(batch_chan$delta_ijc)
  
  gamma_ic_star = gamma_num
  gamma_ic_star$avg = gamma_ic_star$avg / gamma_denom
  #gamma_ic_star$avg = gamma_ic_star$avg/gamma_denom$avg 
  
  ## returns zero if only one slide
  #if(is.na(gamma_ic_star$avg[1])){gamma_ic_star$avg<-0}
  return(gamma_ic_star)
}
update_delta2 = function(batch_chan, beta_c,omega_c,channel,slide_var){
  batch_chan$delta_vals = batch_chan[,channel] #- batch_chan$gamma_ic
  
  #- batch_chan$lambda_ijc)
  delta_num = batch_chan %>%
    group_by_at(slide_var) %>%
    #summarise(avg=mean(delta_vals))
    summarise(avg = sum((delta_vals - gamma_ic)^2),.groups='drop')
    #summarise(avg = sqerr(delta_num),.groups='drop')
  
  delta_denom = batch_chan %>%
    group_by_at(slide_var) %>%
    count()
  
  delta_denom$n = delta_denom$n/2 + omega_c - 1
  delta_num$avg = 0.5*delta_num$avg + beta_c
  
  delta_ijc_star = delta_num
  delta_ijc_star$avg = delta_ijc_star$avg/delta_denom$n
  
  #delta_ijc_star[is.na(delta_ijc_star$avg),]$avg = 0.00001
  
  return(delta_ijc_star)
}

## checking convergence
gamma_conv = function(batch_chan, gamma_stars,slide_var){
  gams = batch_chan[,c(slide_var,'gamma_ic')] %>% distinct()
  return(mean(abs(gams[match(unlist(gamma_stars[,slide_var]),
                             gams[,slide_var]),]$gamma_ic - gamma_stars$avg))) ## MAE
}
delta_conv = function(batch_chan, delta_stars,slide_var){
  dels = batch_chan[,c(slide_var,'delta_ijc')] %>% distinct()
  return(mean(abs(dels[match(unlist(delta_stars[,slide_var]),
                             dels[,slide_var]),]$delta_ijc - delta_stars$avg))) ## MAE
}

## function to combat-adjust for one channel
adjust_vals = function(channel,slide_var,chan,remove_zeroes=TRUE,
                       tol = 0.0001){
  if(remove_zeroes){
    ## remove zeroes if needed
    leftover = chan[chan[,channel] <=0,]
    chan = chan[(chan[,channel] > 0),]
    
  }
  
  ### -------COMBAT EMPIRICAL VALUES-------
  
  ## get alpha (grand mean)
  chan$alpha_c = mean(chan[,channel])
  
  pre_gamma_ic = chan %>% 
    group_by_at(slide_var) %>% 
    summarise(avg=mean(get(channel)), .groups = 'drop')
  chan$pre_gamma_ic = pre_gamma_ic[match(chan[,slide_var],unlist(pre_gamma_ic[,slide_var])),]$avg
  chan$pre_gamma_ic = chan$pre_gamma_ic #- chan$alpha_c
  
  chan[,paste0("Adj_",channel)] = chan[,channel] - chan$alpha_c - chan$pre_gamma_ic
  sigma_c = var(chan[,paste0("Adj_",channel)])
  chan[,channel] = (chan[,channel] - mean(chan[,channel]))/sigma_c
  
  ## get gammas (slide means)
  gamma_ic = chan %>% 
    group_by_at(slide_var) %>% 
    summarise(avg=mean(get(channel)), .groups = 'drop')
  
  chan$gamma_ic = gamma_ic[match(chan[,slide_var],unlist(gamma_ic[,slide_var])),]$avg
  chan$gamma_ic = chan$gamma_ic 
  
  ## get deltas (slide variances)
  chan$delta_ijc = chan[,channel]
  
  delta_ijc = chan %>%
    group_by_at(slide_var) %>%
    summarise(v=sum((get(channel) - gamma_ic)^2), .groups='drop')
  
  chan$delta_ijc = (delta_ijc[match(chan[,slide_var],unlist(delta_ijc[,slide_var])),]$v)
  
  ### -------COMBAT HYPERPARAMETERS-------
  ## slide level mean
  gamma_c = mean(chan$gamma_ic)
  tau_c = var(chan$gamma_ic)
  
  ## slide level variances
  M_c = mean(chan$delta_ijc)
  S_c = var(chan$delta_ijc)
  
  ## is this correct?
  omega_c = (M_c + 2*S_c)/S_c
  beta_c = (M_c^3 + M_c*S_c)/S_c
  
  ### -------CALLING COMBAT BATCH EFFECTS FUNCTIONS-------
  batch_chan = chan ## duplicate the dataframe to iterate
  batch_chan$delta_ijc_inv = 1/batch_chan$delta_ijc
  
  gamma_c; tau_c
  M_c; S_c
  omega_c; beta_c
  
  ### -------COMBAT BATCH EFFECT ADJUSTMENT-------
  
  ## run a single iteration
  ## run delta first
  delta_stars = update_delta2(batch_chan, beta_c, omega_c,channel,slide_var)
  check_delta_conv = delta_conv(batch_chan, delta_stars,slide_var)
  batch_chan$delta_ijc = (delta_stars[match(batch_chan[,slide_var],unlist(delta_stars[,slide_var])),]$avg)
  batch_chan$delta_ijc_inv = 1/batch_chan$delta_ijc
  
  ## now update gamma
  gamma_stars = update_gamma(batch_chan, gamma_c, tau_c,channel,slide_var=slide_var)
  check_gamma_conv = gamma_conv(batch_chan, gamma_stars,slide_var=slide_var)
  batch_chan$gamma_ic = gamma_stars[match(batch_chan[,slide_var],unlist(gamma_stars[,slide_var])),]$avg
  
  total_mae = sum(check_gamma_conv,check_delta_conv)
  iterations = 1
  ## first check of MAE
  #print(paste0('Total MAE after ', iterations,' iterations: ', round(total_mae,8)))
  
  ## run until convergence 
  while(total_mae > tol){ 
    ## run delta first
    delta_stars = update_delta2(batch_chan, beta_c, omega_c,channel,slide_var)
    check_delta_conv = delta_conv(batch_chan, delta_stars,slide_var)
    batch_chan$delta_ijc = (delta_stars[match(batch_chan[,slide_var],unlist(delta_stars[,slide_var])),]$avg)
    batch_chan$delta_ijc_inv = 1/batch_chan$delta_ijc
    
    ## now update gamma
    gamma_stars = update_gamma(batch_chan, gamma_c, tau_c,channel,slide_var=slide_var)
    check_gamma_conv = gamma_conv(batch_chan, gamma_stars,slide_var=slide_var)
    batch_chan$gamma_ic = gamma_stars[match(batch_chan[,slide_var],unlist(gamma_stars[,slide_var])),]$avg
    
    total_mae = sum(check_gamma_conv,check_delta_conv)
    iterations = iterations + 1
    ## final check of MAE
    #print(paste0('Total MAE after ', iterations,' iterations: ', round(total_mae,4)))
  }
  
  ### -------COMBAT BATCH EFFECT RESULTS-------
  
  ## now adjust for the batch effects
  batch_chan$Y_ijc_star = (sigma_c/batch_chan$delta_ijc)*(batch_chan[,channel] - batch_chan$gamma_ic) + batch_chan$alpha_c
  
  ## add zeroes back in if needed
  if(remove_zeroes){
    ## add back in zeroes
    leftover$Y_ijc_star = 0
    leftover[,colnames(batch_chan)[!(colnames(batch_chan) %in% colnames(leftover))]] = NA
    batch_chan = rbind(batch_chan,leftover)
  }
  return(batch_chan$Y_ijc_star)
}