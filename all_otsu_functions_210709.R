## --- packages to run multiotsu
require(reticulate)
skf = import("skimage.filters")

## --- otsu data collection
## based on EBImage::otsu()
get_otsu = function(y, breaks=1024){
  h = hist.default(y, breaks = breaks, plot = FALSE)
  counts = as.double(h$counts)
  mids = as.double(h$mids)
  len = length(counts)
  w1 = cumsum(counts)
  w2 = w1[len] + counts - w1
  cm = counts * mids
  m1 = cumsum(cm)
  m2 = m1[len] + cm - m1
  var = w1 * w2 * (m2/w2 - m1/w1)^2
  # find the left- and right-most maximum and return the threshold value in between
  maxi = which(var == max(var, na.rm = TRUE))
  (mids[maxi[1]] + mids[maxi[length(maxi)]] ) /2
}

get_otsu_sk = function(x,k=4){
  skf$threshold_otsu(image = np_array(x))
}

generate_otsu_data = function(sc,all_vars,methods,k=4){
  ## split data into slides
  slides = split(sc,factor(sc$SlideID))
  ## setup otsu data
  thresholds = data.frame()
  
  for(v in all_vars){
    ## split data by channel
    for(s in slides){
      ## get variable names
      all_iters = paste0(v,"_",methods)
      otsu_vals = unname(sapply(all_iters,function(v){
        get_otsu_sk(x=s[,v])
      }))
      ## get values
      thresholds = rbind(c(unique(s$SlideID),v,otsu_vals),thresholds)
    }
  }
  cols = paste0(methods, "_threshold")
  colnames(thresholds) = c('SlideID','channel',cols)
  for(v in grep("threshold",colnames(thresholds))){
    thresholds[,v] = as.numeric(thresholds[,v])
  } 
  return(thresholds)
}

## --- otsu mse
get_thresholds = function(all_vars,meth,dat=sc){
  channels = paste0(all_vars,"_",meth)
  apply(dat[,channels], 2, get_otsu_sk)
}

get_tdat_mse = function(otsus){
  ## convert otsus df into long format
  otsuLong = reshape(otsus, direction='wide',
                     v.names=colnames(otsus)[-c(1:2)],
                     timevar = 'channel', idvar='SlideID', sep='_')
  
  ## MSEs from Otsu applied to full dataset
  tdat = data.frame()
  for(m in methods){
    tdat = rbind(tdat, round(colMeans(sweep(otsuLong[,(paste0(m,"_threshold_",all_vars))], 2, get_thresholds(all_vars,m),FUN="-")^2),3))
  }
  
  rownames(tdat) = methods
  colnames(tdat) = all_vars
  
  return(tdat)
}

## --- otsu SDs
get_all_sds = function(otsus){
  ## sds by method
  all_sds= data.frame()
  for(m in methods){
    method_sds = otsus %>%
      group_by(channel) %>%
      summarize(!!m:=round(sd(get(paste0(m,"_threshold"))),3))
    all_sds = rbind(all_sds,t(method_sds[,m]))
  }
  
  colnames(all_sds) = all_vars
  
  return(all_sds)
}

## --- otsu misclassification
checkThreshold = function(vec,thr1,thr2){
  tab = caret::confusionMatrix(factor(vec>thr1),factor(vec>thr2),)$table
  
  return(tab/sum(tab) * (1-diag(2)))
}

get_misclassification = function(m, otsus, all_vars1=all_vars, dat=sc){
  slides = split(dat,factor(dat$SlideID))
  slideIDs = unique(dat$SlideID)
  thresholds = get_thresholds(all_vars1,meth=m)
  
  method_checks = data.frame()
  for(v in all_vars1){
    #print(v)
    all_checks= c()
    for(i in 1:length(slides)){
      ## get values for channel
      s = slides[[i]]
      vec = s[,paste0(v,"_",m)] ## data
      
      thr1 = thresholds[paste0(v,"_",m)]
      thr2 = otsus[otsus$SlideID==unique(s$SlideID) & otsus$channel==v,paste0(m,'_threshold')]
      
      all_checks = c(all_checks,sum(checkThreshold(vec,thr1,thr2)))
    }
    check_dat = data.frame(all_checks)
    check_dat$channel = v
    check_dat$method = m
    check_dat$SlideID = slideIDs
    method_checks = rbind(method_checks,check_dat)
  }
  return(method_checks)
}

run_misclassification = function(dat1,otsus1,methods1){
  all_method_checks = data.frame()
  for(m in methods1){
    #print(m)
    mdat = get_misclassification(m=m,dat=dat1,otsus=otsus1)
    all_method_checks = rbind(all_method_checks,mdat)
  }
  all_means = all_method_checks %>% group_by(channel,method) %>% summarise(mean_misclass=mean(all_checks))
  all_means = tidyr::spread(all_means,channel,mean_misclass)
  return(all_means)
}