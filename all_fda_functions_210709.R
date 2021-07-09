require(fda)

## function to calculate weighted mean
weighted.mean.fd = function (x, w, ...) 
{
  if (!inherits(x, "fd")) 
    stop("'x' is not of class 'fd'")
  coef <- x$coefs
  coefd <- dim(coef)
  ndim <- length(coefd)
  basisobj <- x$basis
  nbasis <- basisobj$nbasis
  dropind <- basisobj$dropind
  ndropind <- length(dropind)
  if (ndim == 2) {
    coefmean <- matrix(apply(coef, 1, weighted.mean, w=w), nbasis - ndropind, 
                       1)
    coefnames <- list(dimnames(coef)[[1]], "Mean")
  }
  else {
    nvar <- coefd[3]
    coefmean <- array(0, c(coefd[1], 1, nvar))
    for (j in 1:nvar) coefmean[, 1, j] <- apply(coef[, , 
                                                     j], 1, weighted.mean, w=w)
    coefnames <- list(dimnames(coef)[[1]], "Mean", dimnames(coef)[[3]])
  }
  fdnames <- x$fdnames
  fdnames[[2]] <- "mean"
  fdnames[[3]] <- paste("mean", fdnames[[3]])
  meanfd <- fd(coefmean, basisobj, fdnames)
  meanfd
}

## function to register variables
register_var = function(var, 
                        normedVar,
                        normedVariter,
                        cb_atl,
                        len=512,
                        weighted=TRUE,
                        offset=0.0001){
  ## scale data if needed
  if(sum(cb_atl[,var]<0) > 0){
    #cb_atl[,var] = cb_atl[,var] + -min(cb_atl[,var])
    cb_atl[cb_atl[,var]<0,var] = 0
  }
  
  ## ---- data setup
  x = split(cb_atl,factor(cb_atl$SlideID)) ## split data by slide
  rang = range(cb_atl[which(cb_atl[,var]!=0),var]) ## get range of values != 0
  if(rang[1] < 0){rang[1] = 0}
  argvals = seq(rang[1],rang[2],len=len) ## get evenly spaced values for computation later
  
  ## ---- density of var
  densY = sapply(1:length(x), function(i){ ## calculate density for each nonzero value across slides
    #alt_density(x,i,len,rang,var)
    density(x[[i]][,var][which(x[[i]][,var]!=0)],
            from=rang[1],
            to=rang[2],
            n=len,
            na.rm=TRUE)$y
  })
  
  ## ---- setup basis functions
  fdobj_basis = create.bspline.basis(rangeval = rang, norder = 4,nbasis=21) ## create bspline basis with cubic splines (approx hist)
  #fdobj_basis = create.bspline.basis(rangeval = rang, norder = 4,nbasis=25)
  wbasis = create.bspline.basis(rangeval = rang, norder = 2,nbasis=2) ## create bspline basis with linear (transform data)
  Wfd0   <- fd(matrix(0,wbasis$nbasis,1),wbasis) ## setup `fda` object
  WfdPar <- fdPar(Wfd0, Lfdobj = int2Lfd(0),lambda = 0) ## setup roughness for `fda` object
  
  ## ---- initial registration (warp functions)
  fdobj   <- smooth.basis(argvals, densY, fdobj_basis, fdnames = c("x", "samples", "density"))$fd ## estimates the densities using bsplines 
  #regDens = register.fd(yfd=fdobj, WfdParobj = WfdPar,dbglev = 0,crit=1) ## register densities using roughness penalties
  
  ## ---- reverse registration (inverse warp functions)
  if(weighted==TRUE){
    x1 = split(cb_atl[,var],factor(cb_atl$SlideID)) ## split data by slide
    w = sapply(x1, function(y) sum(y>offset))
    y0s = weighted.mean.fd(fdobj, w=w)
  }else{
    y0s = mean.fd(fdobj) ## get mean curve from registration
  }
  y0s$coefs = do.call(cbind, rep(list(y0s$coefs), ncol(fdobj$coefs)))
  ## crit=2: broken in internal fda code
  regDens = register.fd(y0fd = fdobj, yfd=y0s,WfdParobj = WfdPar, dbglev = 0,crit = 1) ## register to get actual warping functions back to real data
  
  ## ---- register raw data
  xp = lapply(1:length(x), ## register each slide using inverse warp functions
              function(ind){ 
                x[[ind]] = data.frame(x[[ind]])
                x[[ind]][,normedVar] = 0;
                x[[ind]][which(x[[ind]][,var]>0),normedVar] = eval.fd(x[[ind]][which(x[[ind]][,var]>0), var], regDens$warpfd[ind]);
                x[[ind]]})
  
  ##cb_atl = data.frame(data.table::rbindlist(xp))
  
  # fittedY = as.matrix(eval.fd(y0s,argvals))
  # reg = list(regDens,fittedY,argvals)
  # saveRDS(reg,
  #         paste0('~/../..//media/disk2/atlas_mxif/combat/handling_zeroes/method12_0423/registered_files/registration_files',normedVariter,'.Rds'))
  
  
  return(data.frame(data.table::rbindlist(xp))) ## combine data into initial form with registered values included
}

## (1) the unregistered curve (blue dashed line)
## (2) the target curve (red dashed line)
## (3) the registered curve (blue solid line).