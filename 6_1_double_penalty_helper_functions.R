
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

require(splines)
require(VCA)
require(Matrix)

#-----------------#
#### FUNCTIONS ####
#-----------------#
#-----------------------------------------------------------------------------------.
# (I) B-spline basis and Matrix D_p associated to the difference operator  ---------
#-----------------------------------------------------------------------------------.

pspline.comp<- function (x, xl, xr, ndx, bdeg, pord) {
  dx <- (xr - xl)/ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- spline.des(knots = knots, x = x, ord = bdeg + 1, outer.ok = TRUE)$design
  m = ncol(B)
  n = nrow(B)
  D = diff(diag(m), differences = pord)
  list( B = B, D = D, knots = knots)
  
}

#---------------------------------------------------------------.
# (II) B-spline basis (or its derivative) at new values ---------
#---------------------------------------------------------------.

predict.bbase <- function(knots,bdeg,x, deriv = 0) {
  B <- spline.des(knots = knots, x = x, derivs = deriv,
                  ord =  bdeg + 1, outer.ok = TRUE)$design
  B
} 

#---------------------------------------------------------------.
# (III) Lambda Matrices  ---------
#---------------------------------------------------------------.

Lambda.matrices<-function (g) {
  length.eq <- all(sapply(g, function(x) {
    diff(range(unlist(lapply(ifelse(is.list(x), x, list(x)), 
                             dim)))) < .Machine$double.eps^0.5
  }))
  if (length.eq) {
    l <- length(g)
    if (l == 1) {
      if (length(g[[1]]) == 1) {
        res <- g[[1]]
      }
      else {
        res <- do.call("c", lapply(g, function(x) x))
      }
    }
    else {
      dim <- sapply(g, function(x) {
        if (is.list(x)) 
          unlist(lapply(x, ncol))[1]
        else ncol(x)
      })
      end <- cumsum(dim)
      init <- end - dim + 1
      res <- do.call("c", lapply(1:length(g), function(x, 
                                                       g, init, end, dim) {
        temp <- g[[x]]
        if (is.list(temp)) {
          lapply(temp, function(y, x, dim) {
            aux <- diag(rep(0, sum(dim)))
            aux[init[x]:end[x], init[x]:end[x]] <- y
            aux
          }, x = x, dim = dim)
        }
        else {
          aux <- diag(rep(0, sum(dim)))
          aux[init[x]:end[x], init[x]:end[x]] <- temp
          list(aux)
        }
      }, g = g, init = init, end = end, dim = dim))
    }
  }
  else {
    stop("Error in Lambda.matrices")
  }
  res <- lapply(res, function(x) Matrix(x))
  res
}

#-------------------------------------------------------------------------------------------.
# (IV) Double penalty P-splines model using a modified "SOP" for smoothing parameters ----
#-------------------------------------------------------------------------------------------.
SubjectSpecificCurves.DP <- function(x, y, case ,xgrid = NULL, 
                                     smooth.pop = list(nseg = 10, bdeg = 3, pord1 = 2, pord2 = 0),
                                     smooth.ind = list(nseg = 5, bdeg = 3, pord1 = 2, pord2 = 0),
                                     std.err = FALSE, maxit = 200, thr = 1e-3,
                                     trace = TRUE)
{
  #
  # x     : Covariate
  # xgrid : grid for prediction 
  # y     : Observed data 
  # case  : Id indiv
  # smooth.pop : specifies the components of the P-spline for the population
  # smooth.ind : specifies the components of the P-spline for the individuals
  #
  
  if (trace) start <- proc.time()[3]
  
  #-------------- (IV--a): Create the dataframe with the dataset information --------------------------------
  df <- data.frame(x = x, y = y, case = case)
  df <- df[order( df$case, df$x),] # Order the data (subject, time)
  ids <- unique(df$case) # Number of individuals
  n.ind <- length(unique(case))
  if(!is.null(xgrid)) length.pred <- length(xgrid)
  
  #-------------- (IV--b): Define the prediction grid xp -------------
  if(is.null(xgrid))  {
    xp <- seq(min(df$x), max(df$x),length = 100)#, length =length(x))
  }else {
    xp <- xgrid
  }
  
  #-------------- (IV--b.1): Design matrices for the population curve with the first penalty order (prediction grid xp) ----------
  MM.pop.1 <- pspline.comp(xp, min(xp), max(xp), smooth.pop$nseg, smooth.pop$bdeg, smooth.pop$pord1)
  B.pop.1  <- MM.pop.1$B
  
  # Design matrices for the population derivative curve (prediction grid xp) 
  B.pop.d.1 <- predict.bbase(knots = MM.pop.1$knots, bdeg = smooth.pop$bdeg,x = xp, deriv = 1)
  
  #Matrices for the population curve with the second penalty order (prediction grid xp)
  if(smooth.pop$pord2 > 0) { 
    MM.pop.2  <- pspline.comp(xp, min(xp), max(xp), smooth.pop$nseg, smooth.pop$bdeg, smooth.pop$pord2)
    n.pen.pop <- 2  # Number of penalizations
  } else { 
    MM.pop.2  <- 0
    n.pen.pop <- 1 # Number of penalizations
  }
  
  #-------------- (IV--b.2): Design matrices for the individual curves  with the first penalty order (prediction grid xp)----------
  MM.ind.1 <- pspline.comp(xp, min(xp), max(xp), smooth.ind$nseg, smooth.ind$bdeg, smooth.ind$pord1)
  B.ind.1 <- MM.ind.1$B
  
  # Design matrices for the individual derivative curves (prediction grid xp)
  B.ind.d.1 <- predict.bbase(knots = MM.ind.1$knots, bdeg = smooth.ind$bdeg,x = xp, deriv = 1)
  
  # Matrices for the individual curve with the second penalty order (prediction grid xp)
  if(smooth.ind$pord2 > 0) {
    MM.ind.2  <- pspline.comp(xp, min(xp), max(xp), smooth.ind$nseg, smooth.ind$bdeg, smooth.ind$pord2)
    n.pen.ind <- 3 # Number of penalizations 
  } else  { 
    MM.ind.2  <- 0
    n.pen.ind <- 2 # Number of penalizations (Ridge penalty by default)
  }
  
  #-------------- (IV--c): Construct B-spline matrices for the observed data ------------------------------
  pop.m <- list() 
  ind.m <- list()
  pop.m.d <- list()
  ind.m.d <- list()
  
  count_x_id <- vector(mode = "numeric", length = 0)
  jj <- 0
  
  for(j in 1:n.ind) {
    # Get the observed values of the covariate "x"
    xtmp <- df[df$case == ids[jj + 1], "x"]
    # Count the number of observation by ind
    count_x_id[jj+1] <- length(xtmp) 
    
    # Population B-spline bases 
    Btmp.p <- predict.bbase(knots = MM.pop.1$knots, x = xtmp, bdeg = smooth.pop$bdeg)
    pop.m[[j]]  <- Btmp.p
    
    # Population derivative B-splone bases 
    Btmp.p.d <- predict.bbase(knots = MM.pop.1$knots, x = xtmp, bdeg = smooth.pop$bdeg, deriv = 1)
    pop.m.d[[j]] <- Btmp.p.d
    
    
    # Individual B-spline bases (deviation)
    Btmp <- predict.bbase(knots = MM.ind.1$knots, x = xtmp, bdeg = smooth.ind$bdeg)
    ind.m[[jj + 1]] <- Btmp
    
    # Individual derivative B-spline bases (deviation)
    Btmp.d <-  predict.bbase(knots = MM.ind.1$knots, bdeg = smooth.ind$bdeg, x = xtmp, deriv = 1)
    ind.m.d[[jj + 1]] <- Btmp.d
    
    jj <- jj + 1
  }
  
  
  # Bases for the population and individual component 
  B.pop.m <- Matrix:::Matrix(do.call(rbind,pop.m)) #B(x)
  B.ind.m <- Matrix:::bdiag(ind.m) #B_i(x)
  
  # B-spline matrix with the information of all individuals
  B_g <- as.matrix(cbind(B.pop.m, B.ind.m))
  
  # Derivative B-spline base for the population and individual component
  B.pop.m.d <-  Matrix:::Matrix(do.call(rbind,pop.m.d)) #B'(x)
  B.ind.m.d <- Matrix:::bdiag(ind.m.d) #B_i'(x)
  
  # Derivative B-spline matrix with the information of all individuals
  B_g.d <- as.matrix(cbind(B.pop.m.d, B.ind.m.d))
  
  ind.ind.trt <- rep(1, n.ind)
  C.trt = matrix(rep(1, n.ind), ncol = 1)
  
  #-------------- (IV--d): Number of parameters  ----------------------------------------------
  np <- c(ncol(B.pop.1), ncol(B.ind.1))
  
  #-------------- (IV--e): Construct precision matrix:-----------------------------------------
  
  g <- list()
  
  # Smooth main effect: 
  if(is.numeric(MM.pop.2)){
    g[[1]] <- t(MM.pop.1$D)%*%MM.pop.1$D
  } else{
    g[[1]] <- list(t(MM.pop.1$D)%*%MM.pop.1$D, t(MM.pop.2$D)%*%MM.pop.2$D)#diag(MM.pop$d)
  }
  
  # List with the overlaping penalties for the individual curves;
  # Difference penalty; and
  # Ridge penalty
  
  if(is.numeric(MM.ind.2)){
    g[[2]] <- list(kronecker(diag(n.ind), t(MM.ind.1$D)%*%MM.ind.1$D), kronecker(diag(n.ind), diag(ncol(B.ind.1))))
  } else {
    g[[2]] <- list(kronecker(diag(n.ind), t(MM.ind.1$D)%*%MM.ind.1$D), kronecker(diag(n.ind), diag(ncol(B.ind.1))),kronecker(diag(n.ind), t(MM.ind.2$D)%*%MM.ind.2$D))
  }
  
  # Construct the components of the precision matrix (as needed by the algorithm)
  g <- Lambda.matrices(g)
  
  if(is.numeric(MM.pop.2) & is.numeric(MM.ind.2)){
    names(g)  <- c("pop.1","ind.diff.1","ind.ridge") # There is not a second penalty
  } else if(!is.numeric(MM.pop.2) & is.numeric(MM.ind.2)){
    names(g)  <- c("pop.1","pop.2","ind.diff.1","ind.ridge") # There is not a second penalty in the individual term
  } else if(is.numeric(MM.pop.2) & !is.numeric(MM.ind.2)){
    names(g)  <- c("pop.1","ind.diff.1","ind.ridge","ind.diff.2") # There is not a second penalty in the population term
  } else{  
    names(g)  <- c("pop.1","pop.2","ind.diff.1","ind.ridge","ind.diff.2")
  }
  
  #-------------- (IV--f):  Perform the optimization, extract the estimates, compute sd, etc ------------------------------
  
  # Initialize the parameters
  la <- c(1, rep(1, l = length(g)))
  devold <- 1e10
  dla1 <- 0
  
  nobs <- length(y)
  Hinv.warning <- vector()
  
  if(trace) {
    cat("Effective dimensions\n")
    cat("-------------------------\n")
    cat(sprintf("%1$3s %2$12s","It.","dla"), sep = "")
    cat(sprintf("%14s", names(g)), sep = "")
    cat("\n")
  }
  
  # Start  
  z <- df$y 
  # B'B
  V  <- Matrix:::Matrix(crossprod(B_g))
  # B'y
  u <- Matrix:::Matrix(t(B_g)%*%z)
  
  # start the algorithm to optimization (Lambdas and coefficients)
  for (it in 1:maxit) {
    
    # Compute the precision marix G^-1
    Ginv <- 0
    for (i in 1:length(g)) {
      Ginv <- Ginv + (1/la[i+1])*g[[i]]
    }
    
    # Construct G in a cleverer way (the only problem is for the population)
    n.param.pop <- ncol(B.pop.1)
    aux <- 0
    for(i in 1:n.pen.pop) {
      aux <- aux + (1/la[i+1])*g[[i]]
    }
    
    aux <- aux[1:n.param.pop, 1:n.param.pop]
    G1 <-  Matrix:::Matrix(MPinv(aux))
    
    aux <- 0
    for(i in (n.pen.pop + 1):(n.pen.pop + n.pen.ind)) {
      aux <- aux + (1/la[i+1])*g[[i]]
    }
    
    aux <-  Matrix:::Matrix(aux[(n.param.pop+1):ncol(Ginv), (n.param.pop+1):ncol(Ginv)])
    G2 <- solve(aux)      
    
    # Compute G matrix
    G <- Matrix::bdiag(G1, G2)
    
    # Compute Hinverse matrix
    H <- (1/la[1])*V + Ginv
    
    Hinv <- try( Matrix:::Matrix(solve(H)))
    if(class(Hinv)[1] == "try-error") {
      Hinv <- Matrix:::Matrix(MPinv(H))
    }
    
    Hinv.warning[it] <- ifelse(class(Hinv)[1] == "try-error", 1, 0) 
    
    # Estimate the coefficients
    b <- (1/la[1])*Hinv%*%u
    
    # Estimate the smoothing parameters with SOP algorithm 
    aux <- G - Hinv
    updates <- lapply(1:length(g), function(i, g, la, b, aux) {
      g.inv.d <- (1/la[i+1])*g[[i]]
      ed <- sum(colSums(t(aux)*g.inv.d))
      tau <- as.vector((t(b)%*%g[[i]]%*%b)/ed)
      tau <-  ifelse(tau <= 1e-20, 1e-20,tau)
      res <- list(ed = ed, tau = tau)
      
      res
    },	 g = g, la = la, b = b, aux = aux)
    
    
    ed <- as.list(unlist(updates)[seq(1, 2*length(g), by = 2)]) # Effective dimension values
    tau <- as.list(unlist(updates)[seq(2, 2*length(g), by = 2)]) # Variance component values
    
    ssr <- t(z)%*%z - t(b) %*% (2*u - V %*% b)
    ed.tot <- sum(colSums(t(Hinv)*V))*((1/la[1])) 
    phi <- as.numeric((ssr / (length(z) - ed.tot)))
    
    # New variance components and convergence check
    lanew <- c(phi, unlist(tau))
    
    dla <- sqrt(sum((lanew - la)^2)/sum(lanew^2))
    
    if(trace) {
      cat(sprintf("%1$3d %2$12.6f", it, dla), sep = "")
      cat(sprintf("%14.3f", unlist(ed)), sep = "")
      cat('\n')
    }
    if (dla < thr || abs(dla1-dla) < 0.0000001) break #####
    if (sum(ed<0)!=0) break
    
    la <- lanew
    dla1 <- dla
    
  }
  ##  Curve fitted Values
  eta <- B_g%*%b  # fit individual curve in obs data
  
  ##  Pop Curve fitted Values
  eta.2 <- B.pop.m%*%b[1:ncol(B.pop.m)] # fit population curve in obs data
  
  # End 
  
  if(sum(unlist(ed)<0)!=0){
    res <- print(" eds are Negative")
    attr(res, 'ed negativo' ) <- df
    class(res) <- "try-error"
    res
  } else {
    if(is.numeric(MM.pop.2) & is.numeric(MM.ind.2)){
      names(g) <- names(ed) <- names(tau) <- c("pop.1","ind.diff.1","ind.ridge")
    } else if(!is.numeric(MM.pop.2) & is.numeric(MM.ind.2)){
      names(g) <- names(ed) <- names(tau) <- c("pop.1","pop.2","ind.diff.1","ind.ridge")
    }else if(is.numeric(MM.pop.2) & !is.numeric(MM.ind.2)){
      names(g) <- names(ed) <- names(tau) <- c("pop.1","ind.diff.1","ind.ridge","ind.diff.2")
    }else {  
      names(g) <- names(ed) <- names(tau) <- c("pop.1","pop.2","ind.diff.1","ind.ridge","ind.diff.2")
    }
    
    mm.ind.trt = matrix(1, ncol = 1, nrow = 1)
    
    
    #-------------- (IV--g): Return values -------------------------------------------------------
    coeff <- b
    res <- list()
    
    res$data <- df
    res$smooth.pop <- smooth.pop
    res$smooth.ind <- smooth.ind
    res$MM <- list(pop.1 = MM.pop.1,pop.2 = MM.pop.2, ind.1 = MM.ind.1, ind.2 = MM.ind.2)
    res$ed <- unlist(ed)
    res$tot_ed <- sum(unlist(ed))
    res$tau <- unlist(tau)
    res$sigma2 <- phi
    res$coeff <- coeff
    res$grid.x <- xp
    res$Hinv.warning<-Hinv.warning
    res$dla <-dla 
    res$convergence <- ifelse(it < maxit, TRUE, FALSE)
    res$Lambda.matrices <- g
    res$dim <- np
    
    # Predicted curve in the grid xp
    
    # Population curve
    dm.pop <- mm.ind.trt%x%B.pop.1
    eta_pop <- matrix(dm.pop%*%coeff[1:ncol(dm.pop)], ncol = 1) # calculated in a grid
    res$fitted.pop.grid <- data.frame("xgrid"=xp,"fitted"=eta_pop)
    
    # Individual curves
    etap <- cbind(C.trt%x%B.pop.1, diag(n.ind)%x%B.ind.1)%*%coeff #calculated in a grid
    aux <- matrix(etap, ncol = n.ind)
    res$fitted.ind.grid <- data.frame("xgrid"=rep(xp, length(unique(case))), "case"=rep(unique(case), each=length.pred),
                                      "fitted"=c(lapply(split(aux, rep(ind.ind.trt, each = nrow(aux))), 
                                                        function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))[[1]]))
    
    # Individual deviation from population
    eta_pop_tot <- c(eta_pop[,1])
    aux <- matrix(etap - rep(eta_pop_tot,n.ind), ncol = n.ind)
    res$fitted.dev.ind.grid <- data.frame("xgrid" = rep(xp, length(unique(case))), "case" = rep(unique(case), each=length.pred),
                                          "fitted" = c(lapply(split(aux, rep(ind.ind.trt, each = nrow(aux))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))[[1]]))
    
    # Derivative  pop grid
    mm.p.d.pred<-Matrix(cbind(C.trt%x%B.pop.d.1, diag(n.ind)%x%matrix(0, ncol = ncol(B.ind.d.1), nrow = nrow(B.ind.d.1)))) #with xp
    
    fit_deriv.pop.pred <- as.matrix(mm.p.d.pred%*%coeff)
    res$deriv.pop.grid<-data.frame("xgrid" = xp, "fitted.d" = fit_deriv.pop.pred[1:length(xp),]) #length(xp) because is the grid measured in all ind.
    
    # Derivative  indiv grid
    mm_deriv.pred <- Matrix(cbind(C.trt%x%B.pop.d.1, diag(n.ind)%x%B.ind.d.1))
    fit_deriv.id.pred <- as.matrix(mm_deriv.pred%*%coeff)
    res$deriv.ind.grid <- data.frame("xgrid" = rep(xp,length(unique(case))), "fitted.d" = fit_deriv.id.pred, "case" = rep(unique(case), each = nrow(B.ind.d.1)))
    
    # Curve fitted with observations
    ## pop
    res$fitted.pop <- data.frame(x = df$x ,fitted = as.vector(eta.2))  # calculated in the obs
    ## Idiv
    res$fitted.ind <- data.frame(x = df$x, case = df$case, fitted = c(eta[,1]))# unlist(split(mu, rep(unique(case),times=count_x_id))))
    # Individual deviation from population
    res$fitted.dev.ind <- data.frame(x = df$x, case = df$case,res$fitted.ind$fitted-res$fitted.pop$fitted)
    row.names(res$fitted.ind)<-1:nrow(df)
    
    
    # Derivative  pop obs.
    mm.p.d<-cbind(B.pop.m.d,Matrix(0, ncol = ncol(B.ind.m.d), nrow = nrow(B.ind.m.d)))
    
    fit_deriv.pop <- as.matrix(mm.p.d%*%coeff) 
    res$deriv.pop<-data.frame("x"=df$x,"fitted.d"=fit_deriv.pop)
    
    # Derivative  indiv obs.
    fit_deriv.id <- as.matrix(B_g.d%*%coeff)
    res$deriv.ind<-data.frame("x"=df$x,"fitted.d"=fit_deriv.id,"case"=df$case)
    
    
    # Calculate standard error
    if(std.err) {
      #  Individual
      ## se ind. pred 
      mm.pred <- Matrix(cbind(C.trt%x%B.pop.1, diag(n.ind)%x%B.ind.1)) #with xp
      se_fit.pred <- sqrt(diag(mm.pred%*%Hinv%*%t(mm.pred)))
      res$std_err.ind.grid <- data.frame("xgrid"=rep(xp,length(unique(case))),"se_fit"=se_fit.pred,"case"=rep(unique(case),each=nrow(B.ind.1)))
      
      ## se ind. data obs 
      se_fit <- sqrt(diag(B_g%*%Hinv%*%t(B_g))) #Grid
      res$std_err.ind <- data.frame("x"=df$x,"se_fit"=se_fit,"case"=df$case)
      
      # se deriv ind pred
      se_fit_deriv.pred <- sqrt(diag(mm_deriv.pred%*%Hinv%*%t(mm_deriv.pred)))
      res$deriv.ind.grid$se_fit <- se_fit_deriv.pred
      
      # se deriv ind obs.
      se_fit_deriv <- sqrt(diag(B_g.d%*%Hinv%*%t(B_g.d)))
      res$deriv.ind$se_fit<-se_fit_deriv
      
      # Population  
      ## se pop. grid
      mm.p.pred <- Matrix(cbind(C.trt%x%B.pop.1, diag(n.ind)%x%Matrix(0, ncol=ncol(B.ind.1), nrow = nrow(B.ind.1)))) #with xp
      se_fit.p.pred <- sqrt(diag(mm.p.pred%*%Hinv%*%t(mm.p.pred)))[1:length(xp)]
      res$std_err.pop.grid <- data.frame("xgrid"=xp,"se_fit"=se_fit.p.pred)
      
      ## se pop. data obs.
      mm.p<-cbind(B.pop.m,Matrix(0, ncol=ncol(B.ind.m), nrow = nrow(B.ind.m)))
      se_fit.p <- sqrt(diag(mm.p%*%Hinv%*%t(mm.p)))
      res$std_err.pop <- data.frame("x"=x,"se_fit"=se_fit.p)   
      
      # se deriv pop grid
      se_fit.p_deriv.pred <- sqrt(diag(mm.p.d.pred%*%Hinv%*%t(mm.p.d.pred)))
      res$deriv.pop.grid$se_fit = se_fit.p_deriv.pred[1:length(xp)]
      
      # se deriv pop obs.
      se_fit.p_deriv <- sqrt(diag(mm.p.d%*%Hinv%*%t(mm.p.d)))
      res$deriv.pop$se_fit = se_fit.p_deriv
      
    }
    
    if (trace) end <- proc.time()[3]
    
    #-------------- (IV--h): comp.time -------------------------------------------------------
    
    res$comp.time <- data.frame(time = end - start, it = it,
                                  nseg = unlist(smooth.ind)[1], bdeg = unlist(smooth.ind)[2],
                                  pord1 = unlist(smooth.ind)[3],
                                  pord2 = unlist(smooth.ind)[4],
                                  Model = ifelse(smooth.ind$pord2 == "0", "DC", "DP"))
    
    if (trace) {
      end <- proc.time()[3]
      cat("All process", end - start, "seconds\n")
    }
    class(res) <- "SubjectSpecificCurves.SOP"
    res
  }
}