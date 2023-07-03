#' @export
pcaest <- 
  function(# - [ ]  eigen, polarity
    x = data.table(date, retmat),
    iscale = c('cov', 'cor'),
    method = c( 'ML','unbiased'),
    signmethod=c('ordered','reference'),
    rcref='WC-3', #reference series name
    refpol=c(1,-1,-1), #polarity associated
    rollsum=1, #apply rollsum - not used
    verbose=T,
    center=T,
    pcawin=x[,range(date)],
    rotwin=c(x[,sort(date)[2]],pcawin[2]), #default to eliminate first period (SNR motive)
    rotate=T,
    krot=2:3,
    doplot=F
  ){
    iscale <- match.arg(iscale)
    method <- match.arg(method)
    signmethod <- match.arg(signmethod)
    pcawin <- sort(round(as.Date(pcawin)))
    stopifnot(is.Date(pcawin)&all(pcawin%in%x[,date])&length(pcawin)==2)
    rotwin <- sort(round(as.Date(rotwin)))
    stopifnot(is.Date(rotwin)&all(rotwin%in%x[,date])&length(rotwin)==2)
    ipcawin <- setNames(match(pcawin,x[,date]),pcawin)
    irotwin <- setNames(match(rotwin,x[,date]),rotwin)
    nbar <- ncol(x) - 1
    x0 <- rollapply(x[,-'date'],width=rollsum,FUN=sum,partial=T,align='right')
    if(rollsum!=1&verbose) {print(paste0('rollsum=',rollsum,' in pcaest'))}
    x1 <- cov.wt(x0[ipcawin[1]:ipcawin[2],],
                 method = method,
                 center = center,
                 cor = T)
    x2 <- x1[[iscale]]
    x3 <- eigen(x = x2)
    dimnames(x3$vectors) <- dimnames(x1$cov) 
    thresh <- sqrt(.Machine$double.eps)
    if(any(x3$values<thresh)) {
      kbar <- max(which(x3$values>thresh))
      x3$values <- x3$values[1:kbar]
      x3$vectors <- x3$vectors[,1:kbar,drop=F]
    } else {
      kbar <- ncol(x2)
    }
    if(signmethod=='ordered') {
      print('ordered method')
      signfinder <-
        function(evec, pola = c(1, -1, -1)[1:min(3,length(evec))]) {
          #sign by regression on 1st 3 even zero phase cos(x.centred); pola is tgt for last
          n <- length(evec) - 1
          x <-
            data.table(
              y = evec,
              f0 = rep(1, n + 1),
              f1 = -cos(pi * (0:n) / n),
              f3 = cos(-2 * pi * (0:n) / n)
            )
          x1 <- summary(lm(y ~ . - 1, x))$coefficients
          if(any(is.na(x1[,2]))){
            x2 <- sign(sum(x1[, 1] * pola))
          } else {
            x2 <- sign(sum(x1[, 3] * pola))#sum of t-stats (scale-invariant)
          }
          x2
        }
      x4 <- unlist(lapply(data.table(x3$vectors), signfinder))
      if(any(is.na(x4))) {x4 <- rep(1,length(x4))}
    } else {
      stopifnot(rcref%in%rownames(x3$vectors))
      iref <- match(rcref,rownames(x3$vectors))
      jref <- seq_along(refpol)
      x4 <- c(refpol*sign(x3$vectors[iref,jref]),rep(1,ncol(x3$vectors)-length(jref)))
    }
    x3$vectors <- sweep(x3$vectors,
                        STAT = x4,
                        MAR = 2,
                        FUN = `/`)
    x4 <- list(
      x = x,
      xrs = x0, #rs 'rollsum applied'
      date = x[, date],
      sigma=sqrt(diag(x1$cov)),
      xx = x2,
      eig = x3,
      tantheta = rep(0, nbar),
      g = rep(1, nbar),
      par = list(
        method = method,
        iscale = iscale,
        xsect = '',
        rollsum=rollsum,
        kbar=kbar,
        ipcawin=ipcawin,
        irotwin=irotwin
      )
    )
    if(rotate) {
      x4 <- pcarot0(x4)
    }
    if(doplot) {plot(cumsum(pcaz(x4))[,1:3],scr=1,col=1:3)}
    x4
  }

#' @export
rrr2 <-
  function( # - [ ]  #rotate 2D from a real number
    tantheta = .1
  ) {
    th <- atan(tantheta)
    matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, 2)
  }

#' @export
rr3 <-
  function( # - [ ]  rotate 2 out of jbar in the (x1,xjrot) plane
    tantheta = 1,
    jrot = 3,
    jbar = 3
  ) {
    x1 <- diag(jbar)
    if (1 < jrot) {
      x1[matrix(c(c(1, jrot, 1, jrot),
                  c(1, 1, jrot, jrot)), 4, 2)] <- rrr2(tantheta)
    }
    x1
  }

#' @export
tanrot <-
  function( # - [ ]  anglevector -> rotation
    tantheta = rep(0, 3),
    jbar = 5
  ) {
    tantheta[1] <- 0 #all rotations are perp. to plane {1,j}
    x1 <- list(diag(jbar))
    for (j in which(tantheta!=0)) {
      x1[[j]] <- rr3(tantheta = tantheta[j],
                     jrot = j,
                     jbar = jbar)
    }
    Reduce(`%*%`, x1)
  }

#' @export
pcah <-
  function( # - [ ] h holdings (unit variance factor portfolios)
    xest
  ) {
    nbar <- ncol(xest$x[, -'date'])
    kbar <- ncol(xest$eig$vectors)
    if (xest$par$iscale == 'cov') {
      sigmainv <- rep(1, nbar)
    } else if (xest$par$iscale == 'cor') {
      sigmainv <- 1 /xest$sigma
    }
    x3 <- diag(sigmainv) %*%
      xest$eig$vectors %*%
      diag(1 / sqrt(xest$eig$values)) %*%
      tanrot(xest$tantheta,jbar=kbar) %*%
      diag(1 / xest$g[1:kbar])
    rownames(x3) <- names(xest$x[, -'date'])
    colnames(x3) <- zeroprepend(1:ncol(x3),3)
    x3
  }

#' @export
pcab <-
  function( # - [ ] beta : only when g=1 are these covariance-related, otherwise 'sensitivity' to a scaled factor
    xest
  ) {
    nbar <- ncol(xest$x[, -'date'])
    kbar <- ncol(xest$eig$vectors)
    if (xest$par$iscale == 'cov') {
      sigma <- rep(1, nbar)
    } else if (xest$par$iscale == 'cor') {
      sigma <- xest$sigma
    }
    x3 <- diag(sigma) %*%
      xest$eig$vectors %*%
      diag(sqrt(xest$eig$values)) %*%
      tanrot(xest$tantheta,jbar=kbar) %*%
      diag(xest$g[1:kbar])
    rownames(x3) <- names(xest$x[, -'date'])
    colnames(x3) <- zeroprepend(1:ncol(x3),3)
    x3
  }

#' @export
pcaz <-
  function( # - [ ]  factor timeseries
    xest = pcaestd,
    x = xest$x,
    h = pcah(xest)
  ) {
    zoo(as.matrix(x[, -'date']) %*% h, x[, date])
  }

#' @export
pcaobj <-
  function( # - [ ] objective function for rotation
    tantheta = 0,
    xest = y2a,
    j = 2,
    years=20
  ) {
    stopifnot(1 < j) #see definition of rotation
    dstart <- Sys.Date()-years*365.25
    xest$tantheta[j] <- tantheta
    diff(range(cumsum(pcaz(xest)[dstart<=xest$date, j])))
  }

#' @export
pcajscale <-
  function( # - [ ] column scalar
    xest,
    jscale = c('var', 'gross', 'long', 'short'),
    beta1=F #modified default to FALSE
  ) {
    jscale <- match.arg(jscale)
    x1 <- copy(xest)
    x1$g <- rep(1, ncol(xest$x) - 1) #unscaled (=varscaled)
    x2 <- pcah(x1)
    if (jscale ==         'gross') {
      x1$g <- apply(abs(x2), 2, sum)*.5
    } else if (jscale ==  'long')  {
      x1$g <- apply(x2 * (x2 > 0), 2, sum)
    } else if (jscale ==  'short') {
      x1$g <- apply(abs(x2) * (x2 < 0), 2, sum)
    }
    x1$par$jscale <- jscale
    if(beta1) { #mean(beta)=1 
      print(paste0('mean pcab=',mean(pcab(x1)[,1])))
      x1$g[1] <- x1$g[1]/mean(pcab(x1)[,1])
      x1$par$beta1 <- T
    } else {
      x1$par$beta1 <- F
    }
    x1
  }

#' @export
pcarot <- 
  function(
    x,...
  ) {
    pcaest(x=x,...)
  }

#' @export
pcarot0 <- 
  function( # - [ ] class='pcaest' | min-range-rotator taking pcaest as input
    x1, #pcaest object
    krot=2:3,
    irotwin=x1$par$irotwin[1]:x1$par$irotwin[2]
  ) {
    krot <- setdiff(krot,1)
    krot <- krot[krot<=ncol(x1$x[,-'date'])]
    f1=function(rr=0,jj=2,xm,irotwin) { #apply rr with jj
      x1 <- copy(xm)
      x1$tantheta[jj] <- rr
      xx <- diff(range(cumsum(pcaz(xest=x1)[irotwin])[,jj,drop=F]))
      xx
    }
    kbar <- x1$par$kbar 
    initialgrid <- ((-20:20)/21)+.01
    kset <- sort(unique(pmin(krot,kbar)))
    for(i in seq_along(kset)) {
      k <- kset[i]
      start <- initialgrid[which.min(sapply(initialgrid,f1,xm=x1,jj=k,irotwin=irotwin))]
      x1$tantheta[k] <- nlm(f=f1,p=start,j=k,xm=x1,irotwin=irotwin)$estimate
    }
    x1 #with updated tantheta
  }


#' @export
pcadrc <-  # - [ ] class='pcaest' | min-range-rotator taking pcaest as input
  function(
    pca,
    date0=as.Date('1994-12-31'),
    kbar=3
) {
  dates <- c(date0,pca$date)
  x2 <- x1 <- NA
  for(i in 1:nrow(pca$x)) {
    x0 <- data.table(
        y=unlist(t(pca$x[i,-1])),
        x=pcab(pca)[,1:kbar]
      )%>%
      setnames(.,c('y',paste0('b',1:kbar)))%>%
      lm(y~.,.)%>%
      summary(.)
    x1[i] <- x0[['adj.r.squared']]
    x2[i] <- x0[['r.squared']]
  }
  x3 <- data.table(
    days=as.numeric(diff(dates)),       #days
    rbarsq=x1,                          #rbarsq
    rsq=x2,
    r=apply(pcaz(pca)[,2:kbar]^2,1,sum),#k=2,3 variance
    start=dates[-length(dates)],        #period start
    end=dates[-1],                       #period end
    pcaz(pca)[,2:3]%>%
      coredata(.)%>%
      data.table(.)%>%
      setnames(.,c('z2','z3'))
  )
  x3
}
