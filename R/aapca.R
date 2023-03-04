#' @export
pcaest <- function(x = data.table(date, retmat),
									 iscale = c('cov', 'cor'),
									 method = c( 'ML','unbiased'),
									 signmethod=c('ordered','reference'),
									 rcref='WC-3', #reference series name
									 refpol=c(1,-1,-1), #polarity associated
									 rollsum=1, #apply rollsum
									 verbose=T,
									 center=T,
									 pcawin=x[,range(date)],
									 rotwin=c(x[,sort(date)[2]],pcawin[2]), #default to eliminate first period (SNR motive)
									 rotate=T,
									 krot=2:3
									 ){
	# - [ ]  eigen, polarity
  #browser()
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
	dimnames(x3$vectors) <- dimnames(x1$cov) #this is assumed
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
		xrs = x0, #has rollsum applied
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
	plot(cumsum(pcaz(x4))[,1:3],scr=1,col=1:3)
	x4
}

#' @export
rrr2 <- function(#unique rotation from a real number
	tantheta = .1) {
	# - [ ]  #rotate 2D [conflict with aappd::rr2 so renamed from rr2 to rrr2 211126]
	th <- atan(tantheta)
	matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, 2)
}

#' @export
rr3 <- function(tantheta = 1,
								jrot = 3,
								#rotate in the (x1,xjrot) plane
								jbar = 3) {
	# - [ ]  rotate 2 out of jbar
	x1 <- diag(jbar)
	if (1 < jrot) {
		x1[matrix(c(c(1, jrot, 1, jrot),
								c(1, 1, jrot, jrot)), 4, 2)] <- rrr2(tantheta)
	}
	x1
}


#' @export
tanrot <- function(tantheta = rep(0, 3),
									 jbar = 5) {
	# - [ ]  anglevector -> rotation
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
pcah <- function(xest) {
	# - [ ] h holdings (unit variance factor portfolios)
	nbar <- ncol(xest$x[, -'date'])
	kbar <- ncol(xest$eig$vectors)
	if (xest$par$iscale == 'cov') {
		sigmainv <- rep(1, nbar)
	} else if (xest$par$iscale == 'cor') {
		#sigmainv <- 1 / diag(xest$cov$cov)
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
pcab <- function(xest) {
	# - [ ] beta : only when g=1 are these covariance-related, otherwise 'sensitivity' to a scaled factor
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
pcaz <- function(xest = pcaestd,
								 x = xest$x,
								 h = pcah(xest)) {
	# - [ ]  factor timeseries
	zoo(as.matrix(x[, -'date']) %*% h, x[, date])
}


#' @export
pcaobj <- function(tantheta = 0,
									 xest = y2a,
									 j = 2,
									 years=20
									 ) {
	# - [ ]
	stopifnot(1 < j) #see definition of rotation
	dstart <- Sys.Date()-years*365.25
	xest$tantheta[j] <- tantheta
	diff(range(cumsum(pcaz(xest)[dstart<=xest$date, j])))
}

#' @export
pcajscale <- function(xest,
											jscale = c('var', 'gross', 'long', 'short'),
											beta1=F #modified default to FALSE
											) {
	# - [ ]
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
	if(beta1) { #mean(beta)=1 - but really want both b and h to sum to 1 so not satisfactory
		print(paste0('mean pcab=',mean(pcab(x1)[,1])))
		x1$g[1] <- x1$g[1]/mean(pcab(x1)[,1])
		x1$par$beta1 <- T
	} else {
		x1$par$beta1 <- F
	}
	x1
}

#pcarot has been replaced with pcaest(rotate=T)
# pcarot <- function( #this replaces pxmo::f210222b
# 	x, #=dcast(f210310ed[,.(date,cenq=paste0(cen,q),xdot)],date~cenq,value.var='xdot')
# 	krot=2:3,#ncol(x[,-'date']),
# 	rollsum=1, #12, changed 0520
# 	center=T,
# 	signmethod='reference',
# 	rcref=max(names(x[,-'date'])),
# 	window=nrow(x),
# 	method=c('ML','unbiased'),
# 	refpolx=c(1,-1,-1),
# 	startdx=as.Date('1996-12-31'), #added at some point
# 	maxrot=1 #added 220510 but not used
# ) {
# 	# - [ ] class='pcaest' | min-range-rotator: rot(est(pca(f210222ad))) nbar x qbar assets
# 	method <- match.arg(method)
# 	krot <- setdiff(krot,1)
# 	krot <- krot[krot<=ncol(x[,-'date'])]
# 	x1 <- pcaest(
# 		x=x,
# 		iscale='cov',
# 		method=method,
# 		signmethod=signmethod,
# 		rcref=rcref,
# 		refpol=refpolx,
# 		rollsum=rollsum,
# 		center=center,
# 		window=window
# 	)
# 	f1=function(rr=0,jj=2,xm,startd=startdx) { #apply rr with jj
# 		x1 <- copy(xm)
# 		x1$tantheta[jj] <- rr
# 		xx <- diff(range(cumsum(pcaz(xest=x1)[x1[['date']]>=startdx])[,jj,drop=F]))
# 		xx
# 	}
# 	kbar <- x1$par$kbar #x2$par$kbar
# 	initialgrid <- ((-20:20)/21)+.01
# 	for(k in sort(unique(pmin(krot,kbar)))) {
# 		start <- initialgrid[which.min(sapply(initialgrid,f1,xm=x1,jj=k))]
# 		x1$tantheta[k] <- nlm(f=f1,p=start,j=k,xm=x1,startd=startdx)$estimate
# 	}
# 	x1
# }

#' @export
pcarot0 <- function( #design rotation 230223 - normally only to be called by pcaest
	x1, #is a pcaest object
	krot=2:3,#ncol(x[,-'date']),
	irotwin=x1$par$irotwin[1]:x1$par$irotwin[2] #tested and set in pcaest, index into x1
) {
	# - [ ] class='pcaest' | min-range-rotator taking pcaest as input
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
	x1 #is a pcaest object with updated [['tantheta']]
}

