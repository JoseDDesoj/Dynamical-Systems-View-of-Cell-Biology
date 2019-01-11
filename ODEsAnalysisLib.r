################################################################################################################################
library(deSolve)
library(rootSolve)
################################################################################################################################
tarrows <- function(out,ds=0.1,...)
{
 # select point at predefined distance from point (p1,p2)
 p   <- unlist(out[nrow(out),2:3])
 dd  <- (out[,2]-p[1])^2+ (out[,3]-p[2])^2
 dd2 <-  c(dd[-1],dd[1])
 i1<-which(dd<ds&dd2>ds | dd>ds&dd2<ds)

 p   <- unlist(out[1,2:3])
 dd  <- (out[,2]-p[1])^2+ (out[,3]-p[2])^2
 dd2 <-  c(dd[-1],dd[1])
 i2<-which(dd<ds&dd2>ds | dd>ds&dd2<ds)[1]
 ii <- c(i1,i2)
# ii <- which(dd>ds)
# iseq <- seq(1,length (ii),15)
 for (i in ii ) arrows(out[i,2],out[i,3],out[i+1,2],out[i+1,3],length=0.1,lwd=1,...)

}
################################################################################################################################
trajectory <- function(Func, X, parameters, ds=0.1, Col=1) {
times <- seq(0, 100, by = 0.01)
out    <- as.data.frame(ode(X, times, Func, parms = parameters))

matplot(out$X1,out$X2,type="l", lwd=2, add=TRUE, col=Col)
tarrows(out,ds)
}
################################################################################################################################
trajectoryAttractor <- function(Func, X, parameters, Attractors, ds=0.1, Cols) {
	times <- seq(0, 50, by = 0.05)
	out    <- as.data.frame(ode(X, times, Func, parms = parameters))
	#c(X1=X[1], X2=X[2])
	#print(c(X[1], X[2]))
	#print(out)

	Dists <- t(sapply(1:nrow(Attractors), function(i) {
		dist(rbind(out[nrow(out),2:3],Attractors[i,]))
	}))
	#print(Dists)
	#print(which(Dists==min(Dists)))
	#print(Attractors[which(Dists==min(Dists)),])
	Col=Cols[which(Dists==min(Dists))]

	matplot(out$X1,out$X2,type="l", lwd=1, add=TRUE, col=Col)
	tarrows(out,ds,col=Col)
}
################################################################################################################################
trajectoryAttractor2 <- function(Func, X, parameters, Attractors, ds=0.1, Cols) {
	times <- seq(0, 50, by = 0.05)
	out    <- as.data.frame(ode(X, times, Func, parms = parameters))
	#c(X1=X[1], X2=X[2])
	#print(c(X[1], X[2]))
	#print(out)

	if(abs(out[nrow(out),2]-out[nrow(out),3])<=0.1) Col <- Cols[2] #Simetrico
	if(out[nrow(out),2]-out[nrow(out),3]<(-0.1)) Col <- Cols[1]   #Simetrico
	if(out[nrow(out),2]-out[nrow(out),3]>0.1) Col <- Cols[3]      #Simetrico
	#if(out[nrow(out),2]>0.8 & out[nrow(out),3]>0.65) Col <- Cols[2] # Verde
	#if(out[nrow(out),2]<0.8 & out[nrow(out),3]>0.65) Col <- Cols[1]   # Azul
	#if(out[nrow(out),2]>0.8 & out[nrow(out),3]<0.65) Col <- Cols[3]    # Rojo
	matplot(out$X1,out$X2,type="l", lwd=1, add=TRUE, col=Col)
	tarrows(out,ds,col=Col)
}
################################################################################################################################
RandomInitial <- function(Nvariables, Nconditions, Rango) {
Nums <- matrix(runif(Nconditions*Nvariables, Rango[1], Rango[2]), Nconditions, Nvariables)
colnames(Nums) <- c("X1", "X2")
return(Nums)
}
################################################################################################################################
RandomInitialRango <- function(Nconditions, RangoV1, RangoV2) {
	V1 <- runif(Nconditions, RangoV1[1],RangoV1[2])
	V2 <- runif(Nconditions, RangoV2[1], RangoV2[2])
	Inis <- as.matrix(cbind(V1, V2))
	colnames(Inis) <- c("X1", "X2")
return(Inis)
}
################################################################################################################################
GraphEqPointsNewton <- function(Function, parameters, Inits) {
	t(sapply(1:nrow(Inits), function(i) {
		Equil <- stode(y = c(X1=Inits[i,1], X2=Inits[i,2]), fun = Function, parms = parameters, pos = TRUE)$y
		points(Equil[1], Equil[2], col=2, pch=20)
	}))
}
################################################################################################################################
GraphEqPointsSimulation <- function(Function, parameters, Inits) {
	t(sapply(1:nrow(Inits), function(i) {
		Equil <- runsteady(y = c(X1=Inits[i,1], X2=Inits[i,2]), fun = Function, parms = parameters, times = c(0, 1e5))$y
		points(Equil[1], Equil[2], col=2, pch=20)
	}))
}
################################################################################################################################
GraphEqPointsStability <- function(Function, parameters, Inits) {
	t(sapply(1:nrow(Inits), function(i) {
		Equil <- stode(y = Inits[i,], fun = Function, parms = parameters, pos = TRUE)$y
		Jacob <- jacobian.full(y=Equil,func=Function, parms=parameters)
		EigenVal <- eigen(Jacob)$values
		# white:unstable node, black:stable node, grey:saddle
		if (sign(EigenVal[1])>0 & sign(EigenVal[2])>=0) col <- "white"
		if (sign(EigenVal[1])<0 & sign(EigenVal[2])<=0) col <- "black"
		if (sign(EigenVal[1])* sign(EigenVal[2])   <0 ) col <- "grey"
		points(Equil[1], Equil[2],pch=21,cex=2.0, bg=col,col="black")
	}))
}
################################################################################################################################
EqPointsEigenv <- function(Function, parameters, Inits) {
	Equilib <- t(sapply(1:nrow(Inits), function(i) {
		Equil <- stode(y = Inits[i,], fun = Function, parms = parameters, pos = TRUE)$y
	}))
		#print(Equilib)
		EqPoints <- as.matrix(Equilib[which(duplicated(round(Equilib, 2),MARGIN=1)==FALSE), ])
		#print(EqPoints)
		#print(EqPoints)
		EigenVals <- t(sapply(1:nrow(EqPoints), function(i) {
		Jacob <- jacobian.full(y=c(EqPoints[i,]) ,func=Function, parms=parameters)
		EigenVal <- eigen(Jacob)$values
		}))

		return(cbind(EqPoints, EigenVals))
}
################################################################################################################################
EqPoints <- function(Function, parameters, Inits) {
	Equilib <- t(sapply(1:nrow(Inits), function(i) {
		Equil <- stode(y = Inits[i,], fun = Function, parms = parameters, pos = TRUE)$y
	}))
		#print(Equilib)
		if(length(which(duplicated(round(Equilib, 2),MARGIN=1)==FALSE))==1) EqPoints <- t(as.matrix(Equilib[which(duplicated(round(Equilib, 2),MARGIN=1)==FALSE),]))
		else EqPoints <- as.matrix(Equilib[which(duplicated(round(Equilib, 2),MARGIN=1)==FALSE), ])
		#print(EqPoints)
		#print(EqPoints)

		return(cbind(EqPoints))
}
################################################################################################################################
RunAttractors <- function(Function, parameters, Init) {
  Attractors <- t(sapply(1:nrow(Init), function(i) runsteady(y = Init[i,], fun = Function, parms = parameters, times = c(0, 1e5))$y))
  if(length(which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE))==1) Attractors <- t(as.matrix(Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]))
  if(length(which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE))>1) Attractors <- Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]

  return(Attractors)
}
################################################################################################################################
# Input - Output from EqPointsEigenv
# Output - Equilibrium points of type = Type ("Saddle", "Stable", "Unstable")
PointStabilityType <- function(EquStabMatrix, Type) {
if(Type=="Saddle") return(EquStabMatrix[which(t(sapply(1:nrow(EqPointsEig), function(i) sum(EqPointsEig[i,3:4]<0)==1))), 1:2]) 	 # Saddle -,+
if(Type=="Stable") return(EquStabMatrix[which(t(sapply(1:nrow(EqPointsEig), function(i) sum(EqPointsEig[i,3:4]<0)==2))), 1:2]) 	 # Stable -,-
if(Type=="Unstable") return(EquStabMatrix[which(t(sapply(1:nrow(EqPointsEig), function(i) sum(EqPointsEig[i,3:4]<0)==0))), 1:2]) # Instable +,+
}
################################################################################################################################
GraphBifurcationAttractStab <- function(Function, parameters, Inits, ParVal) {
	t(sapply(1:nrow(Inits), function(i) {
		Equil <- stode(y = Inits[i,], fun = Function, parms = parameters, pos = TRUE)$y
		Jacob <- jacobian.full(y=Equil,func=Function, parms=parameters)
		EigenVal <- eigen(Jacob)$values
		# white:unstable node, black:stable node, grey:saddle
		if (sign(EigenVal[1])>0 & sign(EigenVal[2])>=0) Col <- "white"
		if (sign(EigenVal[1])<0 & sign(EigenVal[2])<=0) Col <- "black"
		if (sign(EigenVal[1])* sign(EigenVal[2])   <0 ) Col <- "grey"
		points(ParVal, Equil[2],pch=20,col=Col)
	}))
}
################################################################################################################################
GraphBifurcationAttractStab2 <- function(Function, parameters, Init, ParVal) {
  Attractors <- t(sapply(1:nrow(Init), function(i) runsteady(y = Init[i,], fun = Function, parms = parameters, times = c(0, 1e5))$y))
  if(length(which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE))==1) Attractors <- t(as.matrix(Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]))
  else Attractors <- Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]
  #print(Attractors)
  for(i in 1:nrow(Attractors)) {
    Jacob <- jacobian.full(y=Attractors[i,],func=Function, parms=parameters)
    EigenVal <- eigen(Jacob)$values
    #print(EigenVal)
    #white:unstable node, black:stable node, grey:saddle
    if (sign(EigenVal[1])>0 & sign(EigenVal[2])>=0) Col <- "white"
    if (sign(EigenVal[1])<0 & sign(EigenVal[2])<=0) Col <- "black"
    if (sign(EigenVal[1])* sign(EigenVal[2])   <0 ) Col <- "grey"
    points(ParVal, Attractors[i,2],pch=20, col=Col)
  }
}
################################################################################################################################
GraphBiAttractStab2 <- function(Function, parameters, Init,...) {
  Attractors <- t(sapply(1:nrow(Init), function(i) runsteady(y = Init[i,], fun = Function, parms = parameters, times = c(0, 1e5))$y))
  if(length(which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE))==1) Attractors <- t(as.matrix(Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]))
  if(length(which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE))>1) Attractors <- Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]
  #else Attractors <- Attractors[which(duplicated(round(Attractors, 2),MARGIN=1)==FALSE),]
  #print(Attractors)
  for(i in 1:nrow(Attractors)) {
    print("################################")
    print(Attractors[i,])
    print("################################")
    Jacob <- jacobian.full(y=Attractors[i,],func=Function, parms=parameters)
    EigenVal <- eigen(Jacob)$values
    #print(EigenVal)
    #white:unstable node, black:stable node, grey:saddle
    if (sign(EigenVal[1])>0 & sign(EigenVal[2])>=0) Col <- "white"
    if (sign(EigenVal[1])<0 & sign(EigenVal[2])<=0) Col <- "black"
    if (sign(EigenVal[1])* sign(EigenVal[2])   <0 ) Col <- "grey"
    points(Attractors[i,1], Attractors[i,2], col=Col,...)
  }
}
################################################################################################################################
phasearrows <- function(fun,xlims,ylims,resol=25, col='black', add=F,parms=NULL,jitter=FALSE,...) {
  if (add==F) {
     plot(1,xlim=xlims, ylim=ylims, type='n',...)#xlab="x",ylab="y");
  }
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=T, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=F, resol, resol);
  npts <- resol*resol;
  if(jitter) {
    xspace <- abs(diff(xlims))/(resol*10);
    yspace <- abs(diff(ylims))/(resol*10);
    x <- x + matrix(runif(npts, -xspace, xspace),resol,resol);
    y <- y + matrix(runif(npts, -yspace, yspace),resol,resol);
  }
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  maxx <- max(abs(z1));
  maxy <- max(abs(z2));
  dt <- min( abs(diff(xlims))/maxx, abs(diff(ylims))/maxy)/resol;
  lens <- sqrt(z1^2 + z2^2);
  lens2 <- lens/max(lens);
  arrows(c(x), c(y), c(x+dt*z1/((lens2)+.1)), c(y+dt*z2/((lens2)+.1)),length=.04, col=col);
}
################################################################################################################################
nullclines <- function(fun,xlims, ylims, resol=250, add=F,parms=NULL) {
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1,levels=c(0), drawlabels=F,add=add, lwd=2, col="blue");
  contour(x[,1],y[1,],z2,levels=c(0), drawlabels=F,add=T, lwd=2, col="forestgreen");
  title(main="Blue=x nullcline, Green=y nullcline",cex=0.35);
}
################################################################################################################################
newton=function(func,x0=NULL,parms=NULL,tol=1e-16,niter=40,inc=1e-6,plotit=TRUE) {
  x=x0; #initial x
  if (is.null(x0)){x = locator(n=1); x=c(x$x,x$y)};
  nx = length(x); # length of state vector
  ######### Newton iteration loop: start
  for(i in 1:niter){
   y = func(0,x,parms)[[1]]
   df = matrix(0,nx,nx); # Compute df
   for(j in 1:nx) {
	#Increment vector for estimating derivative wrt jth coordinate
	v=rep(0,nx);
	v[j] = inc;
        df[,j]=  (func(t,x+v,parms)[[1]] - func(t,x-v,parms)[[1]])/(2*inc)
    }
    if (sum(y^2) < tol){  #check for convergence
        if(plotit){
	     ev=eigen(df)$values; pch1=1+as.numeric(Im(ev[1])!=0); pch2=1+as.numeric(max(Re(ev))<0);
	     pchs=matrix( c(2,17,1,16),2,2,byrow=T);
	     points(x[1],x[2],type="p",pch=pchs[pch1,pch2],cex=1.5)
 	}
	cat("Fixed point (x,y) = ",x,"\n");
	cat("Jacobian Df=","\n"); print(df);cat("Eigenvalues","\n"); print(eigen(df)$values); cat("Eigenvectors","\n"); print(eigen(df)$vectors);
	return(list(x=x,df=df))
    }
    x = x - solve(df,y) # one more step if needed
    cat(i, x, "\n") #print out the next iterate
  }
  ######### Newton iteration loop: end
 cat("Convergence failed");
}
################################################################################################################################
DrawManifolds=function(fun.lsoda,parms,x0=NULL,maxtime=100) {
	xbar=newton(fun.lsoda,x0=x0,parms=parms,plotit=FALSE);
	x=xbar$x; df=xbar$df; V=eigen(df)$vectors; ev=eigen(df)$values;
	if (ev[1]*ev[2] > 0) {
	  cat("Fixed point is not a saddle \n");
	}else{
          i1=which(ev>0); i2=which(ev<0);
	  v1=V[,i1]; v2=V[,i2]; eps=1e-3;
	  out1=lsoda(times=seq(0,maxtime,.1),y=x+eps*v1,func=fun.lsoda,parms=parms); points(out1[,2],out1[,3],type="l",lwd=2,col="red");
	  out2=lsoda(times=seq(0,maxtime,.1),y=x-eps*v1,func=fun.lsoda,parms=parms); points(out2[,2],out2[,3],type="l",lwd=2,col="red");
	  out3=lsoda(times=-seq(0,maxtime,.1),y=x+eps*v2,func=fun.lsoda,parms=parms); points(out3[,2],out3[,3],type="l",lwd=2,col="black");
	  out4=lsoda(times=-seq(0,maxtime,.1),y=x-eps*v2,func=fun.lsoda,parms=parms); points(out4[,2],out4[,3],type="l",lwd=2,col="black");
	  title(sub="Black=stable manifold, Red=unstable manifold");
	}
}
################################################################################################################################
