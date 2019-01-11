###############################################################################################
TwoGeneModel = function(t, x, parms) {
	# Since parms has names, it must be called as a list
	with(as.list(c(x, parms)),{
		dx1 = a1*X1^nHill/(Tet^nHill + X1^nHill)   +   b1*Tet^nHill/(Tet^nHill + X2^nHill) - k1*X1
		dx2 = a2*X2^nHill/(Tet^nHill + X2^nHill)   +   b2*Tet^nHill/(Tet^nHill + X1^nHill) - k2*X2
		dX <- c(dx1, dx2)
		list(dX)
	})
}
###############################################################################################
twoGeneCircuit <- function(X1,X2,parms) {
	with(as.list(c(parms)), {
	c(a1*X1^nHill/(Tet^nHill + X1^nHill)   +   b1*Tet^nHill/(Tet^nHill + X2^nHill) - k1*X1,
	  a2*X2^nHill/(Tet^nHill + X2^nHill)   +   b2*Tet^nHill/(Tet^nHill + X1^nHill) - k2*X2)
	})
}
###############################################################################################
# Flower GRN Module ODEs
###############################################################################################
W.AG <- function(EMF1, TFL1, AP2, WUS, LFY, AP1, AG, SEP) {
  return(    LFY*(1-EMF1) * (1-((AP1*AP2*(1-WUS)) * (1-AG*SEP*(1-TFL1)))) + (1-EMF1)*(1-TFL1)*(1-AP2) - (LFY*(1-EMF1)*(1-((AP1*AP2*(1-WUS))*(1-AG*SEP*(1-TFL1))))) * ((1-EMF1)*(1-TFL1)*(1-AP2))     )
}

W.AP1 <- function(AG, TFL1, FT, LFY, PI, AP3) {
  return(  (1-AG)*(1-TFL1* (1-LFY*FT)) )
}

W.FUL <- function(AP1, TFL1) {
  return( (1-AP1) * (1-TFL1) )
}

W.FT <- function(EMF1) {
  return( 1-EMF1 )
}

W.EMF1 <- function(LFY) {
  return( 1-LFY )
}

W.LFY <- function(EMF1, TFL1) {
  return( 1-EMF1*TFL1 )
}

W.AP2 <- function(TFL1) {
  return( 1-TFL1 )
}

W.WUS <- function(WUS, AG, SEP) {
  return(  WUS*(1-AG*SEP)  )
}

W.SEP <- function(LFY) {
  return(  LFY  )
}

W.UFO <- function(UFO) {
  return(  UFO  )
}

W.PI <- function(LFY, AG, AP3, PI, SEP, AP1) {
  return(  (LFY*(AG+AP3-AG*AP3)) + (PI*SEP*AP3*(AG+AP1-AG*AP1)) - (LFY*(AG+AP3-AG*AP3))*(PI*SEP*AP3*(AG+AP1-AG*AP1))   )
}

W.AP3 <- function(AG, LFY, UFO, AP3, PI, SEP, AP1) {
  return(  (LFY*UFO) + (PI*SEP*AP3*(AG+AP1-AG*AP1)) - (LFY*UFO) * (PI*SEP*AP3*(AG+AP1-AG*AP1))   )
}

W.TFL1 <- function(AP1, LFY, EMF1) {
  return(  (1-AP1)*(1-LFY)*EMF1   )
}

FlorModuleODEs = function(t, x, parms) {
  # Since parms has names, it must be called as a list
  with(as.list(c(x, parms)),{
    dFUL = 1/(1 + exp(-2*b*(W.FUL(AP1, TFL1) - wthr))) - kFUL*FUL
    dFT = 1/(1 + exp(-2*b*(W.FT(EMF1) - wthr))) - kFT*FT
    dAP1 = 1/(1 + exp(-2*b*(W.AP1(AG, TFL1, FT, LFY, PI, AP3) - wthr))) - kAP1*AP1
    dEMF1 = 1/(1 + exp(-2*b*(W.EMF1(LFY) - wthr))) - kEMF1*EMF1
    dLFY = 1/(1 + exp(-2*b*( W.LFY(EMF1, TFL1) - wthr))) - kLFY*LFY
    dAP2 = 1/(1 + exp(-2*b*( W.AP2(TFL1) -wthr))) - kAP2*AP2
    dWUS = 1/(1 + exp(-2*b*( W.WUS(WUS, AG, SEP) - wthr))) - kWUS*WUS
    dAG =  1/(1 + exp(-2*b*( W.AG(EMF1, TFL1, AP2, WUS, LFY, AP1, AG, SEP) - wthr))) - kAG*AG
    dPI =  1/(1 + exp(-2*b*( W.PI(LFY, AG, AP3, PI, SEP, AP1) -wthr))) - kPI*PI  
    dSEP =  1/(1 + exp(-2*b*( W.SEP(LFY) - wthr))) - kSEP*SEP
    dAP3 =  1/(1 + exp(-2*b*( W.AP3(AG, LFY, UFO, AP3, PI, SEP, AP1) - wthr))) - kAP3*AP3    
    dUFO =  1/(1 + exp(-2*b*( W.UFO(UFO) - wthr))) - kUFO*UFO    
    dTFL1 =  1/(1 + exp(-2*b*( W.TFL1(AP1, LFY, EMF1) - wthr))) - kTFL1*TFL1    
    dX <- c(dFUL,dFT,dAP1,dEMF1,dLFY,dAP2,dWUS,dAG,dPI,dSEP,dAP3,dUFO,dTFL1)
    list(dX)
  })
}
###############################################################################################


