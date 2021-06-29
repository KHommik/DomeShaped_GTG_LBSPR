#################################################################
# The original GTG LBSPR code is produced by Hordyk et al. 2016 #
# Accessed from https://github.com/AdrianHordyk/GTG_LBSPR      #
#################################################################
#################################################################
# Original GTG LBSPR code was modified by adding the options for# 
# dome-shaped selevtivity (normal & log-normal selectivity).    #
# Modifications to the code were made by Hommik et al.          #
#################################################################

#################################################################
#               Some general instructions:                      #
# The code is modified such that it can fit dome-shaped         #
# selectivity (normal and log-normal).                          #
# User has to pre-specify selectivity parameters derived from   #
# SELECT method (Millar & Fryer, 1999).                         #
# Selectivity parameters are defined in lines: 474-476.         #
# In line 477 minimum landing limit (MLL) can be defined.       #
# If MLL is defined then model assumes that selection           #
# probability for lengths < MLL is 0.                           #
#################################################################

#runmod<-GTGLBSPRSim(StockPars, Fleet, SizeBins, sel)

GTGDomeLBSPRSim <- function(StockPars, FleetPars, SizeBins=NULL,sel)  { 
  #=c("Logistic","Normal","Knife")						  
  sink(stdout(), type="message")  
  # Assign Variables 
  NGTG <- StockPars$NGTG 
  GTGLinfBy <- StockPars$GTGLinfBy 
  if (!exists("GTGLinfBy")) GTGLinfBy <- NA
  if (is.null(GTGLinfBy)) GTGLinfBy <- NA
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  MaxSD <- StockPars$MaxSD 
  MK <- StockPars$MK 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Walpha <- StockPars$Walpha 
  Wbeta <- StockPars$Wbeta 
  FecB <- StockPars$FecB 
  Steepness <- StockPars$Steepness 
  Mpow <- StockPars$Mpow
  R0 <- StockPars$R0 
  
  SL50 <- FleetPars$SL50 # k1 from SELECT
  SL95 <- FleetPars$SL95 #sigma (normal.loc, log-normal) or k2 (normal.sca) from SELECT
  MLLKnife <- FleetPars$MLLKnife
  
  if(sel == "Normal.loc"){                    # normal selectivity with fixed spread
    SLmesh <- FleetPars$SLmesh                # mesh sizes
    MLLNormal <- FleetPars$MLLNormal          # minimum landing limit, if necessary
    if (is.null(MLLNormal)) MLLNormal <- NA
  } else if(sel == "Normal.sca"){             # normal selectivity with proportional spread
    SLmesh <- FleetPars$SLmesh 
    MLLNormal <- FleetPars$MLLNormal
    if (is.null(MLLNormal)) MLLNormal <- NA 
  } else if(sel == "logNorm"){                # lognormal selectivity
    SLmesh <- FleetPars$SLmesh 
    MLLNormal <- FleetPars$MLLNormal
    if (is.null(MLLNormal)) MLLNormal <- NA  
  } else if(sel=="Knife"){                    # Knife-edge selectivity
    MLLKnife <- FleetPars$MLLKnife
  }  
  FM <- FleetPars$FM 
  
  
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
  
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 1
    SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  # Error Catches #
  if (!(exists("NGTG") | exists("GTGLinfBy"))) stop("NGTG or GTGLinfBy must be specified")
  if (!exists("R0")) R0 <- 1E6
  if (is.null(R0)) R0 <- 1E6
  
  # Set up Linfs for the different GTGs
  if (exists("NGTG") & !exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
    GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
  } else  if (!exists("NGTG") & exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
    NGTG <- length(DiffLinfs)
  } else if (exists("NGTG") & exists("GTGLinfBy")) {
    if (!is.na(GTGLinfBy)) {
      DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
      NGTG <- length(DiffLinfs)
    } 
    if (is.na(GTGLinfBy)) {
      DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
      GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
    }  
  } 
  # Distribute Recruits across GTGS 
  RecProbs <- dnorm(DiffLinfs, Linf, sd=SDLinf) / 
    sum(dnorm(DiffLinfs, Linf, sd=SDLinf)) 
  
  # Length Bins 
  if (is.null(ToSize)) ToSize <- max(DiffLinfs, Linf + MaxSD * SDLinf)
  LenBins <- seq(from=0, by=Linc, to=ToSize)
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=(length(LenBins)-1))
  
  Weight <- Walpha * LenMids^Wbeta
  
  # Maturity and Fecundity for each GTG 
  L50GTG <- L50/Linf * DiffLinfs # Maturity at same relative size
  L95GTG <- L95/Linf * DiffLinfs # Assumes maturity age-dependant 
  DeltaGTG <- L95GTG - L50GTG
  MatLenGTG <- sapply(seq_along(DiffLinfs), function (X) 
    1.0/(1+exp(-log(19)*(LenMids-L50GTG[X])/DeltaGTG[X])))
  FecLenGTG <- MatLenGTG * LenMids^FecB # Fecundity across GTGs 
  
  if(sel=="Logistic"){
    VulLen <- 1.0/(1+exp(-log(19)*((LenBins+0.5*Linc)-SL50)/((SL95)-(SL50)))) # Selectivity-at-Length
    
  }else if(sel=="Normal.sca"){ 
    VulLen <- 0
    for (j in seq_along(SLmesh)){
      VulLen <- VulLen + exp(-0.5*(((LenBins+0.5*Linc)-((SL50)*SLmesh[j]))/((SL95)^0.5*SLmesh[j]))^2)
    }
    if(!is.na(MLLNormal)) VulLen[LenBins < MLLNormal] <- 0 
    VulLen <- VulLen/max(VulLen)
    
  }else if(sel=="Normal.loc"){ 
    VulLen <- 0
    for (j in seq_along(SLmesh)){
      VulLen <- VulLen + exp(-0.5*(((LenBins+0.5*Linc)-((SL50)*SLmesh[j]))/((SL95)))^2)
    }
    if(!is.na(MLLNormal)) VulLen[LenBins < MLLNormal] <- 0 
    VulLen <- VulLen/max(VulLen)
    
  }else if(sel=="logNorm"){ 
    VulLen <- 0
    for (j in seq_along(SLmesh)){
      VulLen <- VulLen + exp(-0.5*((log(LenBins+0.5*Linc)-log((SL50)*SLmesh[j]))/(SL95))^2)
    }
    if(!is.na(MLLNormal)) VulLen[LenBins < MLLNormal] <- 0 
    VulLen <- VulLen/max(VulLen)
    
  }else if(sel=="Knife"){    # knife-edge selectivity
    VulLen <- 0
    VulLen[(LenBins+0.5*Linc) < MLLKnife] <- 0
    VulLen[(LenBins+0.5*Linc) > MLLKnife] <- 1
    SL95 <- SL50 <- NA 
  }
  

  # Add F-mortality below MLL
  SelLen <- VulLen # Selectivity is equal to vulnerability currently
  
  # Life-History Ratios 
  MKL <- MK * (Linf/(LenBins+0.5*Linc))^Mpow # M/K ratio for each length class
  # Matrix of MK for each GTG
  # MKMat <- sapply(seq_along(DiffLinfs), function(X) 
  # MKL + Mslope*(DiffLinfs[X] - CentLinf))
  MKMat <- matrix(rep(MKL, NGTG), nrow=length(MKL), byrow=FALSE)
  
  FK <- FM * MK # F/K ratio 
  FKL <- FK * SelLen # F/K ratio for each length class   
  # FkL[Legal == 0] <- FkL[Legal == 0] * DiscardMortFrac 
  ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG
  
  # Set Up Empty Matrices 
  # number-per-recruit at length
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) 
  
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- 
    NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), 
                                                ncol=NGTG) # number per GTG in each length class 
  # Distribute Recruits into first length class
  NPRFished[1, ] <- NPRUnfished[1, ] <- RecProbs * R0 
  for (L in 2:length(LenBins)) { # Calc number at each size class
    NPRUnfished[L, ] <- NPRUnfished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^MKMat[L-1, ]
    NPRFished[L, ] <- NPRFished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^ZKLMat[L-1, ]
    ind <- DiffLinfs  < LenBins[L]
    NPRFished[L, ind] <- 0
    NPRUnfished[L, ind] <- 0
  } 
  NPRUnfished[is.nan(NPRUnfished)] <- 0
  NPRFished[is.nan(NPRFished)] <- 0
  NPRUnfished[NPRUnfished < 0] <- 0
  NPRFished[NPRFished < 0] <- 0
  
  for (L in 1:length(LenMids)) { # integrate over time in each size class
    NatLUnFishedPop[L, ] <- (NPRUnfished[L,] - NPRUnfished[L+1,])/MKMat[L, ]
    NatLFishedPop[L, ] <- (NPRFished[L,] - NPRFished[L+1,])/ZKLMat[L, ]  
    FecGTGUnfished[L, ] <- NatLUnFishedPop[L, ] * FecLenGTG[L, ]
  }
  
  if(sel=="Logistic"){
    VulLen2 <- 1.0/(1+exp(-log(19)*(LenMids-(SL50))/((SL95)-(SL50))))# Selectivity-at-Length
    
  }else if(sel=="Normal.sca"){   # normal selectivity with proportional spread   
    VulLen2 <- 0
    for (j in seq_along(SLmesh)){
      VulLen2 <- VulLen2 + exp(-0.5*((LenMids-((SL50)*SLmesh[j]))/((SL95)^0.5*SLmesh[j]))^2)
    }
    if(!is.na(MLLNormal)) VulLen2[LenMids < MLLNormal] <- 0 
    VulLen2 <- VulLen2/max(VulLen2)
    
  }else if(sel=="Normal.loc"){   # normal selectivity with fixed spread
    VulLen2 <- 0
    for (j in seq_along(SLmesh)){
      VulLen2 <- VulLen2 + exp(-0.5*((LenMids-((SL50)*SLmesh[j]))/(SL95))^2)
    }
    if(!is.na(MLLNormal)) VulLen2[LenMids < MLLNormal] <- 0 
    VulLen2 <- VulLen2/max(VulLen2)
    
  }else if(sel=="logNorm"){   # lognormal selectivity
    VulLen2 <- 0
    for (j in seq_along(SLmesh)){
      VulLen2 <- VulLen2 + exp(-0.5*((log(LenMids)-log((SL50)*SLmesh[j]))/(SL95))^2)
    }
    if(!is.na(MLLNormal)) VulLen2[LenMids < MLLNormal] <- 0 
    VulLen2 <- VulLen2/max(VulLen2)
    
  }else if(sel=="Knife"){   # knife-edge selectivity
    VulLen2 <- 0
    VulLen2[LenMids < MLLKnife] <- 0
    VulLen2[LenMids > MLLKnife] <- 1
    SL95 <- SL50 <- NA
  }
  
  
  
  # print(LenMids)
  # print(c(SL50, SL95))
  # plot(LenMids, VulLen2)
  
  
  # print(cbind(LenMids, VulLen2))
  NatLUnFishedCatch <- NatLUnFishedPop * VulLen2 # Unfished Vul Pop
  NatLFishedCatch <- NatLFishedPop * VulLen2 # Catch Vul Pop
  
  # plot(LenMids, apply(NatLFishedCatch, 1, sum), type="p")
  # matplot(LenMids, (NatLFishedCatch), type="l")
  
  # Expected Length Structure - standardised 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum))
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum))
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))
  
  # Calc SPR
  EPR0 <- sum(NatLUnFishedPop * FecLenGTG) # Eggs-per-recruit Unfished
  EPRf <- sum(NatLFishedPop * FecLenGTG) # Eggs-per-recruit Fished
  SPR <- EPRf/EPR0 
  
  # Equilibrium Relative Recruitment
  recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
  # RelRec/R0 - relative recruitment 
  YPR <- sum(NatLFishedPop  * Weight * VulLen2) * FM 
  Yield <- YPR * RelRec
  
  # Calc Unfished Fitness 
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per Group
  FitPR <- Fit/RecProbs # Fitness per-recruit
  FitPR <- FitPR/median(FitPR)
  ## Debugging
  # plot(FitPR, ylim=c(0,2)) # Should be relatively flat for equal fitness across GTG
  
  # Mslope ignored in this version 
  ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting Mslope 
  Pen <- 0; if (min(MKMat) <= 0 ) Pen <- (1/abs(min(MKMat)))^2 * 1E12 # Penalty for optimising Mslope
  
  ObjFun <- ObjFun + Pen
  # print(cbind(Mslope, ObjFun, Pen))
  
  # Calculate spawning-per-recruit at each size class
  SPRatsize <- cumsum(rowSums(NatLUnFishedPop * FecLenGTG))
  SPRatsize <- SPRatsize/max(SPRatsize)
  
  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$YPR <- YPR
  Output$LCatchFished <- ExpectedLenCatchFished
  Output$LPopFished <- ExpectedLenPopFished
  Output$LCatchUnfished <- ExpectedLenCatchUnfished
  Output$LPopUnfished <- ExpectedLenPopUnfished
  Output$NatLPopFished <- NatLFishedPop
  Output$NatLPopUnFish <- NatLUnFishedPop
  Output$NatLCatchUnFish <- NatLUnFishedCatch
  Output$NatLCatchFish <- NatLFishedCatch
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- DiffLinfs[2] - DiffLinfs[1]
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$Winf <- Walpha * Linf^Wbeta
  Output$FecLen <- FecLenGTG 
  Output$MatLen <- MatLenGTG 
  Output$SelLen <- SelLen
  Output$MKL <- MKL
  Output$MKMat <- MKMat 
  Output$FKL <- FKL 
  Output$ZKLMat <- ZKLMat 
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$FitPR <- FitPR
  Output$Diff <- range(FitPR)[2] - range(FitPR)[1]
  Output$L50GTG <- L50GTG 
  Output$L95GTG <- L95GTG
  Output$SPRatsize <- SPRatsize
  Output$RelRec <- RelRec
  return(Output)
}

##########################
# Optimisation Functions #
##########################

#OptFun calls the length dist simulation and SPR calculations for given M/K, Linf and sigma^2 or CV Linf and trial selectivity parameters plus F/M
#the negative log likelihood (multinomial likelihood) of the observed length distn given life-history and selectivity parameters is returned

OptFunDome <- function(tryFleetPars, LenDat, fixedFleetPars, StockPars, SizeBins=NULL, 
                     mod=c("GTG", "LBSPR"),sel=c("Logistic","Normal.sca","Knife", "logNorm", "Normal.loc")) {
  
  Fleet <- NULL
  if(sel=="Logistic"){
    Fleet$SL50 <- fixedFleetPars$SL50  # exp(tryFleetPars[1]) * StockPars$Linf
    Fleet$SL95 <- fixedFleetPars$SL95  # Fleet$SL50  + (exp(tryFleetPars[2]) * StockPars$Linf)
  }else if(sel=="Normal.sca"){
    Fleet$SL50 <- fixedFleetPars$SL50 
    Fleet$SL95 <- fixedFleetPars$SL95
    Fleet$SLmesh <- fixedFleetPars$SLmesh
    if(!is.null(fixedFleetPars$MLLNormal)) Fleet$MLLNormal <- fixedFleetPars$MLLNormal
  }else if(sel=="Normal.loc"){
    Fleet$SL50 <- fixedFleetPars$SL50 
    Fleet$SL95 <- fixedFleetPars$SL95
    Fleet$SLmesh <- fixedFleetPars$SLmesh
    if(!is.null(fixedFleetPars$MLLNormal)) Fleet$MLLNormal <- fixedFleetPars$MLLNormal
  }else if(sel=="logNorm"){
    Fleet$SL50 <- fixedFleetPars$SL50 
    Fleet$SL95 <- fixedFleetPars$SL95
    Fleet$SLmesh <- fixedFleetPars$SLmesh
    if(!is.null(fixedFleetPars$MLLNormal)) Fleet$MLLNormal <- fixedFleetPars$MLLNormal
  }else if(sel=="Knife"){
    Fleet$MLLKnife <- fixedFleetPars$MLLKnife
  }
  
  
  Fleet$FM <- exp(tryFleetPars[1]) # changed to 1 from 3, as other parameters are fixed
  
  if (mod == "GTG") runMod <-  GTGDomeLBSPRSim(StockPars, Fleet, SizeBins, sel)
  if (mod == "LBSPR") runMod <- DomeLBSPRSim(StockPars, Fleet, SizeBins, sel)
  
  LenDat <- LenDat + 1E-15 # add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat)
  predProb <- runMod$LCatchFished 
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb))

  #Penalty functions are used when model estimates the selectivity  
  # add penalty for SL50
  #trySL50 <- exp(tryFleetPars[1])
  #PenVal <- NLL
  #Pen <- dbeta(trySL50, shape1=5, shape2=0.01) * PenVal  #penalty for trySL50 values close to 1/SL50 close to Linf 
  #if(!is.finite(NLL)) return(1E9 + runif(1, 1E4, 1E5))
  #if (Pen == 0) {Pen <- PenVal * trySL50
}
# plot(xx, dbeta(xx, shape1=5, shape2=0.01) )


#}

#NLL <- NLL+Pen

#return(NLL)
#}

DoOptDome <- function(StockPars, fixedFleetPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR"),
                    sel=c("Logistic","Normal.sca","Knife", "logNorm", "Normal.loc")) {
  
  SDLinf <- StockPars$CVLinf * StockPars$Linf
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 1
    SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
    SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  # Starting guesses
  # sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf
  # sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
  sFM <- 0.5 
  
  if(sel=="Logistic"){
    Start <- log(c(sFM))#sSL50, sDel, sFM))  #tryFleetPars
  }else if(sel=="Normal.sca"){
    Start <- log(c(sFM))  #tryFleetPars
  }else if(sel=="Normal.loc"){
      Start <- log(c(sFM))  #tryFleetPars
  }else if (sel=="logNorm"){
    Start <- log(c(sFM)) #tryFleetPars
  } else if (sel=="Knife"){
    Start <- log(c(sFM)) #tryFleetPars
  }
  
  
  opt <- nlminb(Start, OptFunDome, LenDat=LenDat, 
                fixedFleetPars=fixedFleetPars, StockPars=StockPars, 
                SizeBins=SizeBins, mod=mod, sel=sel, 
                control= list(iter.max=300, eval.max=400, abs.tol=1E-20))
  
  newFleet <- NULL 
  newFleet$FM <- exp(opt$par[1]) #changed to 1, previoulsy 3
  newNLL<-opt$objective

  
  if(sel=="Logistic"){
    newFleet$SL50 <- fixedFleetPars$SL50  # exp(opt$par[1]) * StockPars$Linf 
    newFleet$SL95 <- fixedFleetPars$SL95  # newFleet$SL50 + exp(opt$par[2]) * StockPars$Linf
  }else if(sel=="Normal.sca"){
    newFleet$SL50 <- fixedFleetPars$SL50     # prescribed values, not optimised
    newFleet$SL95 <- fixedFleetPars$SL95     # prescribed values, not optimised
    newFleet$SLmesh <- fixedFleetPars$SLmesh # prescribed values, not optimised
    if(!is.null(fixedFleetPars$MLLNormal)) newFleet$MLLNormal <- fixedFleetPars$MLLNormal
  }else if(sel=="logNorm"){
    newFleet$SL50 <- fixedFleetPars$SL50     # prescribed values, not optimised
    newFleet$SL95 <- fixedFleetPars$SL95     # prescribed values, not optimised
    newFleet$SLmesh <- fixedFleetPars$SLmesh # prescribed values, not optimised
    if(!is.null(fixedFleetPars$MLLNormal)) newFleet$MLLNormal <- fixedFleetPars$MLLNormal
  }else if(sel=="Normal.loc"){
    newFleet$SL50 <- fixedFleetPars$SL50     # prescribed values, not optimised
    newFleet$SL95 <- fixedFleetPars$SL95     # prescribed values, not optimised
    newFleet$SLmesh <- fixedFleetPars$SLmesh # prescribed values, not optimised
    if(!is.null(fixedFleetPars$MLLNormal)) newFleet$MLLNormal <- fixedFleetPars$MLLNormal
  }else if(sel=="Knife"){
    newFleet$MLLKnife <- fixedFleetPars$MLLKnife
  }
  
  if (mod == "GTG") runMod <-  GTGDomeLBSPRSim(StockPars, newFleet, SizeBins, sel)
  if (mod == "LBSPR") runMod <- DomeLBSPRSim(StockPars, newFleet, SizeBins, sel)
  
  Out <- NULL 
  Out$Ests <- c(FM=newFleet$FM, SL50=newFleet$SL50, SL95=newFleet$SL95, 
                SLmesh=newFleet$SLmesh,SPR=runMod$SPR)
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  Out$NLL<-newNLL
  Out$opt_par<-opt$par
  return(Out)
}
