#################################
# Model run with empirical data #
#################################

library(plyr) 
library(dplyr)
library(ggplot2)
library(rlist)

# load GTG LB-SPR routines
source("GTGLBSPR_Dome.R")

# read in the input data ----

# brown trout length data, Glen Lough ====
#setwd("")
trout <- read.csv(file="Trout_Selection.csv", header=TRUE)
catchAtLength <- trout %>% filter(lake=="Glen")  
empLengthData <- catchAtLength$TLcm  # empirical length data of brown trout in Glen Lough
lengthBinWidth <- 1


# visualise length distribution ====
# geom_histogram
pgLHist <- ggplot(data = catchAtLength) + 
  geom_histogram(aes(x = TL), binwidth = lengthBinWidth, boundary = 0, closed = "left", 
                 fill = "grey75", colour = "black") + 
  theme_bw()
pgLHist

# hist
lHist <- hist(empLengthData, breaks = seq(0, ceiling(max(empLengthData)), lengthBinWidth), col="grey",  
     xlab = "Length (cm)", xlim=c(0, 1.2*ceiling(max(empLengthData))), right = FALSE,
     main="Catch length composition", plot = TRUE) 
# samples_glen <- readRDS("samples_glen.rds") # reading in the Linf posterior distribution



# GTG LB-SPR assessment ----

# Brown trout in Glen - stock parameters ====
StockPars <- NULL 
# StockPars$MK <- 1.8
StockPars$NGTG <- 17
StockPars$Linf <- 40.7
#StockPars$Linf <- 120.
StockPars$CVLinf <- 0.1
StockPars$MaxSD <- 2
StockPars$L50 <- 13.834  
StockPars$L95 <- 24.577 
#StockPars$L50 <- 45.31
#StockPars$L95 <- 60
StockPars$FecB <- 3 
StockPars$Walpha <- 0.0084
StockPars$Wbeta <- 3.1115
StockPars$Steepness <- 0.8 
StockPars$Mpow <- 0

# gear selectivity ====
FleetPars <- NULL

# see Table 4 
# "Dome-shaped selectivity in LB-SPR: Length-Based assessment of data-limited inland fish stocks sampled with gillnets"
# https://www.sciencedirect.com/science/article/pii/S0165783620300916

# lognorm
FleetPars$SL1 <- 10.43     # k_mode from SELECT (log-normal)
FleetPars$SL2 <- 0.27      # standard deviation at log scale from SELECT
gearSelectivity <- "logNorm"

# norm.loc
#FleetPars$SL1 <- 9.52
#FleetPars$SL2 <- 7.06
#gearSelectivity <- "Normal.loc"

FleetPars$SLmesh <- c(1.25, 1.55, 1.95, 2.4, 2.9, 3.5, 4.3, 5.5) # used mesh sizes
FleetPars$MLLNormal <- 23   # minimum landing limit (MLL)
#FleetPars$MLLKnife <- 23

# logistic
# FleetPars <- NULL
# gearSelectivity <- "Logistic"
# FleetPars$SL1 <- 75.0
# FleetPars$SL2 <- 90.0


# preliminary visualisation of selectivity ####

# calculate selectivity
lengthFish <- seq(0, StockPars$Linf*(1 + StockPars$CVLinf*StockPars$MaxSD), length.out = 101)
if(!is.null(FleetPars$SLmesh)) meshSize <- FleetPars$SLmesh

if(gearSelectivity == "Logistic"){
  SL50 <- FleetPars$SL1; SL95 <- FleetPars$SL2
  gearSelLen <- 1.0/(1+exp(-log(19)*(lengthFish-SL50)/((SL95)-(SL50)))) # Selectivity-at-Length
}else if(gearSelectivity=="Normal.sca"){
  SLk1 <- FleetPars$SL1; SLk2 <- FleetPars$SL2; MLLNormal <- FleetPars$MLLNormal
  gearSelLen <- 0
  for (j in seq_along(SLmesh)){
    gearSelLen <- gearSelLen + exp(-0.5*((lengthFish-((SLk1)*meshSize[j]))/((SLk2)^0.5*meshSize[j]))^2)
  }
  if(!is.na(MLLNormal)) gearSelLen[lengthFish < MLLNormal] <- 0 
  gearSelLen <- gearSelLen/max(gearSelLen)
  
}else if(gearSelectivity=="Normal.loc"){ 
  SLk <- FleetPars$SL1; SLsigma <- FleetPars$SL2; MLLNormal <- FleetPars$MLLNormal
  gearSelLen <- 0
  for (j in seq_along(meshSize)){
    gearSelLen <- gearSelLen + exp(-0.5*((lengthFish-((SLk)*meshSize[j]))/((SLsigma)))^2)
  }
  if(!is.na(MLLNormal)) gearSelLen[lengthFish < MLLNormal] <- 0 
  gearSelLen <- gearSelLen/max(gearSelLen)
  
}else if(gearSelectivity=="logNorm"){ 
  SLk <- FleetPars$SL1; SLsigma <- FleetPars$SL2; MLLNormal <- FleetPars$MLLNormal
  gearSelLen <- 0
  for (j in seq_along(meshSize)){
    gearSelLen <- gearSelLen + exp(-0.5*((log(lengthFish)-log((SLk)*meshSize[j]))/(SLsigma))^2)
  }
  if(!is.na(MLLNormal)) gearSelLen[lengthFish < MLLNormal] <- 0 
  gearSelLen <- gearSelLen/max(gearSelLen)
  
}else if(gearSelectivity=="Knife"){    # knife-edge selectivity
  MLLKnife <- FleetPars$MLLKnife
  gearSelLen <- 0
  gearSelLen[lengthFish < MLLKnife] <- 0
  gearSelLen[lengthFish > MLLKnife] <- 1
}

# add selectivity line to histogram
lines(lengthFish, max(lHist$counts)*gearSelLen, lty =2, col = "black", lwd = 1.5)
pgLHist <- pgLHist + geom_line(data = data.frame(length = lengthFish, selectivity = max(lHist$counts)*gearSelLen),
                               aes(x = length, y = selectivity), colour = "black", linetype = 2, size = 1.25)
pgLHist

# gear/fleet parameters  are considered "fixed" - not optmised when fitting to data
fixedFleetPars <- FleetPars


# apply LB-SPR: test run ====

# length data - discretisation ####
SizeBins <- NULL
SizeBins$Linc <- lengthBinWidth
SDLinf <- StockPars$CVLinf * StockPars$Linf
SizeBins$ToSize <- StockPars$Linf + SDLinf*StockPars$MaxSD # maximum length of length bins based on GTG pars

LenBins <- seq(from=0, to=SizeBins$ToSize, by=SizeBins$Linc)
LenMids <- seq(from=0.5*SizeBins$Linc, by=SizeBins$Linc,length.out=(length(LenBins)-1))
LenDat <- as.vector(table(cut(unlist(empLengthData), LenBins))) 

# add length data used in LBSPR to histogram
points(LenMids, as.vector(table(cut(unlist(empLengthData), LenBins, right = FALSE))), pch = 1)

# test run LBSPR ####
# estimation routine DoOpt
StockPars$MK <- 1.8
testOpt <- DoOptDome(StockPars, fixedFleetPars, LenDat, SizeBins, "GTG", gearSelectivity)

# estimation outputs (note SLmesh, SL1, SL2 may be fixed)
testOpt$Ests


# per recruit simulation based on FM etc. estimates
FleetPars <- list(FM = testOpt$Ests[["FM"]], SL1 = testOpt$Ests[["SL1"]], SL2 =testOpt$Ests[["SL2"]],
                  SLmesh = FleetPars$SLmesh, MLLNormal = FleetPars$MLLNormal)
prSim <- GTGDomeLBSPRSim(StockPars, FleetPars, SizeBins, sel = gearSelectivity)
sum(prSim$LCatchFished)

# show predicted length
lines(LenMids, max(lHist$counts)*prSim$LCatchFished/max(prSim$LCatchFished), lty = 1, col = "red")
#source("./GTG_LBSPR_CF_1.R")
#DoOpt(StockPars, fixedFleetPars, LenDat, SizeBins, "GTG") # GTGLBSPR_CF_1.R



# LB-SPR: apply to multiple M/K (and Linf) values ====
#nsims <- 1 #how many runs
MKseq <- c(1.5, 1.8, 2.0)
#MKseq <- c(1.8, 2.0, 2.2)
Linfseq  <- c(40.7, 42) # Just using a single estimate of Linf
#Linfseq  <- c(120) # Just using a single estimate of Linf
# Linf <- samples_glen$Linf # posterior distribution of Linf

# table to save output of the model
GTGOut <- array(NA, dim=c(length(Linfseq), 5, length(MKseq)), #dim=c(nsims,5,3,2), 
                dimnames = list(d1 = paste0("Linf", as.character(Linfseq)), #c("Linf40.7", "Linf42"),
                                d2 = c("NLL","FM","SL50", "SL95","SPR"),
                                d3 = c("MK1.5", "MK1.8", "MK2.0")
                                )
                ) 
# table to save the model input values
GTGIn <- array(NA, dim = c(length(Linfseq), 2, length(MKseq))) #c(nsims, length(Linfseq), 2, length(MKseq))


output_pred <- list()
# for (i in 1:nsims){
#   StockPars$Linf <- Linf[i]
for (iLinf in seq_along(Linfseq)){
  StockPars$Linf <- Linfseq[iLinf]
  for (iMK in seq_along(MKseq)){
    StockPars$MK <- MKseq[iMK]
    
    GTGIn[iLinf,1,iMK] <- StockPars$MK
    GTGIn[iLinf,2,iMK] <- StockPars$Linf
    
    runopt <- DoOptDome(StockPars, fixedFleetPars, LenDat=LenDat, SizeBins=SizeBins, mod="GTG", 
                        sel = gearSelectivity) # this is the model; specify the selectivity ("Normal.sca"/"Normal.loc"/logNorm") 
    
    GTGOut[iLinf,1,iMK] <- runopt$NLL[1]
    GTGOut[iLinf,2,iMK] <- runopt$Ests[1]  #FM
    GTGOut[iLinf,3,iMK] <- runopt$Ests[2]  #SL50
    GTGOut[iLinf,4,iMK] <- runopt$Ests[3]  #SL95
    GTGOut[iLinf,5,iMK] <- runopt$Ests["SPR"]  #SPR
    
    output_pred <- list.append(output_pred, runopt$PredLen) #saves predicted lengths for each run
  }
}

results_glen <- adply(GTGOut,c(1,3))
input_glen <- adply(GTGIn,c(1,3))

final_glen <- cbind(results_glen, input_glen)
