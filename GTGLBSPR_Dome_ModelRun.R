#################################
# Model run with empirical data #
#################################

library(plyr) 
library(dplyr)
library(ggplot2)
library(rlist)

# read in the input data
setwd("H:/RESEARCH/SSHEPHARD/IFI_Analysis/LBSPR_Trout/Kristiina/GTG_LBSPR")
trout <- read.csv(file="Trout_Selection.csv", header=TRUE)
glen <- trout %>% filter(lake=="Glen")  
EmpData <- glen$TLcm  # length data of brown trout in Glen Lough

# plotting the length distribution
hist(unlist(EmpData), breaks = 15, col="grey", main="Glen", xlab = "Length (cm)", xlim=c(5, 30)) 
# samples_glen <- readRDS("samples_glen.rds") # reading in the Linf posterior distribution


#Brown trout in Glen

StockPars <- NULL 
# StockPars$MK <- 1.8
StockPars$NGTG <- 17
# StockPars$Linf <-40.7
StockPars$CVLinf <- 0.1
StockPars$MaxSD <- 2
StockPars$L50 <- 13.834  
StockPars$L95 <- 24.577 
StockPars$FecB <- 3 
StockPars$Walpha <- 0.0084
StockPars$Wbeta <- 3.1115
StockPars$Steepness <- 0.8 
StockPars$Mpow <- 0

FleetPars <- NULL
FleetPars$SL50 <- 10.43     # k_mode from SELECT (log-normal)
FleetPars$SL95 <- 0.27      # standard deviation at log scale from SELECT
FleetPars$SLmesh <- c(1.25, 1.55, 1.95, 2.4, 2.9, 3.5, 4.3, 5.5) # used mesh sizes
FleetPars$MLLNormal <- 23   # minimum landing limit (MLL)
#FleetPars$MLLKnife <- 23

SizeBins <- NULL
SizeBins$Linc <- 1
SizeBins$ToSize <- 40.7 * 1.1 # This is StockPars$Linf * 1.1

File <- EmpData

LenBins <- seq(from=0, to=SizeBins$ToSize, by=SizeBins$Linc)
LenMids <- seq(from=0.05*SizeBins$Linc, by=SizeBins$Linc,length.out=(length(LenBins)-1))
LenDat <- as.vector(table(cut(unlist(EmpData), LenBins))) 

nsims <- 1 #how many runs
MKseq <- c(1.5, 1.8, 2.0)
Linfseq  <- c(40.7, 42) # Just using a single estimate of Linf
# Linf <- samples_glen$Linf # posterior distribution of Linf

# table to save output of the model
GTGOut <- array(NA, dim=c(nsims,5,3,2), 
                dimnames = list(d1 = c(1:nsims),d2 = c("NLL","FM","SL50", "sigma","SPR"),
                                d3 = c("MK1.5","MK1.8", "MK2.0"),
                                d4 = c("Linf40.7", "Linf42"))) 
# table to save the model input values
GTGIn <- array(NA, dim=c(nsims,2,6))  


output_pred <- list()
# for (i in 1:nsims){
#   StockPars$Linf <- Linf[i]
for (i in seq_along(Linfseq)){
  StockPars$Linf <- Linfseq[i]
  for (sp in seq_along(MKseq)){
    StockPars$MK <- MKseq[sp]
    
    GTGIn[i,1,sp] <- StockPars$MK
    GTGIn[i,2,sp] <- StockPars$Linf
    
    runopt <- DoOpt(StockPars, LenDat=LenDat, SizeBins=SizeBins, mod="GTG", sel = "Normal") # this is the model; specify the selectivity ("Normal"/"logNorm") 
    
    GTGOut[i,1,sp] <- runopt$NLL[1]
    GTGOut[i,2,sp] <- runopt$Ests[1]  #FM
    GTGOut[i,3,sp] <- runopt$Ests[2]  #SL50
    GTGOut[i,4,sp] <- runopt$Ests[3]  #sigma
    GTGOut[i,5,sp] <- runopt$Ests[4]  #SPR
    
    output_pred <- list.append(output_pred, runopt$PredLen) #saves predicted lengths for each run
  }
}

results_glen <- adply(GTGOut,c(1,3))
input_glen <- adply(GTGIn,c(1,3))

final_glen <- cbind(results_glen, input_glen)
