######################################
# did package
# Callaway and Sant'Anna (2021)
# Difference-in-differences with multiple time periods 
# Journal of Econometrics
# Journal of Statistical Software, 225(2):200-230
######################################
library(did)

######################################
# geex package
# Saul and Hudgens (2020) 
# The calculus of M-estimation in R with geex
# Journal of Statistical Software, 92(2):1-15
######################################
library(geex)

source("UDID.R")
# Zika Virus in Brazil
# The source of the dataset are given below:
# Pre- and Post-treatment Outcomes, Treatment, and log population: 
#       zika_Table2.tab in  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY
# log population density and proportion of female
#       https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados

LongData <- read.csv("Zika_Brazil.csv")
# Y    = Outcome = birth rate
# Trt  = Indicator of whether a municipality belongs to Pernambuco (PE) (Trt=1) or belongs to Rio Grande do Sul (RS) (Trt=0)
# Time = 2014 (Time=0) or 2016 (Time=1)
# X_1  = log population
# X_2  = log population density
# X_3  = proportion of female
# ID   = municipality ID

################################################################################
# PT-based Estimators
################################################################################

PT.OR.Cov <- att_gt(yname="Y",
                    tname="Time",
                    idname="ID",
                    gname="Trt",
                    xformla=~X_1+X_2+X_3,
                    data=LongData,
                    bstrap = F,
                    est_method = "reg")
PT.PS.Cov <- att_gt(yname="Y",
                    tname="Time",
                    idname="ID",
                    gname="Trt",
                    xformla=~X_1+X_2+X_3,
                    data=LongData,
                    bstrap = F,
                    est_method = "ipw")
PT.DR.Cov <- att_gt(yname="Y",
                    tname="Time",
                    idname="ID",
                    gname="Trt",
                    xformla=~X_1+X_2+X_3,
                    data=LongData,
                    bstrap = F) 

PP <- function(ll){
  c(ll$att, ll$att-qnorm(0.975)*ll$se, ll$att+qnorm(0.975)*ll$se)
}

Result.PT <- rbind( c(PP(PT.OR.Cov)), 
                    c(PP(PT.PS.Cov)),
                    c(PP(PT.DR.Cov)))

################################################################################
# OREC-based Estimators (UDiD)
################################################################################

Y0     <- LongData$Y[LongData$Time==0]
Y1     <- LongData$Y[LongData$Time==1]
Trt    <- LongData$Trt[LongData$Time==0]
X_1    <- LongData$X_1[LongData$Time==0]
X_2    <- LongData$X_2[LongData$Time==0]
X_3    <- LongData$X_3[LongData$Time==0]

X_1    <- (X_1-min(X_1))/diff(range(X_1))
X_2    <- (X_2-min(X_2))/diff(range(X_2))

Data <- data.frame(cbind(Y0,Y1,Trt,X_1,X_2,X_3))
mY0 <- mean(Data$Y0)

Data$Y0 <- Data$Y0-mY0
Data$Y1 <- Data$Y1-mY0

Y0  <- Data$Y0
Y1  <- Data$Y1
Trt <- Data$Trt

Crude <- mean(Y1[Trt==1]) - mean(Y1[Trt==0])

Design.Matrix.OdRa <- Design.Matrix.Time0 <- Design.Matrix.Time1 <- 
  model.matrix(~X_1+X_2+X_3,data=Data)

UDID.OR.Cov <- UDID.OR(Y0,
                       Y1,
                       Trt,
                       Design.Matrix.Time0 = Design.Matrix.Time0,
                       Design.Matrix.OdRa  = Design.Matrix.OdRa,
                       Design.Matrix.Time1 = Design.Matrix.Time1)

UDID.PS.Cov <- UDID.PS(Y0,
                       Y1,
                       Trt,
                       Design.Matrix.Time0 = Design.Matrix.Time0,
                       Design.Matrix.OdRa  = Design.Matrix.OdRa,
                       Design.Matrix.Time1 = Design.Matrix.Time1)

UDID.DR.Cov <- UDID.DR(Y0,
                       Y1,
                       Trt,
                       Design.Matrix.Time0 = Design.Matrix.Time0,
                       Design.Matrix.OdRa  = Design.Matrix.OdRa,
                       Design.Matrix.Time1 = Design.Matrix.Time1)
 
PP <- function(ll){
  c(ll$Est, ll$Est-qnorm(0.975)*ll$SE, ll$Est+qnorm(0.975)*ll$SE)
}

Result.OREC <- rbind( c(PP(UDID.OR.Cov)), 
                      c(PP(UDID.PS.Cov)),
                      c(PP(UDID.DR.Cov)))

Result.PT
Result.OREC

################################################################################
# Assessment of Covariate Distribution Invariance
################################################################################

Y0     <- LongData$Y[LongData$Time==0]
Y1     <- LongData$Y[LongData$Time==1]
Trt    <- LongData$Trt[LongData$Time==0]
X_1    <- LongData$X_1[LongData$Time==0]
X_2    <- LongData$X_2[LongData$Time==0]
X_3    <- LongData$X_3[LongData$Time==0]

Data <- data.frame(cbind(Y0,Y1,Trt,X_1,X_2,X_3))
Data$Y0 <- Data$Y0
Data$Y1 <- Data$Y1

Y0  <- Data$Y0
Y1  <- Data$Y1
Trt <- Data$Trt

# Test

Tstat <- matrix(0,3,3)

for(tt in 1:3){
  if(tt==1){
    Xt  <- X_1 
  } else if (tt==2){
    Xt  <- X_2
  } else if (tt==3){
    Xt  <- X_3
  }
  
  
  Xtc <- Xt[Data$Trt==0]
  Y0c <- Y0[Data$Trt==0]
  Y1c <- Y1[Data$Trt==0]
  ND <- data.frame(Xtc=Xtc,
                   Y0c=Y0c,
                   Y1c=Y1c)
  
  est <- function(data){
    Xtc <- data$Xtc
    Y0c <- data$Y0c
    Y1c <- data$Y1c
    function(theta){
      c(Xtc - theta[1] - theta[2]*Y0c - theta[3]*Y0c^2,
        Y0c*(Xtc - theta[1] - theta[2]*Y0c - theta[3]*Y0c^2),
        Y0c^2*(Xtc - theta[1] - theta[2]*Y0c - theta[3]*Y0c^2),
        Xtc - theta[4] - theta[5]*Y1c - theta[6]*Y1c^2,
        Y1c*(Xtc - theta[4] - theta[5]*Y1c - theta[6]*Y1c^2),
        Y1c^2*(Xtc - theta[4] - theta[5]*Y1c - theta[6]*Y1c^2))
    }
  }
  
  C0 <- summary(lm(Xtc~Y0c+I(Y0c^2)))$coefficients
  C1 <- summary(lm(Xtc~Y1c+I(Y1c^2)))$coefficients
  
  ME <- m_estimate(
    estFUN = est,
    data   = ND,
    root_control = setup_root_control(start = c(C0[,1],C1[,1])))
  
  ME@rootFUN_results
  ME@vcov
  
  contrast.int <- c(1,0,0,-1,0,0)
  Num.int <- sum(contrast.int*ME@rootFUN_results$root)
  Den.int <- sqrt(t(contrast.int)%*%(ME@vcov)%*%(contrast.int))
  
  Tstat[1,tt] <- Num.int/Den.int
  
  contrast.Y <- c(0,1,0,0,-1,0)
  Num.Y <- sum(contrast.Y*ME@rootFUN_results$root)
  Den.Y <- sqrt(t(contrast.Y)%*%(ME@vcov)%*%(contrast.Y))
  
  Tstat[2,tt] <- Num.Y/Den.Y
  
  contrast.Y <- c(0,0,1,0,0,-1)
  Num.Y <- sum(contrast.Y*ME@rootFUN_results$root)
  Den.Y <- sqrt(t(contrast.Y)%*%(ME@vcov)%*%(contrast.Y))
  
  Tstat[3,tt] <- Num.Y/Den.Y
  
}

Tstat^2

# plot

MM <- matrix(1:6,3,2,byrow=T)
MM <- rbind(c(7,8), MM)
MM <- cbind(9:12, MM)
W  <- 6

layout(MM,
       heights=c(1,W,W,W),
       widths=c(1,W,W))
MAR <- c(3,3,1,0.5)
par(mar=MAR)

for(tt in 1:3){
  
  YL  <- c(2.5,22.5)
  
  if(tt==1){
    Xt  <- X_1
    XL <- c(7,14.5)
  } else if (tt==2){
    Xt  <- X_2
    XL <- c(0,8.25)
  } else if (tt==3){
    Xt  <- X_3
    XL <- c(0.44,0.54)
  }
  
  Xtc <- Xt[Data$Trt==0]
  Y0c <- Y0[Data$Trt==0]
  Y1c <- Y1[Data$Trt==0]
  
  PP <- function(x,y){
    plot(y,x,xlim=YL,ylim=XL,pch=19,cex=0.5,col=rgb(0,0,0,0.4))
  }
  
  PP(Xtc, Y0c)
  PP(Xtc, Y1c)
  
  
  C0 <- summary(lm(Xtc~Y0c))$coefficients
  C1 <- summary(lm(Xtc~Y1c))$coefficients
  
  cbind(C0[,1],C0[,1]-1.96*C0[,2],C0[,1]+1.96*C0[,2])
  cbind(C1[,1],C1[,1]-1.96*C1[,2],C1[,1]+1.96*C1[,2])
  
}

par(mar=MAR*c(0,1,0,1))
plot.new()
text(0.5,0.5,expression("Pre-treatment Birth Rate ("*Y["0"]*")"),cex=1.5)
plot.new()
text(0.5,0.5,expression("Post-treatment Birth Rate ("*Y["1"]*")"),cex=1.5)

par(mar=MAR*c(0,0,0,0))
plot.new() # blank

par(mar=MAR*c(1,0,1,0))
plot.new()
text(0.3,0.5,expression("Log Population"),cex=1.5,srt=90)
text(0.7,0.5,expression("("*X["1"]*")"),cex=1.5,srt=90)
plot.new()
text(0.3,0.5,expression("Population Density"),cex=1.5,srt=90)
text(0.7,0.5,expression("("*X["2"]*")"),cex=1.5,srt=90)
plot.new()
text(0.3,0.5,expression("Proportion of Females"),cex=1.5,srt=90)
text(0.7,0.5,expression("("*X["3"]*")"),cex=1.5,srt=90)


################################################################################
# OREC-based Estimators (UDiD) with discretized odds ratio function
################################################################################

Y0     <- LongData$Y[LongData$Time==0]
Y1     <- LongData$Y[LongData$Time==1]
Trt    <- LongData$Trt[LongData$Time==0]
X_1    <- LongData$X_1[LongData$Time==0]
X_2    <- LongData$X_2[LongData$Time==0]
X_3    <- LongData$X_3[LongData$Time==0]

X_1    <- (X_1-min(X_1))/diff(range(X_1))
X_2    <- (X_2-min(X_2))/diff(range(X_2))

Data <- data.frame(cbind(Y0,Y1,Trt,X_1,X_2,X_3))
mY0 <- mean(Data$Y0)

Data$Y0 <- Data$Y0-mY0
Data$Y1 <- Data$Y1-mY0

Y0  <- Data$Y0
Y1  <- Data$Y1
Trt <- Data$Trt

Design.Matrix.Time0 <- Design.Matrix.Time1 <- model.matrix(~1+X_1+X_2+X_3,data=Data)
Design.Matrix.OdRa <- model.matrix(~1,data=Data)

num.bin <- 10
Bin.Cut    <- quantile(Y0,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
Disc.Y0    <- apply(matrix(Y0,length(Y0),1),1,function(v){sum(v>Bin.Cut)})
Disc.Y1    <- apply(matrix(Y1,length(Y1),1),1,function(v){sum(v>Bin.Cut)})

Design.Matrix.OdRa.Y0.Full <- 
  Design.Matrix.OdRa.Y1.Full <- matrix(0,length(Disc.Y0),num.bin)
for(jj in 1:num.bin){
  Design.Matrix.OdRa.Y0.Full[Disc.Y0==jj-1,jj] <- 1
  Design.Matrix.OdRa.Y1.Full[Disc.Y1==jj-1,jj] <- 1
}


UDID.OR.Cov.Disc.10 <- UDID.OR.Disc(Y0,
                                    Y1,
                                    Trt,
                                    Design.Matrix.Time0 = Design.Matrix.Time0,
                                    Design.Matrix.OdRa = Design.Matrix.OdRa,
                                    Design.Matrix.Time1 = Design.Matrix.Time1,
                                    num.bin = num.bin)

UDID.PS.Cov.Disc.10 <- UDID.PS.Disc(Y0,Y1,Trt,
                                    Design.Matrix.Time0,
                                    Design.Matrix.OdRa.Y0.Full,
                                    Design.Matrix.Time1,
                                    Design.Matrix.OdRa.Y1.Full,
                                    num.bin=num.bin)

UDID.DR.Cov.Disc.10 <- UDID.DR.Disc(Y0,Y1,Trt,
                                    Design.Matrix.Time0,
                                    Design.Matrix.OdRa,
                                    Design.Matrix.OdRa.Y0.Full,
                                    Design.Matrix.Time1,
                                    Design.Matrix.OdRa.Y1.Full,
                                    num.bin=num.bin)


PP <- function(ll){
  c(ll$Est, ll$Est-qnorm(0.975)*ll$SE, ll$Est+qnorm(0.975)*ll$SE)
}


Result.OREC.Disc <- rbind( c(PP(UDID.OR.Cov.Disc.10)),  
                           c(PP(UDID.PS.Cov.Disc.10)), 
                           c(PP(UDID.DR.Cov.Disc.10)) )


