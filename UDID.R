library(geex)

expit <- function(v){ exp(v)/(1+exp(v)) }

LINK <- exp

######################################
# New
######################################

UDID.OR <- function(Y0,Y1,Trt,
                    Design.Matrix.Time0,
                    Design.Matrix.OdRa,
                    Design.Matrix.Time1){
  
  
  NV.Time0 <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1 <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Trt))
  X.Time1    <- as.matrix(cbind(Design.Matrix.Time1))
  
  data0.OdRa <- data.frame(cbind(Y0,X.Time0.OdRa))
  data1    <- data.frame(cbind(Y1,X.Time1))
  
  colnames(data0.OdRa) <- c("Y0",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  colnames(data1)    <- c("Y1",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  
  LM.Time0.OdRa <- lm(Y0~0+., data=data0.OdRa)
  LM.Time1    <- lm(Y1~0+., data=data1[Trt==0,])
  
  EE <- function(data){
    
    function(theta){
      
      Pred0 <- sum( ( theta[3+1:(NV.Time0+NV.OdRa)] ) * as.numeric( data[,3+1:(NV.Time0+NV.OdRa)] ) )
      Gap   <- sum( ( theta[3+NV.Time0+1:(NV.OdRa)] ) * as.numeric( data[,3+NV.Time0+NV.OdRa+1:NV.OdRa] ) )
      Res0  <- (data$Y0 - Pred0 )
      Pred1 <- sum( ( theta[3+NV.Time0+NV.OdRa+1:NV.Time1] ) * as.numeric( data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1] ) )
      Res1  <- (data$Y1 - Pred1 )
      
      c( (data$Trt)*( data$Y1-Pred1-theta[3]/theta[2]*(Gap) - theta[1] ), 
         (Res0)^2 - theta[2],
         (1-data$Trt)*((Res1)^2 - theta[3]),
         as.numeric(data[,3+1:(NV.Time0+NV.OdRa)])*Res0,
         (1-data$Trt)*as.numeric(data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1])*Res1)
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,
                           Design.Matrix.OdRa*Trt,
                           Design.Matrix.OdRa,
                           Design.Matrix.Time1))
  
  ParaStart <- c( mean(Y1[Trt==1])-mean(Y1[Trt==0])+
                    -(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                  mean((LM.Time0.OdRa$residuals)^2),
                  mean((LM.Time1$residuals)^2),
                  LM.Time0.OdRa$coefficients,
                  LM.Time1$coefficients )
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  # EESolve <- m_estimate(estFUN = EE, 
  #                  roots = ParaStart,
  #                  compute_roots=FALSE,
  #                  data  = data)
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "Var0",
                       "Var1",
                       sprintf("OR0_Base_%0.3d",1:NV.Time0),
                       sprintf("OR0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("OR1_Base_%0.3d",1:NV.Time1))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- names(CEE.Full)
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  return(Result)
  
}


UDID.PS <- function(Y0,Y1,Trt,
                    Design.Matrix.Time0,
                    Design.Matrix.OdRa,
                    Design.Matrix.Time1){
  
  
  
  NV.Time0 <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1 <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Y0))
  data0.OdRa   <- data.frame(cbind(Trt,X.Time0.OdRa))
  colnames(data0.OdRa) <- c("Trt",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  pi.Time0.OdRa <- glm(Trt~0+., data=data0.OdRa, family="binomial")
  
  X.Time1    <- as.matrix(cbind(Design.Matrix.Time1))
  data.Time1      <- data.frame(cbind(Trt,
                                      (Design.Matrix.OdRa*Y1)%*%pi.Time0.OdRa$coefficients[NV.Time0+1:NV.OdRa],
                                      X.Time1))
  colnames(data.Time1)    <- c("Trt","Base",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  pi.Time1 <- glm(Trt~0+., data=data.Time1,family="binomial")
  
  EE <- function(data){
    
    function(theta){
      
      Pred0 <- expit( sum( ( theta[3+1:(NV.Time0+NV.OdRa)] ) * as.numeric( data[,3+1:(NV.Time0+NV.OdRa)] ) ) )
      Res0 <- data$Trt - Pred0
      Gap   <- sum( ( theta[3+NV.Time0+1:(NV.OdRa)] ) * as.numeric( data[,3+NV.Time0+NV.OdRa+1:NV.OdRa] ) )
      OR0 <- exp( data$Y1*(Gap) )
      Res1 <- ( (1-data$Trt)*(1+exp( sum( ( theta[3+NV.Time0+NV.OdRa+1:NV.Time1] ) * 
                                            as.numeric( data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1] ) ))*OR0)-1 )
      IPW <- exp( sum( ( theta[3+NV.Time0+NV.OdRa+1:NV.Time1] ) * 
                         as.numeric( data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1] ) ))*OR0
      
      c( theta[1]-(theta[2]-theta[3]),
         (data$Trt)*(data$Y1 - theta[2]), 
         (1-data$Trt)*(data$Y1)*IPW - (1-data$Trt)*IPW*theta[3],
         as.numeric(data[,3+1:(NV.Time0+NV.OdRa)])*Res0 , 
         as.numeric(data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1])*Res1  )
      
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,
                           Design.Matrix.OdRa*Y0,
                           Design.Matrix.OdRa,
                           Design.Matrix.Time1))
  
  ParaStart <- c( mean(Y1[Trt==1])-mean(Y1[Trt==0])+
                    -(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                  mean(Y1[Trt==1]),
                  mean(Y1[Trt==0])+(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                  pi.Time0.OdRa$coefficients,
                  pi.Time1$coefficients[-1] )
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  # EESolve <- m_estimate(estFUN = EE, 
  #                  roots = ParaStart,
  #                  compute_roots=FALSE,
  #                  data  = data)
  
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("PS0_Base_%0.3d",1:NV.Time0),
                       sprintf("PS0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("PS1_Base_%0.3d",1:NV.Time1))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
}



UDID.DR <- function(Y0,Y1,Trt,
                    Design.Matrix.Time0,
                    Design.Matrix.OdRa,
                    Design.Matrix.Time1,
                    M.type="LinXLinY"){
  
  ## Moment Equation for OdRa
  
  NV.Time0   <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1   <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OR <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Trt))
  X.Time0.PS <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Y0))
  X.Time1.OR <- as.matrix(cbind(Design.Matrix.Time1))
  X.Time1.PS <- as.matrix(cbind(Design.Matrix.Time1))
  
  data.Time0.OR   <- data.frame(cbind(Y0, X.Time0.OR))
  data.Time1.OR   <- data.frame(cbind(Y1, X.Time1.OR))
  data.Time0.PS   <- data.frame(cbind(Trt,X.Time0.PS))
  
  colnames(data.Time0.OR) <- c("Y0",
                               sprintf("OR_%0.3d",1:dim(X.Time0.OR)[2]))
  colnames(data.Time1.OR) <- c("Y1",
                               sprintf("OR_%0.3d",1:dim(X.Time1.OR)[2]))
  colnames(data.Time0.PS) <- c("Trt",
                               sprintf("PS_%0.3d",1:dim(X.Time0.OR)[2]))
  
  OR.Time0  <- lm(Y0~0+.,  data=data.Time0.OR)
  OR.Time1  <- lm(Y1~0+.,  data=data.Time1.OR[Trt==0,])
  PS.Time0  <- glm(Trt~0+., data=data.Time0.PS,family="binomial")
  
  data.Time1.PS      <- data.frame(cbind(Trt,
                                         (Design.Matrix.OdRa*Y1)%*%PS.Time0$coefficients[NV.Time0+1:NV.OdRa],
                                         X.Time1.PS))
  colnames(data.Time1.PS) <- c("Trt","Base",
                               sprintf("PS_%0.3d",1:dim(X.Time1.PS)[2]))
  PS.Time1  <- glm(Trt~0+., data=data.Time1.PS,family="binomial")
  
  #########################################################################
  
  EE.Time0 <- function(data){
    
    jump.pos.t <- 1
    jump.pos.d <- 3
    pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    pos.PS.Time0.Full <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0+NV.OdRa)
    pos.PS.Time0.Base <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0)
    pos.OR.Time0.Base <- 1:(NV.Time0)
    pos.OR.Time1.Base <- (NV.Time0+NV.OdRa)+1:(NV.Time1)
    pos.OdRa          <- (2*NV.Time0+2*NV.OdRa+NV.Time1)+1:(NV.OdRa)
    
    function(theta){
      
      Theta <- c( mean(OR.Time0$residuals^2),
                  OR.Time0$coefficients,
                  OR.Time1$coefficients,
                  PS.Time0$coefficients,
                  theta)
      
      # Theta <- theta
      
      Pred.OR.Time0.Full       <- ( sum( ( Theta[jump.pos.t + pos.OR.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Full] ) ) )
      Pred.PS.Time0.Full       <- expit( sum( ( Theta[jump.pos.t + pos.PS.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Full] ) ) )
      
      Pred.OR.Base0 <- ( sum( ( Theta[jump.pos.t + pos.OR.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Base] ) ) )
      Pred.OR.Time1.Base       <- ( sum( ( Theta[jump.pos.t + pos.OR.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time1.Base] ) ) )
      Pred.PS.Base <- expit( sum( ( Theta[jump.pos.t + pos.PS.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Base] ) ) )
      
      Suff.OdRa           <- ( sum( ( Theta[jump.pos.t + pos.OdRa] ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
      Pred.OdRa.Time0     <- LINK(Suff.OdRa*(data$Y0-Pred.OR.Time1.Base))
      Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
      
      if(M.type=="LinXLinY"){
        M.matrix      <- as.numeric( data[,jump.pos.d + pos.OdRa] )*data$Y0
        Pred.M.Base0  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*Pred.OR.Base0
      } else if (M.type=="LinXExpY"){
        M.matrix      <- (as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(data$Y0)
        Pred.M.Base0  <- (as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(Theta[1]/2 + Pred.OR.Base0 )
      } else if (M.type=="ExpXLinY"){
        M.matrix      <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*(data$Y0)
        Pred.M.Base0  <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*Pred.OR.Base0
      } else if (M.type=="ExpXExpY"){
        M.matrix      <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(data$Y0)
        Pred.M.Base0  <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(Theta[1]/2 + Pred.OR.Base0 )
      } else if (M.type=="LogXLinY"){
        M.matrix      <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1  ))*(data$Y0)
        Pred.M.Base0  <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1  ))*Pred.OR.Base0
      } else if (M.type=="LogXExpY"){ 
        M.matrix      <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1  ))*exp(data$Y0)
        Pred.M.Base0  <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1  ))*exp(Theta[1]/2 + Pred.OR.Base0 )
      }  
      
      
      
      Res.OR.0 <- data$Y0  - Pred.OR.Time0.Full
      Res.PS.0 <- data$Trt - Pred.PS.Time0.Full
      
      Res.M.Base0      <- M.matrix - Pred.M.Base0
      Res.PS.0.Base   <- data$Trt - Pred.PS.Base
      c( ((Res.M.Base0)*( Sandwich.OdRa.Time0 ))*(Res.PS.0.Base)  )
      
      
    }
  }
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,         ## OR Time0
                           Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa
                           Design.Matrix.Time1,         ## OR Time1
                           Design.Matrix.Time0,         ## PS Time0
                           Design.Matrix.OdRa*Y0,       ## PS Time0, OdRa
                           Design.Matrix.OdRa))         ## OdRa
  
  theta <- ParaStart <- PS.Time0$coefficients[NV.Time0+(1:NV.OdRa)]
  
  EESolve.Time0 <- m_estimate(estFUN = EE.Time0,
                              data = data,
                              root_control = setup_root_control(start = ParaStart),
                              compute_vcov = FALSE)
  
  CEE0 <- c(mean(OR.Time0$residuals^2),                # OR Time0 Var
            OR.Time0$coefficients,                     # OR Time0
            PS.Time0$coefficients,
            coef(EESolve.Time0))
  
  
  ###############################################################################
  
  EE.Time1 <- function(data){
    
    jump.pos.t <- 5
    jump.pos.d <- 3
    pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    pos.PS.Time0.Full <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0+NV.OdRa)
    pos.OR.Time0.Base <- 1:(NV.Time0)
    pos.OR.Time1.Base <- NV.Time0+NV.OdRa+1:(NV.Time1)
    pos.PS.Time0.Base <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0)
    pos.PS.Time1.Base <- (2*NV.Time0+2*NV.OdRa+NV.Time1)+1:(NV.Time1)
    pos.OdRa          <- (2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)
    
    function(theta){
      
      Theta <- c(theta[c(1:3)],                                       # ATT, ATT1, ATT0
                 CEE0[1],                                             # Var(Y_{t=0})
                 mean(OR.Time1$residuals^2),                          # Var(Y_{t=1})
                 CEE0[1+1:(NV.Time0+NV.OdRa)],                        # OR.Full(t=0)
                 OR.Time1$coefficients,                               # OR(t=1)
                 CEE0[1+(NV.Time0+NV.OdRa)+1:(NV.Time0+NV.OdRa)],     # PS.Full(t=0)
                 theta[-c(1:3)],                                      # PS(t=1)
                 CEE0[1+(NV.Time0+NV.OdRa)*2+1:(NV.OdRa)])            # OdRa
      
      Pred.OR.Time0.Full       <- ( sum( ( Theta[jump.pos.t + pos.OR.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Full] ) ) )
      Pred.OR.Time1.Base       <- ( sum( ( Theta[jump.pos.t + pos.OR.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time1.Base] ) ) )
      Pred.PS.Time0.Full       <- expit( sum( ( Theta[jump.pos.t + pos.PS.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Full] ) ) )
      Pred.PS.Time1.Base.LinLk <- sum( ( Theta[jump.pos.t + pos.PS.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time1.Base] ) )
      
      Pred.OR.Base0 <- ( sum( ( Theta[jump.pos.t + pos.OR.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Base] ) ) )
      Pred.PS.Base <- expit( sum( ( Theta[jump.pos.t + pos.PS.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Base] ) ) )
      
      Suff.OdRa           <- ( sum( ( Theta[jump.pos.t + pos.OdRa] ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
      Pred.OdRa.Time0     <- LINK(Suff.OdRa*(data$Y0-Pred.OR.Time1.Base ))
      Pred.OdRa.Time1     <- LINK(Suff.OdRa*(data$Y1-Pred.OR.Time1.Base ))
      Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
      
      Res.OR.0 <- data$Y0  - Pred.OR.Time0.Full
      Res.OR.1 <- data$Y1  - Pred.OR.Time1.Base
      Res.PS.0 <- data$Trt - Pred.PS.Time0.Full
      Res.PS.1 <- ( (1-data$Trt)*(1+exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1  )-1 )
      
      if(M.type=="LinXLinY"){
        M.matrix      <- as.numeric( data[,jump.pos.d + pos.OdRa] )*data$Y0
        Pred.M.Base0  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*Pred.OR.Base0
      } else if (M.type=="LinXExpY"){
        M.matrix      <- (as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(data$Y0)
        Pred.M.Base0  <- (as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(Theta[4]/2 + Pred.OR.Base0 )
      } else if (M.type=="ExpXLinY"){
        M.matrix      <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*(data$Y0)
        Pred.M.Base0  <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*Pred.OR.Base0
      } else if (M.type=="ExpXExpY"){
        M.matrix      <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(data$Y0)
        Pred.M.Base0  <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(Theta[4]/2 + Pred.OR.Base0 )
      } else if (M.type=="LogXLinY"){
        M.matrix      <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*(data$Y0)
        Pred.M.Base0  <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*Pred.OR.Base0
      } else if (M.type=="LogXExpY"){
        M.matrix      <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*exp(data$Y0)
        Pred.M.Base0  <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*exp(Theta[4]/2 + Pred.OR.Base0 )
      }  
      
      
      Res.M.Base      <- M.matrix - Pred.M.Base0
      Res.PS.0.Base   <- data$Trt - Pred.PS.Base
      
      mu.Time1.hat <- Pred.OR.Time1.Base + Theta[5]*Suff.OdRa
      
      c( (Theta[2]-Theta[3]) - Theta[1],
         data$Trt*(data$Y1 - Theta[2]),
         ((1-data$Trt)*(exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1)*(data$Y1-mu.Time1.hat) + data$Trt*mu.Time1.hat) - data$Trt*Theta[3],
         as.numeric(data[,jump.pos.d + pos.OR.Time1.Base])*Res.PS.1 )
      
      
    }
  } 
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,         ## OR Time0       6 
                           Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa 6
                           Design.Matrix.Time1,         ## OR Time1       6
                           Design.Matrix.Time0,         ## PS Time0       6
                           Design.Matrix.OdRa*Y0,       ## PS Time0, OdRa 6
                           Design.Matrix.Time1,         ## PS Time1       6
                           Design.Matrix.OdRa))         ## OdRa           6
  
  
  theta <- ParaStart <- c( (mean(Y1[Trt==1])-mean(Y1[Trt==0]))-(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                           (mean(Y1[Trt==1])),
                           (mean(Y1[Trt==0])+(mean(Y0[Trt==1])-mean(Y0[Trt==0])))*mean(Trt),
                           OR.Time1$coefficients)                 # PS Time1
  
  EESolve.Time1 <- m_estimate(estFUN = EE.Time1,
                              data = data,
                              root_control = setup_root_control(start = ParaStart),
                              compute_vcov = FALSE)
  
  # EESolve.Time1 <- m_estimate(estFUN = EE.Time1,
  #                             roots = ParaStart,
  #                             compute_roots=FALSE,
  # data  = data)
  
  CEE1 <- c(coef(EESolve.Time1)[1:3],                            # ATT, ATT1, ATT0
            mean(OR.Time0$residuals^2),                          # Var(Y_{t=0})
            mean(OR.Time1$residuals^2),                          # Var(Y_{t=1})
            OR.Time0$coefficients,                               # OR.Full(t=0)
            OR.Time1$coefficients,                               # OR(t=1)
            PS.Time0$coefficients,                               # PS.Full(t=0)
            coef(EESolve.Time1)[-(1:3)],                         # PS(t=1)
            coef(EESolve.Time0))                                 # OdRa
  
  ###############################################################################
  
  EE <- function(data){
    
    jump.pos.t <- 5
    jump.pos.d <- 3
    pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    pos.PS.Time0.Full <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0+NV.OdRa)
    pos.OR.Time0.Base <- 1:(NV.Time0)
    pos.OR.Time1.Base <- NV.Time0+NV.OdRa+1:(NV.Time1)
    pos.PS.Time0.Base <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0)
    pos.PS.Time1.Base <- (2*NV.Time0+2*NV.OdRa+NV.Time1)+1:(NV.Time1)
    pos.OdRa          <- (2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)
    
    function(theta){
      
      Pred.OR.Time0.Full       <- ( sum( ( theta[jump.pos.t + pos.OR.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Full] ) ) )
      Pred.OR.Time1.Base       <- ( sum( ( theta[jump.pos.t + pos.OR.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time1.Base] ) ) )
      Pred.PS.Time0.Full       <- expit( sum( ( theta[jump.pos.t + pos.PS.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Full] ) ) )
      Pred.PS.Time1.Base.LinLk <- sum( ( theta[jump.pos.t + pos.PS.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time1.Base] ) ) 
      
      Pred.OR.Base0 <- ( sum( ( theta[jump.pos.t + pos.OR.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Base] ) ) )
      Pred.PS.Base <- expit( sum( ( theta[jump.pos.t + pos.PS.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Base] ) ) )
      
      Suff.OdRa          <- ( sum( ( theta[jump.pos.t + pos.OdRa] ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
      Pred.OdRa.Time0     <- LINK(Suff.OdRa*(data$Y0-Pred.OR.Time1.Base ))
      Pred.OdRa.Time1     <- LINK(Suff.OdRa*(data$Y1-Pred.OR.Time1.Base ))
      Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
      
      Res.OR.0 <- data$Y0  - Pred.OR.Time0.Full
      Res.OR.1 <- data$Y1  - Pred.OR.Time1.Base
      Res.PS.0 <- data$Trt - Pred.PS.Time0.Full
      Res.PS.1 <- ( (1-data$Trt)*(1+exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1  )-1 )
      
      if(M.type=="LinXLinY"){
        M.matrix      <- as.numeric( data[,jump.pos.d + pos.OdRa] )*data$Y0
        Pred.M.Base0  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*Pred.OR.Base0
      } else if (M.type=="LinXExpY"){
        M.matrix      <- (as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(data$Y0)
        Pred.M.Base0  <- (as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(theta[4]/2 + Pred.OR.Base0 )
      } else if (M.type=="ExpXLinY"){
        M.matrix      <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*(data$Y0)
        Pred.M.Base0  <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*Pred.OR.Base0
      } else if (M.type=="ExpXExpY"){
        M.matrix      <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(data$Y0)
        Pred.M.Base0  <- exp(as.numeric( data[,jump.pos.d + pos.OdRa] ))*exp(theta[4]/2 + Pred.OR.Base0 )
      } else if (M.type=="LogXLinY"){
        M.matrix      <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*(data$Y0)
        Pred.M.Base0  <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*Pred.OR.Base0
      } else if (M.type=="LogXExpY"){
        M.matrix      <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*exp(data$Y0)
        Pred.M.Base0  <- log(as.numeric( data[,jump.pos.d + pos.OdRa]+1 ))*exp(theta[4]/2 + Pred.OR.Base0 )
      }  
      
      
      Res.M.Base      <- M.matrix - Pred.M.Base0
      Res.PS.0.Base   <- data$Trt - Pred.PS.Base
      
      mu.Time1.hat <- Pred.OR.Time1.Base + theta[5]*Suff.OdRa
      
      c( (theta[2]-theta[3]) - theta[1],
         data$Trt*(data$Y1 - theta[2]),
         ((1-data$Trt)*(exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1)*(data$Y1-mu.Time1.hat) + data$Trt*mu.Time1.hat) - data$Trt*theta[3],
         (Res.OR.0)^2 - theta[4],
         (1-data$Trt)*((Res.OR.1)^2 - theta[5]),
         as.numeric(data[,jump.pos.d + pos.OR.Time0.Full])*Res.OR.0, 
         (1-data$Trt)*as.numeric(data[,jump.pos.d + pos.OR.Time1.Base])*Res.OR.1, 
         as.numeric(data[,jump.pos.d + pos.PS.Time0.Full])*Res.PS.0,
         as.numeric(data[,jump.pos.d + pos.PS.Time1.Base])*Res.PS.1,
         (Res.M.Base)*(Sandwich.OdRa.Time0)*(Res.PS.0.Base)  )
      
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,         ## OR Time0
                           Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa
                           Design.Matrix.Time1,         ## OR Time1
                           Design.Matrix.Time0,         ## PS Time0
                           Design.Matrix.OdRa*Y0,       ## PS Time0, OdRa
                           Design.Matrix.Time1,         ## PS Time1
                           Design.Matrix.OdRa))         ## OdRa
  
  ParaStart <- CEE1
  
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  # EESolve <- m_estimate(estFUN = EE,
  #                       roots = ParaStart,
  #                       compute_roots=FALSE,
  #                       data  = data)
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       "Var0",
                       "Var1",
                       sprintf("OR0_Base_%0.3d",1:NV.Time0),
                       sprintf("OR0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("OR1_Base_%0.3d",1:NV.Time1),
                       sprintf("PS0_Base_%0.3d",1:NV.Time0),
                       sprintf("PS0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("PS1_Base_%0.3d",1:NV.Time1),
                       sprintf("OddRatio_%0.3d",1:NV.OdRa))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full 
  
  return(Result)
  
}


















































UDID.OR.Sensitivity <- function(Y0,Y1,Trt,
                                Design.Matrix.Time0,
                                Design.Matrix.OdRa,
                                Design.Matrix.Time1,
                                Original.Est,
                                Original.Var,
                                Sensitivity,
                                type="Est"){
  
  
  NV.Time0 <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1 <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Trt))
  X.Time1    <- as.matrix(cbind(Design.Matrix.Time1))
  
  data0.OdRa <- data.frame(cbind(Y0,X.Time0.OdRa))
  data1    <- data.frame(cbind(Y1,X.Time1))
  
  colnames(data0.OdRa) <- c("Y0",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  colnames(data1)    <- c("Y1",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  
  LM.Time0.OdRa <- lm(Y0~0+., data=data0.OdRa)
  LM.Time1    <- lm(Y1~0+., data=data1[Trt==0,])
  
  if(type=="Est"){
    Sen.Para <- Sensitivity*Original.Est[3+NV.Time0+(1:NV.OdRa)]
  } else if (type=="SE") {
    Sen.Para <- Sensitivity*(sqrt(diag(Original.Var))[3+NV.Time0+(1:NV.OdRa)])
  } else if (type=="Abs") {
    Sen.Para <- Sensitivity
  }
  
  
  EE <- function(data){
    
    function(theta){
      
      Pred0 <- sum( ( theta[3+1:(NV.Time0+NV.OdRa)] ) * as.numeric( data[,3+1:(NV.Time0+NV.OdRa)] ) )
      Gap   <- sum( ( theta[3+NV.Time0+1:(NV.OdRa)] + Sen.Para ) * as.numeric( data[,3+NV.Time0+NV.OdRa+1:NV.OdRa] ) )
      Res0  <- (data$Y0 - Pred0 )
      Pred1 <- sum( ( theta[3+NV.Time0+NV.OdRa+1:NV.Time1] ) * as.numeric( data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1] ) )
      Res1  <- (data$Y1 - Pred1 )
      
      c( (data$Trt)*( data$Y1-Pred1-theta[3]/theta[2]*(Gap) - theta[1] ), 
         (Res0)^2 - theta[2],
         (1-data$Trt)*((Res1)^2 - theta[3]),
         as.numeric(data[,3+1:(NV.Time0+NV.OdRa)])*Res0,
         (1-data$Trt)*as.numeric(data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1])*Res1)
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,
                           Design.Matrix.OdRa*Trt,
                           Design.Matrix.OdRa,
                           Design.Matrix.Time1))
  
  ParaStart <- Original.Est
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  # EESolve <- m_estimate(estFUN = EE, 
  #                  roots = ParaStart,
  #                  compute_roots=FALSE,
  #                  data  = data)
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "Var0",
                       "Var1",
                       sprintf("OR0_Base_%0.3d",1:NV.Time0),
                       sprintf("OR0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("OR1_Base_%0.3d",1:NV.Time1))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  return(Result)
  
}




UDID.PS.Sensitivity <- function(Y0,Y1,Trt,
                                Design.Matrix.Time0,
                                Design.Matrix.OdRa,
                                Design.Matrix.Time1,
                                Original.Est,
                                Original.Var,
                                Sensitivity,
                                type="Est"){
  
  
  
  NV.Time0 <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1 <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Y0))
  data0.OdRa   <- data.frame(cbind(Trt,X.Time0.OdRa))
  colnames(data0.OdRa) <- c("Trt",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  pi.Time0.OdRa <- glm(Trt~0+., data=data0.OdRa, family="binomial")
  
  X.Time1    <- as.matrix(cbind(Design.Matrix.Time1))
  data.Time1      <- data.frame(cbind(Trt,
                                      (Design.Matrix.OdRa*Y1)%*%pi.Time0.OdRa$coefficients[NV.Time0+1:NV.OdRa],
                                      X.Time1))
  colnames(data.Time1)    <- c("Trt","Base",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  pi.Time1 <- glm(Trt~0+., data=data.Time1,family="binomial")
  
  if(type=="Est"){
    Sen.Para <- Sensitivity*Original.Est[3+NV.Time0+1:NV.OdRa]
  } else if (type=="SE") {
    Sen.Para <- Sensitivity*(sqrt(diag(Original.Var))[3+NV.Time0+1:NV.OdRa])
  } else if (type=="Abs") {
    Sen.Para <- Sensitivity
  }
  
  EE <- function(data){
    
    function(theta){
      
      Pred0 <- expit( sum( ( theta[3+1:(NV.Time0+NV.OdRa)] ) * as.numeric( data[,3+1:(NV.Time0+NV.OdRa)] ) ) )
      Res0 <- data$Trt - Pred0
      Gap   <- sum( ( theta[3+NV.Time0+1:(NV.OdRa)] + Sen.Para ) * as.numeric( data[,3+NV.Time0+NV.OdRa+1:NV.OdRa] ) )
      OR0 <- exp( data$Y1*(Gap) )
      Res1 <- ( (1-data$Trt)*(1+exp( sum( ( theta[3+NV.Time0+NV.OdRa+1:NV.Time1] ) * 
                                            as.numeric( data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1] ) ))*OR0)-1 )
      IPW <- exp( sum( ( theta[3+NV.Time0+NV.OdRa+1:NV.Time1] ) * 
                         as.numeric( data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1] ) ))*OR0
      
      c( theta[1]-(theta[2]-theta[3]),
         (data$Trt)*(data$Y1 - theta[2]), 
         (1-data$Trt)*(data$Y1)*IPW - (1-data$Trt)*IPW*theta[3],
         as.numeric(data[,3+1:(NV.Time0+NV.OdRa)])*Res0 , 
         as.numeric(data[,3+NV.Time0+2*NV.OdRa+1:NV.Time1])*Res1  )
      
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,
                           Design.Matrix.OdRa*Y0,
                           Design.Matrix.OdRa,
                           Design.Matrix.Time1))
  
  ParaStart <- Original.Est
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  # EESolve <- m_estimate(estFUN = EE, 
  #                  roots = ParaStart,
  #                  compute_roots=FALSE,
  #                  data  = data)
  
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("PS0_Base_%0.3d",1:NV.Time0),
                       sprintf("PS0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("PS1_Base_%0.3d",1:NV.Time1))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
}

UDID.DR.Sensitivity <- function(Y0,Y1,Trt,
                                Design.Matrix.Time0,
                                Design.Matrix.OdRa,
                                Design.Matrix.Time1,
                                Original.Est,
                                Original.Var,
                                Sensitivity,
                                type="Est"){
  
  ## Moment Equation for OdRa
  
  NV.Time0   <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1   <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OR <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Trt))
  X.Time0.PS <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Y0))
  X.Time1.OR <- as.matrix(cbind(Design.Matrix.Time1))
  X.Time1.PS <- as.matrix(cbind(Design.Matrix.Time1))
  
  data.Time0.OR   <- data.frame(cbind(Y0, X.Time0.OR))
  data.Time1.OR   <- data.frame(cbind(Y1, X.Time1.OR))
  data.Time0.PS   <- data.frame(cbind(Trt,X.Time0.PS))
  
  colnames(data.Time0.OR) <- c("Y0",
                               sprintf("OR_%0.3d",1:dim(X.Time0.OR)[2]))
  colnames(data.Time1.OR) <- c("Y1",
                               sprintf("OR_%0.3d",1:dim(X.Time1.OR)[2]))
  colnames(data.Time0.PS) <- c("Trt",
                               sprintf("PS_%0.3d",1:dim(X.Time0.OR)[2]))
  
  OR.Time0  <- lm(Y0~0+.,  data=data.Time0.OR)
  OR.Time1  <- lm(Y1~0+.,  data=data.Time1.OR[Trt==0,])
  PS.Time0  <- glm(Trt~0+., data=data.Time0.PS,family="binomial")
  
  data.Time1.PS      <- data.frame(cbind(Trt,
                                         (Design.Matrix.OdRa*Y1)%*%PS.Time0$coefficients[NV.Time0+1:NV.OdRa],
                                         X.Time1.PS))
  colnames(data.Time1.PS) <- c("Trt","Base",
                               sprintf("PS_%0.3d",1:dim(X.Time1.PS)[2]))
  PS.Time1  <- glm(Trt~0+., data=data.Time1.PS,family="binomial")
  
  if(type=="Est"){
    Sen.Para <- Sensitivity*Original.Est[5+(2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)]
  } else if (type=="SE") {
    Sen.Para <- Sensitivity*(sqrt(diag(Original.Var))[5+(2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)])
  } else if (type=="Abs") {
    Sen.Para <- Sensitivity
  }
  
  
  ###############################################################################
  
  # EE.Time1 <- function(data){
  # 
  #   jump.pos.t <- 5
  #   jump.pos.d <- 3
  #   pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
  #   pos.PS.Time0.Full <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0+NV.OdRa)
  #   pos.OR.Time0.Base <- 1:(NV.Time0)
  #   pos.OR.Time1.Base <- NV.Time0+NV.OdRa+1:(NV.Time1)
  #   pos.PS.Time0.Base <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0)
  #   pos.PS.Time1.Base <- (2*NV.Time0+2*NV.OdRa+NV.Time1)+1:(NV.Time1)
  #   pos.OdRa          <- (2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)
  # 
  #   function(theta){
  # 
  #     Theta <- c(theta[c(1:3)],                                                           # ATT, ATT1, ATT0
  #                Original.Est[c(4,5,jump.pos.t+1:(NV.Time0+NV.OdRa))],                    # Var(Y_{t=0}), Var(Y_{t=1}), OR.Full(t=0)
  #                OR.Time1$coefficients,                                                   # OR(t=1)
  #                Original.Est[jump.pos.t+(NV.Time0+NV.OdRa)+1:(NV.Time0+NV.OdRa)],        # PS.Full(t=0)
  #                theta[-c(1:3)],                                                          # PS(t=1)
  #                Original.Est[jump.pos.t+(2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)])               # OdRa
  # 
  #     Pred.OR.Time0.Full       <- ( sum( ( Theta[jump.pos.t + pos.OR.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Full] ) ) )
  #     Pred.OR.Time1.Base       <- ( sum( ( Theta[jump.pos.t + pos.OR.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time1.Base] ) ) )
  #     Pred.PS.Time0.Full       <- expit( sum( ( Theta[jump.pos.t + pos.PS.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Full] ) ) )
  #     Pred.PS.Time1.Base.LinLk <- sum( ( Theta[jump.pos.t + pos.PS.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time1.Base] ) )
  # 
  #     Pred.OR.Base <- ( sum( ( Theta[jump.pos.t + pos.OR.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Base] ) ) )
  #     Pred.PS.Base <- expit( sum( ( Theta[jump.pos.t + pos.PS.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Base] ) ) )
  # 
  #     Suff.OdRa0          <- ( sum( ( Theta[jump.pos.t + pos.OdRa]            ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
  #     Suff.OdRa1          <- ( sum( ( Theta[jump.pos.t + pos.OdRa] + Sen.Para ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
  #     
  #     Pred.OdRa.Time0     <- LINK(Suff.OdRa0*(data$Y0-Pred.OR.Time1.Base ))
  #     Pred.OdRa.Time1     <- LINK(Suff.OdRa1*(data$Y1-Pred.OR.Time1.Base ))
  #     Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
  # 
  #     Res.OR.0 <- data$Y0  - Pred.OR.Time0.Full
  #     Res.OR.1 <- data$Y1  - Pred.OR.Time1.Base
  #     Res.PS.0 <- data$Trt - Pred.PS.Time0.Full
  #     Res.PS.1 <- ( (1-data$Trt)*(1+exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1  )-1 )
  # 
  #     M.matrix     <- as.numeric( data[,jump.pos.d + pos.OdRa] )*data$Y0
  #     Pred.M.Base  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*Pred.OR.Base
  # 
  #     Res.M.Base      <- M.matrix - Pred.M.Base
  #     Res.PS.0.Base   <- data$Trt - Pred.PS.Base
  # 
  #     mu.Time1.hat <- Pred.OR.Time1.Base + Theta[5]*Suff.OdRa1
  # 
  #     c( (Theta[2]-Theta[3]) - Theta[1],
  #        data$Trt*(data$Y1 - Theta[2]),
  #        ((1-data$Trt)*(exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1)*(data$Y1-mu.Time1.hat) + data$Trt*mu.Time1.hat) - data$Trt*Theta[3],
  #        as.numeric(data[,jump.pos.d + pos.OR.Time1.Base])*Res.PS.1 )
  # 
  # 
  #   }
  # }
  # 
  # data <- data.frame(cbind(Y0,Y1,Trt,
  #                          Design.Matrix.Time0,         ## OR Time0       6
  #                          Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa 6
  #                          Design.Matrix.Time1,         ## OR Time1       6
  #                          Design.Matrix.Time0,         ## PS Time0       6
  #                          Design.Matrix.OdRa*Y0,       ## PS Time0, OdRa 6
  #                          Design.Matrix.Time1,         ## PS Time1       6
  #                          Design.Matrix.OdRa))         ## OdRa           6
  # 
  # 
  # theta <- ParaStart <- c( (mean(Y1[Trt==1])-mean(Y1[Trt==0]))-(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
  #                          (mean(Y1[Trt==1])),
  #                          (mean(Y1[Trt==0])+(mean(Y0[Trt==1])-mean(Y0[Trt==0])))*mean(Trt),
  #                          PS.Time1$coefficients[-1])                 # PS Time1
  # 
  # EESolve.Time1 <- m_estimate(estFUN = EE.Time1,
  #                             data = data,
  #                             root_control = setup_root_control(start = ParaStart),
  #                             compute_vcov = FALSE)
  # 
  # # EESolve.Time1 <- m_estimate(estFUN = EE.Time1,
  # #                             roots = ParaStart,
  # #                             compute_roots=FALSE,
  # # data  = data)
  # 
  # CEE1 <- c(coef(EESolve.Time1)[1:3],                            # ATT, ATT1, ATT0
  #           mean(OR.Time0$residuals^2),                          # Var(Y_{t=0})
  #           mean(OR.Time1$residuals^2),                          # Var(Y_{t=1})
  #           OR.Time0$coefficients,                               # OR.Full(t=0)
  #           OR.Time1$coefficients,                               # OR(t=1)
  #           PS.Time0$coefficients,                               # PS.Full(t=0)
  #           coef(EESolve.Time1)[-(1:3)],                         # PS(t=1)
  #           coef(EESolve.Time0))                                 # OdRa
  
  
  EE <- function(data){
    
    jump.pos.t <- 5
    jump.pos.d <- 3
    pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    pos.PS.Time0.Full <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0+NV.OdRa)
    pos.OR.Time0.Base <- 1:(NV.Time0)
    pos.OR.Time1.Base <- NV.Time0+NV.OdRa+1:(NV.Time1)
    pos.PS.Time0.Base <- (NV.Time0+NV.OdRa+NV.Time1)+1:(NV.Time0)
    pos.PS.Time1.Base <- (2*NV.Time0+2*NV.OdRa+NV.Time1)+1:(NV.Time1)
    pos.OdRa          <- (2*NV.Time0+2*NV.OdRa+NV.Time1*2)+1:(NV.OdRa)
    
    function(theta){
      
      Pred.OR.Time0.Full       <- ( sum( ( theta[jump.pos.t + pos.OR.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Full] ) ) )
      Pred.OR.Time1.Base       <- ( sum( ( theta[jump.pos.t + pos.OR.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time1.Base] ) ) )
      Pred.PS.Time0.Full       <- expit( sum( ( theta[jump.pos.t + pos.PS.Time0.Full] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Full] ) ) )
      Pred.PS.Time1.Base.LinLk <- sum( ( theta[jump.pos.t + pos.PS.Time1.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time1.Base] ) ) 
      
      Pred.OR.Base <- ( sum( ( theta[jump.pos.t + pos.OR.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.OR.Time0.Base] ) ) )
      Pred.PS.Base <- expit( sum( ( theta[jump.pos.t + pos.PS.Time0.Base] ) * as.numeric( data[,jump.pos.d + pos.PS.Time0.Base] ) ) )
      
      Suff.OdRa0          <- ( sum( ( theta[jump.pos.t + pos.OdRa]            ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
      Suff.OdRa1          <- ( sum( ( theta[jump.pos.t + pos.OdRa] + Sen.Para ) * as.numeric( data[,jump.pos.d + pos.OdRa] )))
      Pred.OdRa.Time0     <- LINK(Suff.OdRa0*(data$Y0-Pred.OR.Time1.Base ))
      Pred.OdRa.Time1     <- LINK(Suff.OdRa1*(data$Y1-Pred.OR.Time1.Base ))
      Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
      
      Res.OR.0 <- data$Y0  - Pred.OR.Time0.Full
      Res.OR.1 <- data$Y1  - Pred.OR.Time1.Base
      Res.PS.0 <- data$Trt - Pred.PS.Time0.Full
      Res.PS.1 <- ( (1-data$Trt)*(1+exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1  )-1 )
      
      M.matrix     <- as.numeric( data[,jump.pos.d + pos.OdRa] )*data$Y0
      Pred.M.Base  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*Pred.OR.Base
      # M.matrix     <- as.numeric( data[,jump.pos.d + pos.OdRa] )*exp(data$Y0/sqrt(Theta[5]))
      # Pred.M.Base  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*exp( 1/2 + Pred.OR.Base/sqrt(Theta[5]) )
      # M.matrix     <- as.numeric( data[,jump.pos.d + pos.OdRa] )*exp(data$Y0)
      # Pred.M.Base  <- as.numeric( data[,jump.pos.d + pos.OdRa] )*exp( Theta[1]/2 + Pred.OR.Base )
      
      Res.M.Base      <- M.matrix - Pred.M.Base
      Res.PS.0.Base   <- data$Trt - Pred.PS.Base
      
      mu.Time1.hat <- Pred.OR.Time1.Base + theta[5]*Suff.OdRa1
      
      c( (theta[2]-theta[3]) - theta[1],
         data$Trt*(data$Y1 - theta[2]),
         ((1-data$Trt)*(exp( Pred.PS.Time1.Base.LinLk )*Pred.OdRa.Time1)*(data$Y1-mu.Time1.hat) + data$Trt*mu.Time1.hat) - data$Trt*theta[3],
         (Res.OR.0)^2 - theta[4],
         (1-data$Trt)*((Res.OR.1)^2 - theta[5]),
         as.numeric(data[,jump.pos.d + pos.OR.Time0.Full])*Res.OR.0, 
         (1-data$Trt)*as.numeric(data[,jump.pos.d + pos.OR.Time1.Base])*Res.OR.1, 
         as.numeric(data[,jump.pos.d + pos.PS.Time0.Full])*Res.PS.0,
         as.numeric(data[,jump.pos.d + pos.PS.Time1.Base])*Res.PS.1,
         (Res.M.Base)*(Sandwich.OdRa.Time0)*(Res.PS.0.Base)  )
      
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,         ## OR Time0
                           Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa
                           Design.Matrix.Time1,         ## OR Time1
                           Design.Matrix.Time0,         ## PS Time0
                           Design.Matrix.OdRa*Y0,       ## PS Time0, OdRa
                           Design.Matrix.Time1,         ## PS Time1
                           Design.Matrix.OdRa))         ## OdRa
  
  ParaStart <- theta <- Original.Est
  
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  # EESolve <- m_estimate(estFUN = EE,
  #                       roots = ParaStart,
  #                       compute_roots=FALSE,
  #                       data  = data)
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       "Var0",
                       "Var1",
                       sprintf("OR0_Base_%0.3d",1:NV.Time0),
                       sprintf("OR0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("OR1_Base_%0.3d",1:NV.Time1),
                       sprintf("PS0_Base_%0.3d",1:NV.Time0),
                       sprintf("PS0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("PS1_Base_%0.3d",1:NV.Time1),
                       sprintf("OddRatio_%0.3d",1:NV.OdRa))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full 
  
  return(Result)
  
}



























library(nnet)




##################################################

UDID.OR.Disc <- function(Y0,Y1,Trt,
                         Design.Matrix.Time0,
                         Design.Matrix.OdRa,
                         Design.Matrix.Time1,
                         num.bin=5){
  
  
  NV.Time0 <- dim(Design.Matrix.Time0)[2]
  NV.OdRa    <- dim(Design.Matrix.OdRa)[2]
  NV.Time1 <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Trt))
  X.Time1    <- as.matrix(cbind(Design.Matrix.Time1))
  
  Bin.Cut    <- quantile(Y0,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
  Disc.Y0    <- apply(matrix(Y0,length(Y0),1),1,function(v){sum(v>Bin.Cut)})
  Disc.Y1    <- apply(matrix(Y1,length(Y1),1),1,function(v){sum(v>Bin.Cut)})
  
  Disc.Y0.Mat <- Disc.Y1.Mat <- matrix(0,length(Disc.Y0),num.bin)
  for(jj in 1:num.bin){
    Disc.Y0.Mat[Disc.Y0==jj-1,jj] <- 1
    Disc.Y1.Mat[Disc.Y1==jj-1,jj] <- 1
  }
  
  data0.OdRa <- data.frame(cbind(as.factor(Disc.Y0),X.Time0.OdRa))
  data1      <- data.frame(cbind(as.factor(Disc.Y1),X.Time1))
  
  colnames(data0.OdRa) <- c("Y0" ,sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  colnames(data1)      <- c("Y1" ,sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  
  GLM.Time0.OdRa <- multinom(Y0~0+., data=data0.OdRa)
  GLM.Time1      <- multinom(Y1~0+., data=data1[Trt==0,])
  LM.Time1       <- lm(Y1[Trt==0]~0+X.Time1[Trt==0,])
  
  Prob.Ctr <- cbind(1,exp( cbind(Design.Matrix.Time0,Design.Matrix.OdRa*0)%*%t(coef(GLM.Time0.OdRa)) ))
  Prob.Trt <- cbind(1,exp( cbind(Design.Matrix.Time0,Design.Matrix.OdRa*1)%*%t(coef(GLM.Time0.OdRa)) ))
  OdRa0.Vec <- (Prob.Trt/Prob.Ctr)
  
  Var1.Linear <- sum(LM.Time1$residuals^2)/sum(1-Trt)
  
  Pred1.Linear <- Design.Matrix.Time1%*%as.numeric( LM.Time1$coefficients )
  
  mu1 <- rep(0,length(Y0))
  
  for(ii in 1:length(Y0)){
    
    QN1 <- pnorm(c(Bin.Cut,Inf),Pred1.Linear[ii],Var1.Linear^(0.5))
    QN1.R <- (c(QN1[1],diff(QN1)))
    QN2 <- -0.5*exp(-(c(Bin.Cut,Inf)-Pred1.Linear[ii])^2/2/Var1.Linear)
    QN2.R <- (c(QN2[1],diff(QN2)))
    
    Numer <- Pred1.Linear[ii]*OdRa0.Vec[ii,]*QN1.R + sqrt(2*Var1.Linear/pi)*OdRa0.Vec[ii,]*QN2.R
    Denom <- OdRa0.Vec[ii,]*QN1.R
    
    mu1[ii] <- sum(Numer)/sum(Denom)
    
  }
  
  
  
  EE <- function(data){
    
    t.pos <- 4
    d.pos <- 2+2*num.bin+1
    
    function(theta){
      
      Theta.OR.Time0 <- (matrix( theta[t.pos+1:((num.bin-1)*(NV.Time0+NV.OdRa))], (NV.Time0+NV.OdRa),num.bin-1 ))
      Pred0.Exp      <- exp(c(0,apply( Theta.OR.Time0 * as.numeric( data[,d.pos+1:(NV.Time0+NV.OdRa)] ),2,sum )))
      Pred0.Prob     <- Pred0.Exp/(sum(Pred0.Exp ))
      
      Design.Control <- c( as.numeric(data[,d.pos+1:(NV.Time0)]),
                           as.numeric(data[,d.pos+(NV.Time0+NV.OdRa+NV.Time1)+1:(NV.OdRa)])*0 )
      
      Design.Trt     <- c( as.numeric(data[,d.pos+1:(NV.Time0)]),
                           as.numeric(data[,d.pos+(NV.Time0+NV.OdRa+NV.Time1)+1:(NV.OdRa)])*1 )
      
      Pred0.Exp.Control <- exp(c(0,apply( Theta.OR.Time0 * Design.Control,2,sum )))
      Pred0.Exp.Trt     <- exp(c(0,apply( Theta.OR.Time0 * Design.Trt    ,2,sum )))
      
      Pred0.Prob.Control <- Pred0.Exp.Control/sum(Pred0.Exp.Control)
      Pred0.Prob.Trt     <- Pred0.Exp.Trt/sum(Pred0.Exp.Trt)
      OdRa0.Vec <- (Pred0.Prob.Trt/Pred0.Prob.Trt[1])/(Pred0.Prob.Control/Pred0.Prob.Control[1])
      
      
      Theta.OR.Time1 <- (matrix( theta[t.pos+((num.bin-1)*(NV.Time0+NV.OdRa))+1:((num.bin-1)*(NV.Time1))], 
                                 (NV.Time1),num.bin-1 ))
      Pred1.Exp  <- exp(c(0,apply( Theta.OR.Time1 * as.numeric( data[,d.pos+(NV.Time0+NV.OdRa)+1:(NV.Time1)] ),2,sum )))
      Pred1.Prob <- Pred1.Exp/(sum(Pred1.Exp ))
      
      Res0 <- Res1 <- rep(0,num.bin-1)
      for(jj in 1:(num.bin-1)){
        Res0[jj] <- data[,2+jj+1]-Pred0.Prob[jj+1]
        Res1[jj] <- (data[,2+num.bin+jj+1]-Pred1.Prob[jj+1])*(1-data$Trt)
      }
      
      Pred1.Linear <- sum( ( theta[t.pos+((num.bin-1)*(NV.Time0+NV.OdRa+NV.Time1))+(1:NV.Time1)] ) *
                             as.numeric( data[,d.pos+(NV.Time0+NV.OdRa)+1:(NV.Time1)] ) )
      Res1.Linear  <- (data$Y1 - Pred1.Linear )
      
      QN1 <- pnorm(c(Bin.Cut,Inf),Pred1.Linear,theta[4]^(0.5))
      QN1.R <- (c(QN1[1],diff(QN1)))
      QN2 <- -0.5*exp(-(c(Bin.Cut,Inf)-Pred1.Linear)^2/2/theta[4])
      QN2.R <- (c(QN2[1],diff(QN2)))
      
      Numer <- Pred1.Linear*OdRa0.Vec*QN1.R + sqrt(2*theta[4]/pi)*OdRa0.Vec*QN2.R
      Denom <- OdRa0.Vec*QN1.R
      
      mu1 <- sum(Numer)/sum(Denom)
      
      c( theta[2]-theta[3]-theta[1],
         data$Trt*(data$Y1-theta[2]),
         data$Trt*(mu1-theta[3]),
         (1-data$Trt)*(Res1.Linear^2-theta[4]),
         rep(Res0,each=(NV.Time0+NV.OdRa))*rep(as.numeric( data[,d.pos+1:(NV.Time0+NV.OdRa)] ),num.bin-1), # OR parameter
         rep(Res1,each=(NV.Time1))*rep(as.numeric( data[,d.pos+(NV.Time0+NV.OdRa)+1:(NV.Time1)] ),num.bin-1),
         (1-data$Trt)*Res1.Linear*as.numeric( data[,d.pos+(NV.Time0+NV.OdRa)+1:(NV.Time1)] ))
      
    }
  }
  
  data <- data.frame(cbind(Y0,Y1,Disc.Y0.Mat,Disc.Y1.Mat,Trt,
                           Design.Matrix.Time0,
                           Design.Matrix.OdRa*Trt,
                           Design.Matrix.Time1,
                           Design.Matrix.OdRa))
  
  ParaStart <- theta <- c( mean(Y1[Trt==1])-mean(mu1[Trt==1]),
                           mean(Y1[Trt==1]),
                           mean(mu1[Trt==1]),
                           sum(LM.Time1$residuals^2)/sum(1-Trt),
                           as.vector(t(coef(GLM.Time0.OdRa))),
                           as.vector(t(coef(GLM.Time1))),
                           as.numeric( LM.Time1$coefficients ))
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  # EESolve <- m_estimate(estFUN = EE,
  #                  roots = ParaStart,
  #                  compute_roots=FALSE,
  #                  data  = data)
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       "Var1",
                       sprintf("OR0_Base_%0.3d",1:NV.Time0),
                       sprintf("OR0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("OR1_Base_%0.3d",1:NV.Time1),
                       sprintf("OdRa_%0.3d",1:NV.Time1))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- names(CEE.Full)
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  return(Result)
  
}














##################################################

UDID.PS.Disc <- function(Y0,Y1,Trt,
                         Design.Matrix.Time0,
                         Design.Matrix.OdRa.Y0.Full,
                         Design.Matrix.Time1,
                         Design.Matrix.OdRa.Y1.Full,
                         num.bin=5){
  
  Design.Matrix.OdRa.Y0 <- Design.Matrix.OdRa.Y0.Full[,-1]
  Design.Matrix.OdRa.Y1 <- Design.Matrix.OdRa.Y1.Full[,-1]
  
  NV.Time0 <- dim(Design.Matrix.Time0)[2]
  NV.OdRa  <- dim(Design.Matrix.OdRa.Y0)[2]
  NV.Time1 <- dim(Design.Matrix.Time1)[2]
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,
                                  Design.Matrix.OdRa.Y0))
  data0.OdRa   <- data.frame(cbind(Trt,X.Time0.OdRa))
  colnames(data0.OdRa) <- c("Trt",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  
  pi.Time0.OdRa <- glm(Trt~0+., data=data0.OdRa, family="binomial")
  
  X.Time1         <- as.matrix(cbind(Design.Matrix.Time1))
  data.Time1      <- data.frame(cbind(Trt,
                                      (Design.Matrix.OdRa.Y1)%*%pi.Time0.OdRa$coefficients[NV.Time0+1:NV.OdRa],
                                      X.Time1))
  colnames(data.Time1)    <- c("Trt","Base",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  pi.Time1 <- glm(Trt~0+., data=data.Time1,family="binomial")
  
  EE <- function(data){
    
    t.pos <- 3
    d.pos <- 3
    
    function(theta){
      
      Pred0 <- expit( sum( ( theta[t.pos+1:(NV.Time0+NV.OdRa)] ) * as.numeric( data[,d.pos+1:(NV.Time0+NV.OdRa)] ) ) )
      Res0 <- data$Trt - Pred0
      Gap1   <- sum( ( theta[t.pos+NV.Time0+1:(NV.OdRa)] ) * as.numeric( data[,d.pos+NV.Time0+NV.OdRa+NV.Time1+1:NV.OdRa] ) )
      OR1 <- exp( Gap1 )
      Res1 <- ( (1-data$Trt)*(1+exp( sum( ( theta[t.pos+NV.Time0+NV.OdRa+1:NV.Time1] ) * 
                                            as.numeric( data[,d.pos+NV.Time0+NV.OdRa+1:NV.Time1] ) ))*OR1)-1 )
      IPW <- exp( sum( ( theta[t.pos+NV.Time0+NV.OdRa+1:NV.Time1] ) * 
                         as.numeric( data[,d.pos+NV.Time0+NV.OdRa+1:NV.Time1] ) ))*OR1
      
        c( theta[1]-(theta[2]-theta[3]),
         (data$Trt)*(data$Y1 - theta[2]), 
         (1-data$Trt)*(data$Y1)*IPW - (1-data$Trt)*IPW*theta[3],
         as.numeric(data[,d.pos+1:(NV.Time0+NV.OdRa)])*Res0 , 
         as.numeric(data[,d.pos+NV.Time0+2*NV.OdRa+1:NV.Time1])*Res1  )
      
    }
  }
  
  
  data <- data.frame(cbind(Y0,Y1,Trt,
                           Design.Matrix.Time0,
                           Design.Matrix.OdRa.Y0,
                           Design.Matrix.Time1,
                           Design.Matrix.OdRa.Y1))
  
  ParaStart <- theta <- c( mean(Y1[Trt==1])-mean(Y1[Trt==0])+
                             -(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                           mean(Y1[Trt==1]),
                           mean(Y1[Trt==0])+(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                           pi.Time0.OdRa$coefficients,
                           pi.Time1$coefficients[-1] )
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("PS0_Base_%0.3d",1:NV.Time0),
                       sprintf("PS0_OdRa_%0.3d",1:NV.OdRa),
                       sprintf("PS1_Base_%0.3d",1:NV.Time1))
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
}



























UDID.DR.Disc <- function(Y0,Y1,Trt,
                         Design.Matrix.Time0,
                         Design.Matrix.OdRa,
                         Design.Matrix.OdRa.Y0.Full,
                         Design.Matrix.Time1,
                         Design.Matrix.OdRa.Y1.Full,
                         num.bin=5){
  
  Design.Matrix.OdRa.Y0 <- Design.Matrix.OdRa.Y0.Full[,-1]
  Design.Matrix.OdRa.Y1 <- Design.Matrix.OdRa.Y1.Full[,-1]
  ## Moment Equation for OdRa
  
  NV.Time0    <- dim(Design.Matrix.Time0)[2]
  NV.OdRa     <- dim(Design.Matrix.OdRa)[2]
  NV.OdRa.Y   <- dim(Design.Matrix.OdRa.Y0)[2]
  NV.Time1    <- dim(Design.Matrix.Time1)[2]
  
  ## Moment Equation for OR
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,Design.Matrix.OdRa*Trt))
  X.Time1    <- as.matrix(cbind(Design.Matrix.Time1))
  
  data0.OdRa <- data.frame(cbind(Y0,X.Time0.OdRa))
  data1    <- data.frame(cbind(Y1,X.Time1))
  
  colnames(data0.OdRa) <- c("Y0",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  colnames(data1)    <- c("Y1",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  
  LM.Time0.OdRa <- lm(Y0~0+., data=data0.OdRa)
  LM.Time1      <- lm(Y1~0+., data=data1[Trt==0,])
  
  
  Bin.Cut    <- quantile(Y0,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
  Disc.Y0    <- apply(matrix(Y0,length(Y0),1),1,function(v){sum(v>Bin.Cut)})
  Disc.Y1    <- apply(matrix(Y1,length(Y1),1),1,function(v){sum(v>Bin.Cut)})
  
  Disc.Y0.Mat <- Disc.Y1.Mat <- matrix(0,length(Disc.Y0),num.bin)
  for(jj in 1:num.bin){
    Disc.Y0.Mat[Disc.Y0==jj-1,jj] <- 1
    Disc.Y1.Mat[Disc.Y1==jj-1,jj] <- 1
  }
  
  data0.OdRa <- data.frame(cbind(as.factor(Disc.Y0),X.Time0.OdRa))
  data1      <- data.frame(cbind(as.factor(Disc.Y1),X.Time1))
  
  colnames(data0.OdRa) <- c("Y0" ,sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  colnames(data1)      <- c("Y1" ,sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  
  
  GLM.Time0.OdRa <- multinom(Y0~0+., data=data0.OdRa)
  GLM.Time1      <- multinom(Y1~0+., data=data1[Trt==0,])
  
  Prob.Ctr <- cbind(1,exp( cbind(Design.Matrix.Time0,Design.Matrix.OdRa*0)%*%t(coef(GLM.Time0.OdRa)) ))
  Prob.Trt <- cbind(1,exp( cbind(Design.Matrix.Time0,Design.Matrix.OdRa*1)%*%t(coef(GLM.Time0.OdRa)) ))
  OdRa0.Vec <- (Prob.Trt/Prob.Ctr)
  
  Var1.Linear <- sum(LM.Time1$residuals^2)/sum(1-Trt)
  
  Pred1.Linear <- Design.Matrix.Time1%*%as.numeric( LM.Time1$coefficients )
  
  mu1 <- rep(0,length(Y0))
  
  for(ii in 1:length(Y0)){
    
    QN1 <- pnorm(c(Bin.Cut,Inf),Pred1.Linear[ii],Var1.Linear^(0.5))
    QN1.R <- (c(QN1[1],diff(QN1)))
    QN2 <- -0.5*exp(-(c(Bin.Cut,Inf)-Pred1.Linear[ii])^2/2/Var1.Linear)
    QN2.R <- (c(QN2[1],diff(QN2)))
    
    Numer <- Pred1.Linear[ii]*OdRa0.Vec[ii,]*QN1.R + sqrt(2*Var1.Linear/pi)*OdRa0.Vec[ii,]*QN2.R
    Denom <- OdRa0.Vec[ii,]*QN1.R
    
    mu1[ii] <- sum(Numer)/sum(Denom)
    
  }
  
  ## Moment Equation for PS
  
  X.Time0.OdRa <- as.matrix(cbind(Design.Matrix.Time0,
                                  Design.Matrix.OdRa.Y0))
  data0.OdRa   <- data.frame(cbind(Trt,X.Time0.OdRa))
  colnames(data0.OdRa) <- c("Trt",sprintf("X_%0.3d",1:dim(X.Time0.OdRa)[2]))
  
  pi.Time0.OdRa <- glm(Trt~0+., data=data0.OdRa, family="binomial")
  
  X.Time1         <- as.matrix(cbind(Design.Matrix.Time1))
  data.Time1      <- data.frame(cbind(Trt,
                                      (Design.Matrix.OdRa.Y1)%*%pi.Time0.OdRa$coefficients[NV.Time0+1:NV.OdRa.Y],
                                      X.Time1))
  colnames(data.Time1)    <- c("Trt","Base",sprintf("X_%0.3d",1:dim(X.Time1)[2]))
  pi.Time1 <- glm(Trt~0+., data=data.Time1,family="binomial")
  
  
  # X_001      X_002      X_003      X_004      X_005      X_001 
  # -0.8211651  2.0094501  2.8306152 -4.1972019  0.9468275  2.6141973  3.7101870  5.4076056 -4.8404992 
  
  #########################################################################
  
  
  
  
  data <- data.frame(cbind(Y0,Y1,
                           Disc.Y0.Mat,Disc.Y1.Mat,
                           Trt,
                           Design.Matrix.Time0,         ## OR Time0
                           Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa
                           Design.Matrix.Time0,         ## PS Time0
                           Design.Matrix.OdRa.Y0,       ## PS Time0, OdRa.Y0
                           Design.Matrix.OdRa.Y0))      ## OdRa.Y0
  colnames(data) <- c("Y0","Y1",
                      sprintf("DiscY0_%0.3d",1:num.bin),
                      sprintf("DiscY1_%0.3d",1:num.bin),
                      "Trt",
                      sprintf("OR_Base%0.3d",1:NV.Time0),
                      sprintf("OR_OdRa%0.3d",1:NV.OdRa),
                      sprintf("PS_Base%0.3d",1:NV.Time0),
                      sprintf("PS_OdRa%0.3d",1:NV.OdRa.Y),
                      sprintf("OREC_OdRaY%0.3d",1:NV.OdRa.Y))
  
  theta <- ParaStart <- c( pi.Time0.OdRa$coefficients[-(1:NV.Time0)] )
  
  
  EE.Time0 <- function(data){
    
    t.pos <- 0
    d.pos <- 3+2*(num.bin)
    
    t.pos.OR.Mult.Time0.Full <- 1:((num.bin-1)*(NV.Time0+NV.OdRa))
    d.pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    d.pos.OR.Time0.Base <- 1:(NV.Time0)
    d.pos.OR.Time0.OdRa <- NV.Time0+1:(NV.OdRa)
    
    t.pos.PS.Time0.Full <- ((num.bin-1)*(NV.Time0+NV.OdRa))+1:(NV.Time0+NV.OdRa.Y)
    t.pos.PS.Time0.Base <- ((num.bin-1)*(NV.Time0+NV.OdRa))+1:(NV.Time0)
    d.pos.PS.Time0.Full <- (NV.Time0+NV.OdRa)+1:(NV.Time0+NV.OdRa.Y)
    d.pos.PS.Time0.Base <- (NV.Time0+NV.OdRa)+1:(NV.Time0)
    
    t.pos.OREC <- ((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0+NV.OdRa.Y)+1:NV.OdRa.Y
    d.pos.OREC <- (NV.Time0+NV.OdRa)+(NV.Time0+NV.OdRa.Y)+1:NV.OdRa.Y
    
    function(theta){
      
      Theta <- c( as.numeric(t(coef(GLM.Time0.OdRa))),
                  pi.Time0.OdRa$coefficients,
                  theta)
      
      names(Theta) <- c(sprintf("OR_Mult%0.3d",1:length(t.pos.OR.Mult.Time0.Full)),
                        sprintf("PS_Base%0.3d",1:NV.Time0),
                        sprintf("PS_OdRa%0.3d",1:NV.OdRa.Y),
                        sprintf("OREC_OdRa%0.3d",1:NV.OdRa.Y))
      
      
      Theta.OR.Time0 <- (matrix( Theta[t.pos+t.pos.OR.Mult.Time0.Full], (NV.Time0+NV.OdRa),num.bin-1 ))
      Pred0.Exp      <- exp(c(0,apply( Theta.OR.Time0 * as.numeric( data[,d.pos+d.pos.OR.Time0.Full] ),2,sum )))
      Pred0.Prob     <- Pred0.Exp/(sum(Pred0.Exp ))
      
      Design.Control <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos+d.pos.OR.Time0.OdRa])*0 )
      
      Design.Trt     <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos+d.pos.OR.Time0.OdRa])*1 )
      
      Pred0.Exp.Control <- exp(c(0,apply( Theta.OR.Time0 * Design.Control,2,sum )))
      Pred0.Exp.Trt     <- exp(c(0,apply( Theta.OR.Time0 * Design.Trt    ,2,sum )))
      
      Pred0.Prob.Control <- Pred0.Exp.Control/sum(Pred0.Exp.Control)
      Pred0.Prob.Trt     <- Pred0.Exp.Trt/sum(Pred0.Exp.Trt)
      # OdRa0.Vec <- (Pred0.Prob.Trt/Pred0.Prob.Trt[1])/(Pred0.Prob.Control/Pred0.Prob.Control[1])
      
      Mult.Res0 <- rep(0,num.bin-1)
      for(jj in 1:(num.bin-1)){
        Mult.Res0[jj] <- data[,2+jj+1]-Pred0.Prob[jj+1]
      }
      
      Moment.Y.Mult <- rep(Mult.Res0,each=(NV.Time0+NV.OdRa))*rep(as.numeric( data[,d.pos+d.pos.OR.Time0.Full] ),num.bin-1)
      
      ####################################################
      
      Pred.PS.Time0.Full   <- expit( sum( ( Theta[t.pos + t.pos.PS.Time0.Full] ) * as.numeric( data[,d.pos + d.pos.PS.Time0.Full] ) ) )
      Pred.PS.Base         <- expit( sum( ( Theta[t.pos + t.pos.PS.Time0.Base] ) * as.numeric( data[,d.pos + d.pos.PS.Time0.Base] ) ) )
      
      Moment.PS <- (data$Trt - Pred.PS.Time0.Full)*as.numeric( data[,d.pos + d.pos.PS.Time0.Full]  )
      
      
      
      Suff.OdRa           <- ( sum( ( Theta[t.pos + t.pos.OREC] ) * as.numeric( data[,d.pos + d.pos.OREC] )))
      Pred.OdRa.Time0     <- LINK(Suff.OdRa)
      Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
      
      ####################################################
      
      M.matrix      <- as.numeric( data[,d.pos+d.pos.OR.Time0.OdRa] )*as.numeric( data[,2+1:(num.bin-1)] )
      Pred.M.Base0  <- as.numeric( data[,d.pos+d.pos.OR.Time0.OdRa] )*Pred0.Prob.Control[-num.bin]
      
      
      Res.M.Base0      <- M.matrix - Pred.M.Base0
      Res.PS.0.Base    <- data$Trt - Pred.PS.Base
      
      c( # Moment.Y.Mult,
        # Moment.PS,
        ((Res.M.Base0)*( Sandwich.OdRa.Time0 ))*(Res.PS.0.Base) )
      
      
    }
  }
  
  
  
  EESolve.Time0 <- m_estimate(estFUN = EE.Time0,
                              data = data,
                              root_control = setup_root_control(start = ParaStart),
                              compute_vcov = FALSE)
  
  ###############################################################################
  
  
  
  
  data <- data.frame(cbind(Y0,Y1,
                           Disc.Y0.Mat,Disc.Y1.Mat,
                           Trt,
                           Design.Matrix.Time0,         ## OR Time0
                           Design.Matrix.OdRa*Trt,      ## OR Time0, OdRa
                           Design.Matrix.Time0,         ## PS Time0
                           Design.Matrix.OdRa.Y0,       ## PS Time0, OdRa.Y0
                           Design.Matrix.OdRa.Y0,       ## OdRa.Y0
                           Design.Matrix.Time1,         ## OR, PS Time1
                           Design.Matrix.OdRa.Y1,       ## OR, PS OdRa.Y1
                           Design.Matrix.OdRa))         ## OR, PS OdRa 
  
  colnames(data) <- c("Y0","Y1",
                      sprintf("DiscY0_%0.3d",1:num.bin),
                      sprintf("DiscY1_%0.3d",1:num.bin),
                      "Trt",
                      sprintf("OR_BaseT0%0.3d",1:NV.Time0),
                      sprintf("OR_OdRaT0%0.3d",1:NV.OdRa),
                      sprintf("PS_BaseT0%0.3d",1:NV.Time0),
                      sprintf("PS_OdRaT0%0.3d",1:NV.OdRa.Y),
                      sprintf("OREC_OdRaT0Y%0.3d",1:NV.OdRa.Y),
                      sprintf("BaseT1%0.3d",1:NV.Time1),
                      sprintf("OREC_OdRaT1Y%0.3d",1:NV.OdRa.Y),
                      sprintf("OREC_OdRaT1%0.3d",1:NV.OdRa))
  
  theta <- ParaStart <- c( (mean(Y1[Trt==1])-mean(Y1[Trt==0]))-(mean(Y0[Trt==1])-mean(Y0[Trt==0])),
                           (mean(Y1[Trt==1])),
                           (mean(Y1[Trt==0])+(mean(Y0[Trt==1])-mean(Y0[Trt==0]))),
                           mean(LM.Time1$residuals^2),
                           as.numeric(t(coef(GLM.Time1))),
                           (LM.Time1$coefficients),
                           pi.Time1$coefficients[-1] )
  
  
  EE.Time1 <- function(data){
    
    t.pos <- 4
    d.pos <- 3+2*(num.bin)
    
    t.pos.OR.Mult.Time0.Full <- 1:((num.bin-1)*(NV.Time0+NV.OdRa))
    d.pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    d.pos.OR.Time0.Base <- 1:(NV.Time0)
    d.pos.OR.Time0.OdRa <- NV.Time0+1:(NV.OdRa)
    
    t.pos.PS.Time0.Full <- ((num.bin-1)*(NV.Time0+NV.OdRa))+1:(NV.Time0+NV.OdRa.Y)
    t.pos.PS.Time0.Base <- ((num.bin-1)*(NV.Time0+NV.OdRa))+1:(NV.Time0)
    t.pos.PS.Time0.OdRa <- ((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0)+(1:NV.OdRa.Y)
    
    d.pos.PS.Time0.Full <- (NV.Time0+NV.OdRa)+1:(NV.Time0+NV.OdRa.Y)
    d.pos.PS.Time0.Base <- (NV.Time0+NV.OdRa)+1:(NV.Time0)
    d.pos.PS.Time0.OdRa <- (NV.Time0+NV.OdRa)+(NV.Time0)+(1:NV.OdRa.Y)
    
    t.pos.OREC <- ((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0+NV.OdRa.Y)+1:NV.OdRa.Y
    d.pos.OREC <- (NV.Time0+NV.OdRa)+(NV.Time0+NV.OdRa.Y)+1:NV.OdRa.Y
    
    t.pos.Time1 <- t.pos+((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0+NV.OdRa.Y+NV.OdRa.Y)
    d.pos.Time1 <- d.pos+(NV.Time0+NV.OdRa)+(NV.Time0+NV.OdRa.Y)+NV.OdRa.Y
    
    t.pos.OR.Mult.Time1.Base <- 1:((num.bin-1)*(NV.Time1))
    d.pos.OR.Mult.Time1.Base <- 1:(NV.Time1)
    
    t.pos.OR.Lin.Time1.Base <- ((num.bin-1)*(NV.Time1))+(1:NV.Time1)
    d.pos.OR.Lin.Time1.Base <- (1:NV.Time1)
    
    t.pos.PS.Time1.Base <- ((num.bin-1)*(NV.Time1))+NV.Time1+(1:NV.Time1)
    d.pos.PS.Time1.Base <- (1:NV.Time1)
    
    d.pos.PS.Time1.OdRa <- (NV.Time1)+(1:NV.OdRa.Y)
    
    d.pos.Basis.Time1.OdRa <- (NV.Time1)+(NV.OdRa.Y)+(1:NV.OdRa)
    
    function(theta){
      
      Theta <- c(theta[c(1:3)],                                       # ATT, ATT1, ATT0
                 theta[c(4)],                                         # Var(Y_{t=1})
                 as.numeric(t(coef(GLM.Time0.OdRa))),
                 pi.Time0.OdRa$coefficients,
                 coef(EESolve.Time0),
                 theta[-c(1:4)] )
      names(Theta) <- c("ATT","ATT1","ATT0",
                        "Var1",
                        sprintf("Y.Mult.T0_%0.3d",1:((num.bin-1)*(NV.Time0+NV.OdRa))),
                        sprintf("PS.T0_Base_%0.3d",1:NV.Time0),
                        sprintf("PS.T0_OdRa_%0.3d",1:NV.OdRa.Y),
                        sprintf("OdRa.T0_OdRa_%0.3d",1:NV.OdRa.Y),
                        sprintf("Y.Mult.T1_%0.3d",1:((num.bin-1)*(NV.Time1))),
                        sprintf("Y.LM.T1_%0.3d",1:NV.Time1),
                        sprintf("PS.T1_Base_%0.3d",1:NV.Time1))
      
      ## Time 1 Mult 
      
      Theta.OR.Time0 <- (matrix( Theta[t.pos+t.pos.OR.Mult.Time0.Full], (NV.Time0+NV.OdRa),num.bin-1 )) 
      Pred0.Exp      <- exp(c(0,apply( Theta.OR.Time0 * as.numeric( data[,d.pos+d.pos.OR.Time0.Full] ),2,sum )))
      Pred0.Prob     <- Pred0.Exp/(sum(Pred0.Exp ))
      
      Design.Control <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos.Time1+d.pos.Basis.Time1.OdRa])*0 )
      
      Design.Trt     <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos.Time1+d.pos.Basis.Time1.OdRa])*1 )
      
      Pred0.Exp.Control <- exp(c(0,apply( Theta.OR.Time0 * Design.Control,2,sum )))
      Pred0.Exp.Trt     <- exp(c(0,apply( Theta.OR.Time0 * Design.Trt    ,2,sum )))
      
      Pred0.Prob.Control <- Pred0.Exp.Control/sum(Pred0.Exp.Control)
      Pred0.Prob.Trt     <- Pred0.Exp.Trt/sum(Pred0.Exp.Trt)
      OdRa0.Vec <- (Pred0.Prob.Trt/Pred0.Prob.Trt[1])/(Pred0.Prob.Control/Pred0.Prob.Control[1])
      
      Theta.OR.Time1 <- (matrix( Theta[t.pos.Time1+t.pos.OR.Mult.Time1.Base], (NV.Time1),num.bin-1 ))
      Pred1.Exp  <- exp(c(0,apply( Theta.OR.Time1 * as.numeric( data[,d.pos.Time1+d.pos.OR.Mult.Time1.Base] ),2,sum )))
      Pred1.Prob <- Pred1.Exp/(sum(Pred1.Exp ))
      
      Mult.Res1 <- rep(0,num.bin-1)
      for(jj in 1:(num.bin-1)){
        Mult.Res1[jj] <- (data[,2+num.bin+jj+1]-Pred1.Prob[jj+1])*(1-data$Trt)
      }
      
      Moment.Y.Mult.Time1 <- rep(Mult.Res1,each=(NV.Time1))*rep(as.numeric( data[,d.pos.Time1+d.pos.OR.Lin.Time1.Base] ),num.bin-1)
      
      # print(Moment.Y.Mult.Time1)
      
      ## Linear model 1
      
      Pred1.Linear <- sum( ( Theta[t.pos.Time1+t.pos.OR.Lin.Time1.Base] ) * as.numeric( data[,d.pos.Time1+d.pos.OR.Lin.Time1.Base] ) )
      Res1.Linear  <- (data$Y1 - Pred1.Linear)
      Moment.Y.Lin.Time1 <- Res1.Linear*as.numeric( data[,d.pos.Time1+d.pos.OR.Lin.Time1.Base] )*(1-data$Trt)
      
      QN1 <- pnorm(c(Bin.Cut,Inf),Pred1.Linear,Theta[4]^(0.5))
      QN1.R <- (c(QN1[1],diff(QN1)))
      QN2 <- -dnorm(c(Bin.Cut,Inf),Pred1.Linear,Theta[4]^(0.5))
      QN2.R <- (c(QN2[1],diff(QN2)))
      
      Numer <- Pred1.Linear*OdRa0.Vec*QN1.R + Theta[4]*OdRa0.Vec*QN2.R
      Denom <- OdRa0.Vec*QN1.R
      
      
      mu1 <- sum(Numer)/sum(Denom)
      
      ## Propensity time 1
      
      Gap1    <- sum( ( Theta[t.pos + t.pos.PS.Time0.OdRa] ) * as.numeric( data[,d.pos.Time1+d.pos.PS.Time1.OdRa] ) )
      OR1 <- exp( Gap1 )
      Res1.PS <- ( (1-data$Trt)*(1+exp( sum( ( Theta[t.pos.Time1 + t.pos.PS.Time1.Base ] ) *
                                               as.numeric( data[,d.pos.Time1+ d.pos.PS.Time1.Base] ) ))*OR1)-1 )
      
      Moment.PS.Time1 <- as.numeric( data[,d.pos.Time1+ d.pos.PS.Time1.Base] )*Res1.PS
      
      IPW.f <- exp( sum( ( Theta[t.pos.Time1 + t.pos.PS.Time1.Base] ) *
                           as.numeric( data[,d.pos.Time1 + d.pos.PS.Time1.Base] ) ))*OR1
      
      
      c( Theta[1]-(Theta[2]-Theta[3]),
         (data$Trt)*(data$Y1 - Theta[2]), 
         (1-data$Trt)*IPW.f*(data$Y1-mu1) + (data$Trt*mu1) - as.numeric(data$Trt)*Theta[3],
         # (1-data$Trt)*IPW*(data$Y1) - (1-data$Trt)*IPW*Theta[3],
         (1-data$Trt)*((Res1.Linear)^2 - Theta[4]),
         Moment.Y.Mult.Time1,
         Moment.Y.Lin.Time1,
         Moment.PS.Time1)
      
      
    }
  } 
  
  EESolve.Time1 <- m_estimate(estFUN = EE.Time1,
                              data = data,
                              root_control = setup_root_control(start = ParaStart),
                              compute_vcov = FALSE)
  
  ParaStart <- CEE1 <- c(coef(EESolve.Time1)[c(1:3)],                                       # ATT, ATT1, ATT0
                         coef(EESolve.Time1)[c(4)],                                         # Var(Y_{t=1})
                         as.numeric(t(coef(GLM.Time0.OdRa))),
                         pi.Time0.OdRa$coefficients,
                         coef(EESolve.Time0),
                         coef(EESolve.Time1)[-c(1:4)] )
  
  ###############################################################################
  
  EE <- function(data){
    
    t.pos <- 4
    d.pos <- 3+2*(num.bin)
    
    t.pos.OR.Mult.Time0.Full <- 1:((num.bin-1)*(NV.Time0+NV.OdRa))
    d.pos.OR.Time0.Full <- 1:(NV.Time0+NV.OdRa)
    d.pos.OR.Time0.Base <- 1:(NV.Time0)
    d.pos.OR.Time0.OdRa <- NV.Time0+1:(NV.OdRa)
    
    t.pos.PS.Time0.Full <- ((num.bin-1)*(NV.Time0+NV.OdRa))+1:(NV.Time0+NV.OdRa.Y)
    t.pos.PS.Time0.Base <- ((num.bin-1)*(NV.Time0+NV.OdRa))+1:(NV.Time0)
    t.pos.PS.Time0.OdRa <- ((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0)+(1:NV.OdRa.Y)
    
    d.pos.PS.Time0.Full <- (NV.Time0+NV.OdRa)+1:(NV.Time0+NV.OdRa.Y)
    d.pos.PS.Time0.Base <- (NV.Time0+NV.OdRa)+1:(NV.Time0)
    d.pos.PS.Time0.OdRa <- (NV.Time0+NV.OdRa)+(NV.Time0)+(1:NV.OdRa.Y)
    
    t.pos.OREC <- ((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0+NV.OdRa.Y)+1:NV.OdRa.Y
    d.pos.OREC <- (NV.Time0+NV.OdRa)+(NV.Time0+NV.OdRa.Y)+1:NV.OdRa.Y
    
    t.pos.Time1 <- t.pos+((num.bin-1)*(NV.Time0+NV.OdRa))+(NV.Time0+NV.OdRa.Y+NV.OdRa.Y)
    d.pos.Time1 <- d.pos+(NV.Time0+NV.OdRa)+(NV.Time0+NV.OdRa.Y)+NV.OdRa.Y
    
    t.pos.OR.Mult.Time1.Base <- 1:((num.bin-1)*(NV.Time1))
    d.pos.OR.Mult.Time1.Base <- 1:(NV.Time1)
    
    t.pos.OR.Lin.Time1.Base <- ((num.bin-1)*(NV.Time1))+(1:NV.Time1)
    d.pos.OR.Lin.Time1.Base <- (1:NV.Time1)
    
    t.pos.PS.Time1.Base <- ((num.bin-1)*(NV.Time1))+NV.Time1+(1:NV.Time1)
    d.pos.PS.Time1.Base <- (1:NV.Time1)
    
    d.pos.PS.Time1.OdRa <- (NV.Time1)+(1:NV.OdRa.Y)
    
    d.pos.Basis.Time1.OdRa <- (NV.Time1)+(NV.OdRa.Y)+(1:NV.OdRa)
    
    function(theta){
      
      Theta <- theta
      
      
      
      Theta.OR.Time0 <- (matrix( Theta[t.pos+t.pos.OR.Mult.Time0.Full], (NV.Time0+NV.OdRa),num.bin-1 ))
      Pred0.Exp      <- exp(c(0,apply( Theta.OR.Time0 * as.numeric( data[,d.pos+d.pos.OR.Time0.Full] ),2,sum )))
      Pred0.Prob     <- Pred0.Exp/(sum(Pred0.Exp ))
      
      Design.Control <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos+d.pos.OR.Time0.OdRa])*0 )
      
      Design.Trt     <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos+d.pos.OR.Time0.OdRa])*1 )
      
      Pred0.Exp.Control <- exp(c(0,apply( Theta.OR.Time0 * Design.Control,2,sum )))
      Pred0.Exp.Trt     <- exp(c(0,apply( Theta.OR.Time0 * Design.Trt    ,2,sum )))
      
      Pred0.Prob.Control <- Pred0.Exp.Control/sum(Pred0.Exp.Control)
      Pred0.Prob.Trt     <- Pred0.Exp.Trt/sum(Pred0.Exp.Trt)
      # OdRa0.Vec <- (Pred0.Prob.Trt/Pred0.Prob.Trt[1])/(Pred0.Prob.Control/Pred0.Prob.Control[1])
      
      Mult.Res0 <- rep(0,num.bin-1)
      for(jj in 1:(num.bin-1)){
        Mult.Res0[jj] <- data[,2+jj+1]-Pred0.Prob[jj+1]
      }
      
      Moment.Y.Mult <- rep(Mult.Res0,each=(NV.Time0+NV.OdRa))*rep(as.numeric( data[,d.pos+d.pos.OR.Time0.Full] ),num.bin-1)
      
      ####################################################
      
      Pred.PS.Time0.Full   <- expit( sum( ( Theta[t.pos + t.pos.PS.Time0.Full] ) * as.numeric( data[,d.pos + d.pos.PS.Time0.Full] ) ) )
      Pred.PS.Base         <- expit( sum( ( Theta[t.pos + t.pos.PS.Time0.Base] ) * as.numeric( data[,d.pos + d.pos.PS.Time0.Base] ) ) )
      
      Moment.PS <- (data$Trt - Pred.PS.Time0.Full)*as.numeric( data[,d.pos + d.pos.PS.Time0.Full]  )
      
      
      
      Suff.OdRa           <- ( sum( ( Theta[t.pos + t.pos.OREC] ) * as.numeric( data[,d.pos + d.pos.OREC] )))
      Pred.OdRa.Time0     <- LINK(Suff.OdRa)
      Sandwich.OdRa.Time0 <- Pred.OdRa.Time0^(-data$Trt)
      
      ####################################################
      
      M.matrix      <- as.numeric( data[,d.pos+d.pos.OR.Time0.OdRa] )*as.numeric( data[,2+1:(num.bin-1)] )
      Pred.M.Base0  <- as.numeric( data[,d.pos+d.pos.OR.Time0.OdRa] )*Pred0.Prob.Control[-num.bin]
      
      
      Res.M.Base0      <- M.matrix - Pred.M.Base0
      Res.PS.0.Base    <- data$Trt - Pred.PS.Base
      
      
      
      ## Time 1 Mult 
      
      Theta.OR.Time0 <- (matrix( Theta[t.pos+t.pos.OR.Mult.Time0.Full], (NV.Time0+NV.OdRa),num.bin-1 ))
      Pred0.Exp      <- exp(c(0,apply( Theta.OR.Time0 * as.numeric( data[,d.pos+d.pos.OR.Time0.Full] ),2,sum )))
      Pred0.Prob     <- Pred0.Exp/(sum(Pred0.Exp ))
      
      Design.Control <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos.Time1+d.pos.Basis.Time1.OdRa])*0 )
      
      Design.Trt     <- c( as.numeric(data[,d.pos+d.pos.OR.Time0.Base]),
                           as.numeric(data[,d.pos.Time1+d.pos.Basis.Time1.OdRa])*1 )
      
      Pred0.Exp.Control <- exp(c(0,apply( Theta.OR.Time0 * Design.Control,2,sum )))
      Pred0.Exp.Trt     <- exp(c(0,apply( Theta.OR.Time0 * Design.Trt    ,2,sum )))
      
      Pred0.Prob.Control <- Pred0.Exp.Control/sum(Pred0.Exp.Control)
      Pred0.Prob.Trt     <- Pred0.Exp.Trt/sum(Pred0.Exp.Trt)
      OdRa0.Vec <- (Pred0.Prob.Trt/Pred0.Prob.Trt[1])/(Pred0.Prob.Control/Pred0.Prob.Control[1])
      
      Theta.OR.Time1 <- (matrix( Theta[t.pos.Time1+t.pos.OR.Mult.Time1.Base], (NV.Time1),num.bin-1 ))
      Pred1.Exp  <- exp(c(0,apply( Theta.OR.Time1 * as.numeric( data[,d.pos.Time1+d.pos.OR.Mult.Time1.Base] ),2,sum )))
      Pred1.Prob <- Pred1.Exp/(sum(Pred1.Exp ))
      
      Mult.Res1 <- rep(0,num.bin-1)
      for(jj in 1:(num.bin-1)){
        Mult.Res1[jj] <- (data[,2+num.bin+jj+1]-Pred1.Prob[jj+1])*(1-data$Trt)
      }
      
      Moment.Y.Mult.Time1 <- rep(Mult.Res1,each=(NV.Time1))*rep(as.numeric( data[,d.pos.Time1+d.pos.OR.Lin.Time1.Base] ),num.bin-1)
      
      # print(Moment.Y.Mult.Time1)
      
      ## Linear model 1
      
      Pred1.Linear <- sum( ( Theta[t.pos.Time1+t.pos.OR.Lin.Time1.Base] ) * as.numeric( data[,d.pos.Time1+d.pos.OR.Lin.Time1.Base] ) )
      Res1.Linear  <- (data$Y1 - Pred1.Linear)
      Moment.Y.Lin.Time1 <- Res1.Linear*as.numeric( data[,d.pos.Time1+d.pos.OR.Lin.Time1.Base] )*(1-data$Trt)
      
      QN1 <- pnorm(c(Bin.Cut,Inf),Pred1.Linear,Theta[4]^(0.5))
      QN1.R <- (c(QN1[1],diff(QN1)))
      QN2 <- -dnorm(c(Bin.Cut,Inf),Pred1.Linear,Theta[4]^(0.5))
      QN2.R <- (c(QN2[1],diff(QN2)))
      
      Numer <- Pred1.Linear*OdRa0.Vec*QN1.R + Theta[4]*OdRa0.Vec*QN2.R
      Denom <- OdRa0.Vec*QN1.R
      
      mu1 <- sum(Numer)/sum(Denom)
      
      ## Propensity time 1
      
      Gap1    <- sum( ( Theta[t.pos + t.pos.PS.Time0.OdRa] ) * as.numeric( data[,d.pos.Time1+d.pos.PS.Time1.OdRa] ) )
      OR1 <- exp( Gap1 )
      Res1.PS <- ( (1-data$Trt)*(1+exp( sum( ( Theta[t.pos.Time1 + t.pos.PS.Time1.Base ] ) *
                                               as.numeric( data[,d.pos.Time1+ d.pos.PS.Time1.Base] ) ))*OR1)-1 )
      
      Moment.PS.Time1 <- as.numeric( data[,d.pos.Time1+ d.pos.PS.Time1.Base] )*Res1.PS
      
      IPW <- exp( sum( ( Theta[t.pos.Time1 + t.pos.PS.Time1.Base] ) *
                         as.numeric( data[,d.pos.Time1 + d.pos.PS.Time1.Base] ) ))*OR1
      
      c( Theta[1]-(Theta[2]-Theta[3]),
         (data$Trt)*(data$Y1 - Theta[2]), 
         ((1-data$Trt)*IPW*(data$Y1-mu1) + data$Trt*mu1) - data$Trt*Theta[3],
         (1-data$Trt)*((Res1.Linear)^2 - Theta[4]),
         Moment.Y.Mult,
         Moment.PS,
         ((Res.M.Base0)*( Sandwich.OdRa.Time0 ))*(Res.PS.0.Base),
         Moment.Y.Mult.Time1,
         Moment.Y.Lin.Time1,
         Moment.PS.Time1)
      
      
    }
  } 
  
  EESolve <- m_estimate(estFUN = EE,
                        data = data,
                        root_control = setup_root_control(start = ParaStart))
  
  
  
  # Gap1    <- as.matrix(data[,d.pos.Time1+d.pos.PS.Time1.OdRa])%*%Theta[t.pos + t.pos.PS.Time0.OdRa]
  # OR1 <- exp( Gap1 )
  # Res1.PS <- ( (1-data$Trt)*(1+exp( as.matrix(data[,d.pos.Time1+ d.pos.PS.Time1.Base])%*%
  #                                     Theta[t.pos.Time1 + t.pos.PS.Time1.Base ])*OR1)-1 )
  # 
  # Moment.PS.Time1 <- as.matrix( data[,d.pos.Time1+ d.pos.PS.Time1.Base] )*Res1.PS
  # 
  # IPW <- exp(  as.matrix( data[,d.pos.Time1 + d.pos.PS.Time1.Base] ) %*%Theta[t.pos.Time1 + t.pos.PS.Time1.Base] )*OR1
  # 
  # mean(IPW*(1-data$Trt))
  
  # EESolve <- m_estimate(estFUN = EE,
  #                       roots = ParaStart,
  #                       compute_roots=FALSE,
  #                       data  = data)
  
  CEE.Full <- coef(EESolve)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0")
  
  CVar.Full <- vcov(EESolve)
  colnames(CVar.Full) <- colnames(CEE.Full)
  
  
  Result <- list()
  Result$Est <- as.numeric( coef(EESolve)[1] )
  Result$SE  <- as.numeric( sqrt(vcov(EESolve)[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full 
  
  return(Result)
  
}



