######################################################################
##                    Sat Dec  2 01:24:15 2023                       #
##                  A.C. GUIDOUM and K. BOUKHETALA                   #
######################################################################
##                Supplementary material for                         #
##  Exact higher-order moments for linear non-homogeneous            # 
##              stochastic differential equation                     #                        
######################################################################
######################################################################
##                 readers can be replicate the results              #
######################################################################

## load "Sim.DiffProc" package and others packages

if (!require("Sim.DiffProc")) {
    install.packages("Sim.DiffProc")
    library("Sim.DiffProc")
}

if (!require("ttutils")) {
    install.packages("ttutils")
    library("ttutils")
}

if (!require("Deriv")) {
    install.packages("Deriv")
    library("Deriv")
}

if (!require("EQL")) {
    install.packages("EQL")
    library("EQL")
}


if (!require("ggsci")) {
    install.packages("ggsci")
    library("ggsci")
}


###########################################
### Simulation model Eq.(2): R function ###
###########################################

Sim_mod <- function(x0,gt,k=NULL,Alpha=1,Beta=1,Gamma=1,T=1,N=1000, M=100000)
           {
   if (Gamma <= 0) stop("'Gamma' must be numeric > 0")
   if (!is.expression(gt)) stop(" 'g(t)' must be expressions in 't'")
   if (is.null(k)){ Gt <- gt}else{
   if (ttutils::isInteger(k)==FALSE) stop("k must be non-negative integer.")
   Gt  <- Deriv::Deriv(gt,"t",nderiv = k)}
   G_t  <- function(t) eval(Gt)
   W <- Sim.DiffProc::BM(T = T, N = N, M = M)
   X <- (G_t(time(W)))^(Alpha)*(x0*G_t(0)^(-Alpha)+Beta*time(W)+Gamma*W)
   name <- "X"
   name <- if(M > 1) paste("X",1:M,sep="")
   X <- ts(X, start = 0, deltat = deltat(W), names=name)
   return(invisible(X))
}

###################################################
### Exact high-order moments Eq.(5): R function ###
###################################################


HO_mom <- function(x0,gt,k=NULL,Alpha=1,Beta=1,Gamma=1,n=1,t)
        {
    t0 <- t[which(t==0)]
    t  <- t[which(t!=0)]
    if (ttutils::isInteger(n)==FALSE) stop("n must be non-negative integer.")
    if (Gamma <= 0) stop("'Gamma' must be numeric > 0")
    if (!is.expression(gt)) stop(" 'g(t)' must be expressions in 't'")
    if (is.null(k)){Gt <- gt}else{
    if (ttutils::isInteger(k)==FALSE) stop("k must be non-negative integer.")
    Gt  <- Deriv::Deriv(gt,"t",nderiv = k)}
    G_t <- function(t) eval(Gt)
    Psi   <- function(t) 2^(-0.5*n)*((Gamma*G_t(t)^(Alpha))/1i)^n
    Theta <- function(t) ((G_t(0)^(-Alpha)*x0+Beta*t )*1i)/(Gamma*sqrt(2*t))
    H     <- EQL::hermite(x=Theta(t),n=n,prob=FALSE)
    result <- Re(Psi(t)*H*t^(0.5*n))
    if (length(t0)==1 & t0==0) {result <- c(x0^n,result)} else { result <- result}
    return(invisible(result))
}
