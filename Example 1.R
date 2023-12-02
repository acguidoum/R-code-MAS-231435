###############################################################
###############################################################
###############################################################
### Example 1: case g(t) is linear time-dependent function ####
###############################################################
###############################################################
###############################################################
####################
## A: alpha
## B: beta
## G: gamma
## g_t: g(t)
####################

set.seed(412321)

## g(t) == 4*t+35
X0=0; g_t <- expression(4*t+35); K=0
A=-2; B=1; G=4 

## Simulation model (8), t \in [0,50]

mod1 <- Sim_mod(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,T=50,N=1000, M=100000)
Time <- as.vector(time(mod1))

## Monte-Carlo moments approximation

MC_mom_n1 <- apply(mod1,1,Sim.DiffProc::moment,order=1,center=FALSE)
MC_mom_n2 <- apply(mod1,1,Sim.DiffProc::moment,order=2,center=FALSE)
MC_mom_n3 <- apply(mod1,1,Sim.DiffProc::moment,order=3,center=FALSE)
MC_mom_n4 <- apply(mod1,1,Sim.DiffProc::moment,order=4,center=FALSE)
MC_mom_n5 <- apply(mod1,1,Sim.DiffProc::moment,order=5,center=FALSE)
MC_mom_n6 <- apply(mod1,1,Sim.DiffProc::moment,order=6,center=FALSE)
MC_mom_n7 <- apply(mod1,1,Sim.DiffProc::moment,order=7,center=FALSE)
MC_mom_n8 <- apply(mod1,1,Sim.DiffProc::moment,order=8,center=FALSE)
MC_mom_n9 <- apply(mod1,1,Sim.DiffProc::moment,order=9,center=FALSE)
MC_mom_n10<- apply(mod1,1,Sim.DiffProc::moment,order=10,center=FALSE)


MC_mom_scaled <- data.frame(MC_mom_n1,MC_mom_n2^(1/2),MC_mom_n3^(1/3),
                            MC_mom_n4^(1/4),MC_mom_n5^(1/5),MC_mom_n6^(1/6),
                            MC_mom_n7^(1/7),MC_mom_n8^(1/8),MC_mom_n9^(1/9),
                            MC_mom_n10^(1/10))

## Exact moments

Ex_mom_n1 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=1,t=Time)
Ex_mom_n2 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=2,t=Time)
Ex_mom_n3 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=3,t=Time)
Ex_mom_n4 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=4,t=Time)
Ex_mom_n5 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=5,t=Time)
Ex_mom_n6 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=6,t=Time)
Ex_mom_n7 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=7,t=Time)
Ex_mom_n8 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=8,t=Time)
Ex_mom_n9 <- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=9,t=Time)
Ex_mom_n10<- HO_mom(x0=X0,gt=g_t,k=K,Alpha=A,Beta=B,Gamma=G,n=10,t=Time)

Ex_mom_scaled <- data.frame(Ex_mom_n1,Ex_mom_n2^(1/2),Ex_mom_n3^(1/3),
                            Ex_mom_n4^(1/4),Ex_mom_n5^(1/5),Ex_mom_n6^(1/6),
                            Ex_mom_n7^(1/7),Ex_mom_n8^(1/8),Ex_mom_n9^(1/9),
                            Ex_mom_n10^(1/10))

## Plots Fig 01

Color <- ggsci::pal_jco("default",alpha=0.8)(10)
Expr <- expression(M[1](t),M[2](t),M[3](t),M[4](t),M[5](t),M[6](t),
        M[7](t),M[8](t),M[9](t),M[10](t))
par(mar = c(4, 3, 1, 1) + 0.1)
matplot(Time,MC_mom_scaled,type="l",lty=1,las=1,col=Color,lwd=2,
        ylab="Moment",xlab="Time t")
matplot(Time,Ex_mom_scaled,type="l",lty=2,add=TRUE,col=1,lwd=2)
legend("topright",Expr,inset = .01,fill=Color,lty=NA,
        lwd=2,cex=1.35)
legend("topright",expression(m[n](t),E(X[t]^{n})),inset = c(0.32,0),
       col=1,lty=c(1,2),lwd=2,cex=1.35,horiz = FALSE, merge = TRUE,bty="n")

