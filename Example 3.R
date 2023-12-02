########################################################################
########################################################################
########################################################################
### Example 3: case g(t) multiple nonlinear time-dependent function ####
########################################################################
########################################################################
########################################################################
####################
## A: alpha
## B: beta
## G: gamma
## g_t: g(t)
####################

set.seed(412321)

## Simulation model (10), t \in [0,50]

X0=0;T=50
A=-0.5;B=0.5;G=8

## the second derivative of g(t)

g_t <- expression(2*(t <= T/2)+(2*t*(t^2 + 3*t + 3)/(t + 1)^3)*(t > T/2))
mod3 <- Sim_mod(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,T=50,N=1000, M=100000)
Time <- as.vector(time(mod3))

## Monte-Carlo moments approximation

MC_mom_n1 <- apply(mod3,1,Sim.DiffProc::moment,order=1,center=FALSE)
MC_mom_n2 <- apply(mod3,1,Sim.DiffProc::moment,order=2,center=FALSE)
MC_mom_n3 <- apply(mod3,1,Sim.DiffProc::moment,order=3,center=FALSE)
MC_mom_n4 <- apply(mod3,1,Sim.DiffProc::moment,order=4,center=FALSE)
MC_mom_n5 <- apply(mod3,1,Sim.DiffProc::moment,order=5,center=FALSE)
MC_mom_n6 <- apply(mod3,1,Sim.DiffProc::moment,order=6,center=FALSE)
MC_mom_n7 <- apply(mod3,1,Sim.DiffProc::moment,order=7,center=FALSE)
MC_mom_n8 <- apply(mod3,1,Sim.DiffProc::moment,order=8,center=FALSE)
MC_mom_n9 <- apply(mod3,1,Sim.DiffProc::moment,order=9,center=FALSE)
MC_mom_n10<- apply(mod3,1,Sim.DiffProc::moment,order=10,center=FALSE)

MC_mom_scaled <- data.frame(MC_mom_n1,MC_mom_n2^(1/2),MC_mom_n3^(1/3),
                            MC_mom_n4^(1/4),MC_mom_n5^(1/5),MC_mom_n6^(1/6),
                            MC_mom_n7^(1/7),MC_mom_n8^(1/8),MC_mom_n9^(1/9),
                            MC_mom_n10^(1/10))

## Exact moments

Ex_mom_n1 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=1,t=Time)
Ex_mom_n2 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=2,t=Time)
Ex_mom_n3 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=3,t=Time)
Ex_mom_n4 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=4,t=Time)
Ex_mom_n5 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=5,t=Time)
Ex_mom_n6 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=6,t=Time)
Ex_mom_n7 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=7,t=Time)
Ex_mom_n8 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=8,t=Time)
Ex_mom_n9 <- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=9,t=Time)
Ex_mom_n10<- HO_mom(x0=X0,gt=g_t,Alpha=A,Beta=B,Gamma=G,n=10,t=Time)

Ex_mom_scaled <- data.frame(Ex_mom_n1,Ex_mom_n2^(1/2),Ex_mom_n3^(1/3),
                            Ex_mom_n4^(1/4),Ex_mom_n5^(1/5),Ex_mom_n6^(1/6),
                            Ex_mom_n7^(1/7),Ex_mom_n8^(1/8),Ex_mom_n9^(1/9),
                            Ex_mom_n10^(1/10))

## Plots Fig 03

Color <- ggsci::pal_jco("default",alpha=0.8)(10)
Expr <- expression(M[1](t),M[2](t),M[3](t),M[4](t),M[5](t),M[6](t),
        M[7](t),M[8](t),M[9](t),M[10](t))
par(mar = c(4, 3, 1, 1) + 0.1)
matplot(Time,MC_mom_scaled,type="l",lty=1,las=1,col=Color,lwd=2,
        ylab="Moment",xlab="Time t")
matplot(Time,Ex_mom_scaled,type="l",lty=2,add=TRUE,col=1,lwd=2)
legend("topleft",Expr,inset = .01,fill=Color,lty=NA,
        lwd=2,cex=1.35)
legend("topleft",expression(m[n](t),E(X[t]^{n})),inset = c(0.32,0),
       col=1,lty=c(1,2),lwd=2,cex=1.35,horiz = FALSE, merge = TRUE,bty="n")
axis(1, c(50/2), expression(T/2),padj=0.1,font.lab=3)
segments(25, 0, 25, 58,lwd=1,lty=4)


## explicitly first, second, and third moments of (10)

plot(Time,Ex_mom_n1,type="l")
trueE1 <- function(t) (sqrt(2)/4)*(t*(t <= T/2)+(sqrt((t*(t+1)^3)/(t^2 +3*t+3)))*(t >T/2))
curve(trueE1,xname="t",add=TRUE,col=2)

plot(Time,Ex_mom_n2,type="l")
trueE2 <- function(t) ((t+256)/8)*(t*(t <= T/2)+(((t+1)^3)/(t^2 +3*t+3))*(t >T/2))
curve(trueE2,xname="t",add=TRUE,col=2)

plot(Time,Ex_mom_n3,type="l")
trueE3 <- function(t) ((sqrt(2)*t*(t+768))/32)*(t*(t <= T/2)+sqrt((t+1)^9 /(t*(t^2+3*t+3)^3))*(t >T/2))
curve(trueE3,xname="t",add=TRUE,col=2)

