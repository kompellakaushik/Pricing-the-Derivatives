## Black-Scholes model parameters
## s0  - current stock price
## k   - strike price
## sig - pricing volatility
## t   - time, in years
## q   - dividend yield if any, default is no dividend yield
## typ - typ of option: 'c' for call and 'p' for put
## 
## Comment
## - The variables s0, k, sig, t and q can be given as a scalar or a 
##   vector; however, only one of the inputs can be a vector and the
##   rest must then be scalars

## BSd computes d1 and d2 and the cummulative probability given each.
## This is used by the functions below.
BSd <- function(s0,k,sig,r,t,q=0) {
  d1 <- (log(s0/k)+(r-q+sig^2/2)*t)/(sig*sqrt(t))
  d2 <- d1 - sig*sqrt(t)
  return(data.frame(d1=d1,d2=d2,pn1=pnorm(d1),
                    pn2=pnorm(d2),pn3=pnorm(-d2),pn4=pnorm(-d1)))
}

## BSN1 computes the first derivative of the cummulative probability
## and is used by the functions below
BSN1 <- function(x) {
  return(exp(-(x^2)/2)/sqrt(2*pi))
}

## BS computes the price of call and put options
BS <- function(s0,k,sig,r,t=1,q=0,typ='c') {
  d <- BSd(s0,k,sig,r,t,q)
  if(typ=='c')
    return( s0*exp(-q*t)*d$pn1 - k*exp(-r*t)*d$pn2 )
  else
    return( k*exp(-r*t)*d$pn3 - s0*exp(-q*t)*d$pn4 )
}

## BSdelta computes the DELTA of call and put options
BSdelta <- function(s0=50,k=50,sig=.3,r=.05, t=1,q=0,typ='c') {
  if(typ=='c')
    return(BSd(s0,k,sig,r,t,q)$pn1*exp(-q*t))
  else
    return((BSd(s0,k,sig,r,t,q)$pn1-1)*exp(-q*t))
}

## BSgamma computes the GAMMA of call and put options  
BSgamma <- function(s0=50,k=50,sig=.3,r=.05,t=1,q=0,typ='c') {
  return(BSN1(BSd(s0,k,sig,r,t,q)$d1)*exp(-q*t)/(s0*sig*sqrt(t)))
}

## BSvega computes the VEGA of call and put options
BSvega <- function(s0=50,k=50,sig=.3,r=.05,t=1,q=0,typ='c') {
  return(s0*BSN1(BSd(s0,k,sig,r,t,q)$d1)*sqrt(t)*exp(-q*t)/100)
}

## BStheta computes the THETA of call and put options
BStheta <- function(s0=50,k=50,sig=.3,r=.05,t=1,q=0,typ='c') {
  d <- BSd(s0,k,sig,r,t,q)
  if(typ=='c') {
    return(q*s0*d$pn1*exp(-q*t) - r*k*exp(-r*t)*d$pn2
           -s0*BSN1(d$d1)*sig*exp(-q*t)/(2*sqrt(t)))
  } else {
    return(-q*s0*d$pn4*exp(-q*t) + r*k*exp(-r*t)*d$pn3
           -s0*BSN1(d$d1)*sig*exp(-q*t)/(2*sqrt(t)))
  }
}

## BSrho computes the RHO of call and put options
BSrho <- function(s0=50,k=50,sig=.3,r=.05,t=1,q=0,typ='c') {
  d <- BSd(s0,k,sig,r,t,q)
  if(typ=='c')
    return(k*t*exp(-r*t)*d$pn2)
  else
    return(-k*t*exp(-r*t)*d$pn3)  
}

## BSpsi computes the PSI of call and put options
BSpsi <- function(s0=50,k=50,sig=.3,r=.05,t=1,q=0,typ='c') {
  d <- BSd(s0,k,sig,r,t,q)
  if(typ=='c')
    return(-s0*t*exp(-q*t)*d$pn2)
  else
    return(s0*t*exp(-q*t)*d$pn3)  
}

## ImpVol - immplied volatility
##          'p' is the first input and is the price of the option, don't include
##          a volatility input because that's an output of this function.
##          Function handles both Black-Scholes and binomial pricing models.
##          Add another input 'sty' which is used in the case of computing 
##          implied volatilities using the binomial model.  'A' is option is 
##          American and European for anything else.
ImpVol <- function(p,s0=50,k=50,r=.05,t=1,q=0,typ='c',N=NA,mod='base',sty='A') {
  if(is.na(N)) {
    f <- function(sig,p,s0,k,r,t,q,typ) {
      return((p-BS(s0=s0,k=k,sig=sig,r=r,t=t,q=q,typ=typ))^2)
      }
    return(optimize(f,c(0,5),p=p,s0=s0,k=k,r=r,t=t,q=q,typ=typ)$minimum)
  } else {
    if(sty=='A') {
      f <- function(sig,p,s0,k,r,t,q,typ,N,mod) {
        return((p-BT(s0=s0,k=k,sig=sig,r=r,t=t,q=q,typ=typ,N=N,
                     mod=mod)$Res$PriceA)^2)
        }
    } else {
      f <- function(sig,p,s0,k,r,t,q,typ,N,mod) {
        return((p-BT(s0=s0,k=k,sig=sig,r=r,t=t,q=q,typ=typ,N=N,
                     mod=mod)$Res$PriceE)^2)
      }
    }
    return(optimize(f,c(0,5),p=p,s0=s0,k=k,r=r,t=t,q=q,typ=typ,
                    N=N,mod=mod)$minimum)
  }
}

## Function inputs:
## s0     Current asset spot price
## K      Strike price
## sig    Annualized volatility given as a decimal, i.e., 30% = 0.3
## r      Annualized interest rate continuously compounded as a decimal
## t      Time in years, e.g., 3 months, t=0.25, 5 days, t = 5/365
## q      Annualized dividend yield given as a decimal
## typ    Option type: 'c' if a call option, otherwise it's a put option
## N      Binomial moodel only - number of time steps, h = t/N is the amount
##        time for each step
## mod    Binomial model only - forward stock price diffusion process:
##            base -  exponentiate the forward movements by r-q times h plus
##                    volatility times the squared root of h for up movements
##                    and less the incremental volatility for down movements.
##            crr -   Cox-ross-Rubinstein - up movements are just the 
##                    expontiated volatility times the squared root of h, and
##                    down movements are just 1/up
##            Lognormal - Jarrow & Rudd, similar to the base, except the 
##                    the up movement is adjusted by r-q less half the variance
## Function outputs: BT returns a list consisting of the results, the 
##                    price tree given as the upper diagonal of the matrix,
##                    the respective European and American option values on each
##                    node.  Plot tree takes as an input the output list
##                    from BT to generate a graphical respresenation of the
##                    result for steps N < 31.  Other functions return prices
##                    or other metrics depending on the function.

BT <- function(s0,k,sig,r,t,q=0,typ='c',N=50,mod='base') {
  tr <- valE  <- valA <- array(numeric(0),c(N+1,N+1))
  dt <- t/N
  ## Pricing a call (typ = 'c'), or put (typ = 'p')
  sgn <- ifelse(typ == 'c', 1, -1)
  ## Determine the volatility model to use: 'base', 'ccr', or Jarrow Rudd
  ## for any other input
  if(mod=='base') {
    u <- exp((r-q)*dt+sig*sqrt(dt))
    d <- exp((r-q)*dt-sig*sqrt(dt))
  } else if(mod=='crr') {  ## Cox-Ross-Rubinstein
    u <- exp(sig*sqrt(dt))
    d <- 1/u
  } else {                 ## Lognormal Tree: Jarrow and Rudd
    u <- exp((r-q-sig^2/2)*dt+sig*sqrt(dt))
    d <- exp((r-q-sig^2/2)*dt-sig*sqrt(dt))
  }
  tr[1,] <- s0*u^(0:N)
  for(i in 2:(N+1)) {
    tr[2:i,i] <- tr[1:(i-1),i-1]*d
  }
  a <- exp((r-q)*dt)
  pu <- (a-d)/(u-d)
  pd <- 1 - pu
  pv <- exp(-r*dt)
  valE[,N+1] <- valA[,N+1] <- ifelse(sgn*(tr[,N+1]-k) > 0, sgn*(tr[,N+1]-k),0)
  for(i in N:1) {
    valE[1:i,i] <- pv * (pu * valE[1:i,i+1] + pd * valE[2:(i+1),i+1])
    tval <- pv * (pu * valA[1:i,i+1] + pd * valA[2:(i+1),i+1])
    ex <- sgn*(tr[1:i,i] - k)
    valA[1:i,i] <- ifelse(ex > tval, ex, tval)
  }
  ## Model returns results including pricing assumptions, in addition, it 
  ## returns asset price, European and American value pricing trees
  return(list(Res=data.frame(PriceE=valE[1,1],PriceA=valA[1,1], u=u,d=d, pu=pu,
                             pd=pd,pv=pv,K=k,T=t,S=s0,D=q,sig=sig,r=r,
                             typ=typ),
              Prices=tr,ValE=valE,ValA=valA))
}

PlotTree <- function(M,Title='Option Value',val='Option',digits=3,X11=T,
                     typ='A') {
  ## M      - Data object from BT (binomail tree), containing all data used 
  ## Title  - User can input another base title
  ## val    - Is the tree for option price tree or underlying asset price
  ## digits - number of digits to report on node prices
  ## X11    - T if you wish to open a graph windown (recommended), or F
  ##          if you want to plot within GUI program (not recommended)
  ## typ    - 'A' if you want to report the value of the American options,
  ##          anything else if you wish to report the prices for Europeans
  if(X11) 
    x11()
  omai <- par()$mai
  par(mai=c(.7,.5,.7,0.2))
  if(val=='Option') {
    if(typ=='A')
      Tree <- round(M$ValA,digits)
    else
      Tree <- round(M$ValE,digits)
    cex <- 1.6/log(1+ncol(Tree))
  } else {
    Tree <- round(M$Prices,digits)
    cex <- 1.6/log(1+ncol(Tree))
  }
  cex <- 1.6/log(1+ncol(Tree))
  if(ncol(Tree)>31)
    stop("Maximum number of steps is 30")
  if(ncol(Tree)<2)
    stop("Minimum number of steps is 2")
  dx <- -0.025
  dy <- ifelse(ncol(Tree)/30>0.4,ncol(Tree)/30,0.4)
  depth  <- ncol(Tree)
  K <- M$Res$K
  sgn <- ifelse(M$Res$typ == 'c', 1, -1)
  cs <- format(M$Res$S,digits=2)
  ck <- format(K,digits=2)
  ct <- format(M$Res$T,digits=4)
  cr <- format(M$Res$r,digits=2)
  csig <- format(M$Res$sig,digits=2)
  cq <- format(M$Res$D,digits=2)
  stitle <- substitute(paste(Title,': S'[0],' = ',cs,', K = ',ck,', T = ',
                             ct,', ', sigma, ' = ', csig,', r = ',cr,', ',
                             delta,' = ',cq))
  plot(x = c(0, depth-1), y = c(-depth + 1, depth - 1), type = "n", col = 0,
       main=stitle, ylab='',xlab='',yaxt='n',
       ylim=c(-depth+1.2,depth-.4),xlim=c(-.1,depth-1))
  if(val=='Option') {
    if(typ=='A') {
      mtext('American Option Value', side=2, line=1, cex=1)
    } else 
      mtext('European Option Value', side=2, line=1, cex=1)
  } else
    mtext('Stock Price', side=2, line=1, cex=1)
  mtext('Step', side=1, line=2, cex=1)
  points(x = 0, y = 0, pch=19, col='darkgoldenrod',cex=.8)
  text(0 + dx, 0 + dy, deparse(Tree[1,1]), cex = cex,col='darkblue',font=2)
  for (i in 1:(depth - 1)) {
    y = (-i):i
    x = rep(c(i + 1, i), times = 2 * i)[1:length(y)]
    lines(x-1, y, col = 'darkgoldenrod',cex=1.2)
    y = seq(from = -i, by = 2, length = i + 1)
    x = rep(i, times = length(y)) + 1
    points(x-1, y, pch=19,col='darkgoldenrod',cex=.8)
    for (j in 1:length(x))
      text(x[j]-1 + dx, y[j] + dy, 
           deparse(Tree[length(x) + 1 - j, i + 1]), cex = cex,
           col=ifelse(sgn*(M$Prices[length(x)+1-j,i+1]-K)>0,'forestgreen',
                      'brown1'),font=2)
  }
  par(mai=omai)
}

