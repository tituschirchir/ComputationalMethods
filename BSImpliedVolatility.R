library(curl)
library(quantmod)
library(scatterplot3d)

today <- "2017-02-10"
tolerance <- 0.0001

#Black Scholes option pricer ----
BlackScholesPricer <- function(isCall, S, K, tau, r, vol)
{
  d1 <-   getD1(S, K, tau, r, vol)
  d2 <- d1 - vol * sqrt(tau)
  ifelse(
    isCall,
    S * pnorm(d1) - K * exp(-r * tau) * pnorm(d2),
    K * exp(-r * tau) * pnorm(-d2) - S * pnorm(-d1)
  )
}
#Put-Call Parity Checker ----
PutCallParityCheck <- function(S, K, tau, r, vol)
{
  C <- BlackScholesPricer(T, S, K, tau, r, vol)
  P <- BlackScholesPricer(F, S, K, tau, r, vol)
  isMet <- C - P == S - K * exp(-r * tau)
  equalswithinTolerance(C - P, S - K * exp(-r * tau))
}

#Implied Volatility Calculator: Will decide between bisection and secant based on flag ----
calcImpliedVols <-
  function(isCall,
           isBisection,
           Cmkt,
           S,
           K,
           tau,
           r,
           tolerance)
  {
    ifelse(
      isBisection,
      bisectionImpliedVols(isCall, Cmkt, S, K, tau, r, tolerance, 2000),
      secantImpliedVols(
        isCall,
        Cmkt,
        S,
        K,
        tau,
        r,
        c(0, 3),
        BlackScholesPricer(isCall, S, K, tau, r, 1)
      )
    )
  }

#Calculate implied volatilities using bisection method ----
bisectionImpliedVols <-
  function(isCall,
           Cmkt,
           S,
           K,
           tau,
           r,
           tolerance,
           maxIter)
  {
    a <- 0
    b <- 5
    CHigh <- Cmkt - BlackScholesPricer(isCall, S, K, tau, r, b)
    CLow <- Cmkt - BlackScholesPricer(isCall, S, K, tau, r, a)
    impliedVol <- 0.0
    midp <- 0
    midCdiff <-
      Cmkt - BlackScholesPricer(isCall, S, K, tau, r, (a + b) / 2)
    if (CHigh * CLow > 0.0)
    {
      impliedVol <- NA
    } else{
      while (midCdiff ^ 2 > tolerance ^ 2)
      {
        midp <- (a + b) / 2
        midCdiff <-
          Cmkt - BlackScholesPricer(isCall, S, K, tau, r, midp)
        ifelse(midCdiff > 0.0, a <- midp, b <- midp)
        impliedVol <- midp
      }
    }
    impliedVol
  }
#Calculate implied volatilities using secant method ----
secantImpliedVols <-
  function(isCall, Cmkt, S, K, tau, r, vols, fn_0) {
    fn_1 <- Cmkt - BlackScholesPricer(isCall, S, K, tau, r, vols[2])
    vol3 <- vols[2] - fn_1 * (vols[2] - vols[1]) / (fn_1 - fn_0)
    ifelse (
      equalswithinTolerance(vol3, vols[2]),
      vol3,
      secantImpliedVols(isCall, Cmkt, S, K, tau, r, c(vols[2], vol3), fn_1)
    )
  }

#Calculate d1 for B-S formulation ----
getD1 <- function(S, K, tau, r, vol)
{
  d1 <- (log(S / K) + (r + 0.5 * vol ^ 2) * tau) / (vol * sqrt(tau))
}
#Calculate Delta exactly ----
delta <- function(S, K, tau, r, vol)
{
  pnorm(getD1(S, K, tau, r, vol))
}
#Calculate Delta approximately ----
CallDeltaApproximation <- function(S, K, tau, r, vol, del)
{
  (BlackScholesPricer(T, S + del, K, tau, r, vol) - BlackScholesPricer(T, S, K, tau, r, vol)) /
    del
}
#Calculate Vega Exactly ----
vega <- function(S, K, tau, r, vol) {
  ndprime(getD1(S, K, tau, r, vol)) * S * sqrt(tau)
}
#Calculate Vega approximately ----
CallVegaApproximation <- function(S, K, tau, r, vol, del)
{
  (BlackScholesPricer(T, S, K, tau, r, vol + del) - BlackScholesPricer(T, S, K, tau, r, vol)) /
    del
}
#Calculate Gamma Exactly ----
gamma <- function(S, K, tau, r, vol) {
  ndprime(getD1(S, K, tau, r, vol)) / (S * vol * sqrt(tau))
}
#Calculate Gamma approximately ----
CallGammaApproximation <- function(S, K, tau, r, vol, del)
{
  (
    BlackScholesPricer(T, S + del, K, tau, r, vol) -
      2 * BlackScholesPricer(T, S, K, tau, r, vol) +
      BlackScholesPricer(T, S - del, K, tau, r, vol)
  ) / (del ^ 2)
}
ndprime <- function(d1)
{
  exp(-d1 ^ 2 / 2) / sqrt(2 * pi)
}
#Checks for equality between two variables within 0.0001 ----
equalswithinTolerance <- function(x, y)
{
  (x - y) ^ 2 <= tolerance ^ 2
}

#Generates Implied Volatilities by 2 methods (Bisection and Secant)
# returns a dataframe of Strike, Time to Maturity, Implied Vols ----
GetImpliedVolatilities <- function(r, security, matDates, today)
{
  calls <- getOptionsData(security, matDates , today)
  Imp.Vol.Bisection <- c(0)
  Imp.Vol.Secant <- c(0)
  for (i in 1:nrow(calls))
  {
    rowq <- calls[i,]
    isCall <- rowq[, "Type"] == "Call"
    Imp.Vol.Bisection[i] <-
      calcImpliedVols(
        isCall <- rowq[, "Type"] == "Call",
        isBisection = T,
        Cmkt = rowq[, "Market.Price"],
        S = 130,
        K = rowq[, "Strike"],
        tau = rowq[, "Time"],
        r = 0.05,
        0.0001
      )
    Imp.Vol.Secant[i] <-
      calcImpliedVols(
        isCall <- rowq[, "Type"] == "Call",
        isBisection = F,
        Cmkt = rowq[, "Market.Price"],
        S = 130,
        K = rowq[, "Strike"],
        tau = rowq[, "Time"],
        r = 0.05,
        0.0001
      )
  }
  calls["Imp.Vol.Bisection"] <- Imp.Vol.Bisection
  calls["Imp.Vol.Secant"] <- Imp.Vol.Secant
  calls["Moneyness"] <-
    Moneyness(130, calls[, 1], calls[, "Time"], rep(r, nrow(calls)))
  calls <-
    subset(calls,!is.na(Imp.Vol.Bisection) &
             !is.na(Imp.Vol.Secant))
  calls[, c(1, 2, 4, 5, 6, 7)]
}
#Accepts a dataframe with volatilities, strikes and maturity date and evaluates
#Calculates Moneyness
Moneyness <- function(S, K, tau, r) {
  K * exp(-r * tau) / S
}
#Gamma, Vega and Delta ----
PopulateGreeks <- function(dataframe, S, r)
{
  df <-
    data.frame("Strike" = dataframe[, 1],
               "Time" = dataframe[, 2],
               "Imp.Vol" = dataframe[, 4])
  df[, "Delta"] <-
    delta(S, dataframe[, 1], dataframe[, 2], r, dataframe[, 4])
  df[, "Vega"] <-
    vega(S, dataframe[, 1], dataframe[, 2], r, dataframe[, 4])
  df[, "Gamma"] <-
    gamma(S, dataframe[, 1], dataframe[, 2], r, dataframe[, 4])
  df
}

#Utility to download options data to Csv from Yahoo! Finance ----
DownloadOptionsDataToCsv <- function(security, matDates)
{
  for (i in 1:length(matDates)) {
    df <-
      getOptionChain(security, src = "yahoo", Exp = matDates[i])$calls
    write.table(
      df,
      file = paste0(security, "_", matDates[i], ".csv"),
      sep = ",",
      qmethod = "double"
    )
  }
}

#Loads options data for a security from .csv file; Format: security_matdate.csv e.g. AAPL_2015-04-17.csv ----
getOptionsData <- function(security, matDate, currentDate)
{
  df <- read.csv(paste0(security, "_", matDate, ".csv"))
  dataframe <- data.frame("Strike" = df[, 1])
  timeToMaturity <-
    as.numeric(difftime(matDate, currentDate), units = "days") / 360
  dataframe["Time"] <- rep(timeToMaturity, nrow(df))
  dataframe["Type"] <- rep("Call", nrow(df))
  dataframe["Market.Price"] <- (df[, "Bid"] + df[, "Ask"]) / 2
  dataframe
}
#Calculate area under a function f between a and b with N steps using trapezoid Rule ----
TrapezoidRule <- function(f, a, b, N)
{
  h <- (b - a) / N
  theSum = 0.5 * h * (f(a) + f(a + h))
  for (i in 1:(N - 1)) {
    theSum <- theSum + 0.5 * h * (f(a + i * h) + f(a + (i + 1) * h))
  }
  theSum
}
#Calculate area under a function f between a and b with N steps using Simpson's Rule ----
SimpsonRule <- function(f, a, b, N)
{
  h <- (b - a) / N
  xa <- seq(from = a, to = b - h, by = h)
  xb <- seq(from = a + h, to = b, by = h)
  sum((xb - xa) / 6 * (f(xa) + 4 * f((xa + xb) / 2) + f(xb)))
}

#Gets Truncation Error for area calculation from Simpson Rule and Trapezoid Rule ----
TruncationError <- function(f, as, Ns) {
  a <- as[length(as)]
  n <- Ns[length(Ns)]
  TrapError <- TrapezoidRule(f, -a, a, Ns[1]) - pi
  SimpError <- SimpsonRule(f, -a, a, Ns[1]) - pi
  dfN <-
    data.frame(
      "a" = a,
      "N" = Ns[1],
      "Simpson Error" = SimpError,
      "Trapezoid Error" = TrapError
    )
  for (i in 2:length(Ns))
  {
    dfN <-
      rbind(dfN, c(
        a,
        Ns[i],
        SimpsonRule(f, -a, a, Ns[i]) - pi,
        TrapezoidRule(f, -a, a, Ns[i]) - pi
      ))
  }
  dfA <-
    data.frame(
      "a" = as[1],
      "N" = n,
      "Simpson Error" = SimpError,
      "Trapezoid Error" = TrapError
    )
  for (i in 2:length(as))
  {
    dfA <-
      rbind(dfA, c(
        as[i],
        n,
        SimpsonRule(f, -as[i], as[i], n) - pi,
        TrapezoidRule(f, -as[i], as[i], n) - pi
      ))
  }
  list(dfN, dfA)
}
#Iterative approximation of area under a function f using Simpson's Rule ----
SimpsonIterativeApproximation <- function(func, a, b, epsilon)
{
  begin <- Sys.time()
  k <- 3
  v1 <- SimpsonRule(func, a, b, 2)
  v2 <- SimpsonRule(func, a, b, 4)
  while (abs(v2 - v1) > epsilon)
  {
    v1 <- v2
    v2 <- SimpsonRule(func, a, b, k ^ 2)
    k = k + 1
  }
  c(
    "Value" = v2,
    "Iterations" = (k - 3),
    "Time" = (Sys.time() - begin)
  )
}

#Iterative approximation of area under a function f using Trapezoid Rule ----
TrapezoidIterativeApproximation <- function(func, a, b, epsilon)
{
  begin <- Sys.time()
  k <- 3
  v1 <- TrapezoidRule(func, a, b, 2)
  v2 <- TrapezoidRule(func, a, b, 4)
  while (abs(v2 - v1) > epsilon)
  {
    v1 <- v2
    v2 <- TrapezoidRule(func, a, b, k ^ 2)
    k = k + 1
  }
  c(
    "Value" = v2,
    "Iterations" = (k - 3),
    "Time" = (Sys.time() - begin)
  )
}

#Heston Model Implementation ----
Call.Heston <-
  function(theta = 0.1,
           sigma = 0.2,
           S0 = 1,
           lambda = 0,
           rho = -0.3,
           V0 = 0.1,
           r = 0,
           q = 0,
           tau = 5,
           Ks = c(0.5, 0.75, 1, 1.25, 1.5),
           Kappas = c(4, 2, 1)) {
    HestonIntegrand <- function(j, phi) {
      a <- kappa * theta
      b <- kappa + lambda
      if (j == 1)
      {
        b <- b - rho * sigma
        u <- 0.5
      }
      else
      {
        u <- -0.5
      }
      
      d <-
        sqrt((rho * sigma * phi * 1i - b) ^ 2 - (sigma ^ 2) * (2 * u * phi * 1i - phi ^ 2))
      g <-
        (b - rho * sigma * phi * 1i + d) / (b - rho * sigma * phi * 1i - d)
      C <-
        (r - q) * phi * 1i * tau + (kappa * theta / (sigma ^ 2)) * ((b - rho * sigma * phi * 1i + d) * tau - 2 * log((1 - g * exp(d * tau)) / (1 - g)))
      D <-
        ((b - rho * sigma * phi * 1i + d) / (sigma ^ 2)) * ((1 - exp(d * tau)) / (1 - g * exp(d * tau)))
      psi <- exp(C + D * V0 + 1i * phi * log(S0 * exp(r * tau)))
      Re((exp(-1i * phi * log(K)) * psi) / (1i * phi))
    }
    P1Function <- function(psi)
    {
      HestonIntegrand(1, psi)
    }
    P2Function <- function(psi)
    {
      HestonIntegrand(2, psi)
    }
    mat <- matrix(nrow = length(Ks), ncol = length(Kappas))
    cNames <- c(0)
    for (m in 1:length(Ks))
    {
      K <- Ks[m]
      rowM <- c(0)
      for (n in 1:length(Kappas))
      {
        kappa <- Kappas[n]
        P1 <-
          0.5 + SimpsonIterativeApproximation(P1Function, 10 ^ -20, 743, 10 ^ (-4))["Value"] / pi
        P2 <-
          0.5 + SimpsonIterativeApproximation(P2Function, 10 ^ -20, 743, 10 ^ (-4))["Value"] / pi
        C <- S0 * P1 - K * exp(-(r - q) * tau) * P2
        mat[m, n] <- C
      }
      colnames(mat) <- paste("Kappa =", Kappas)
      rownames(mat) <- paste("K =", Ks)
    }
    mat
  }
#Utility to check the deviation of one value from another as a percentage ----
deviationFrom <- function(a, b) {
  abs(b - a) * 100 / a
}
