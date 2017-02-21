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