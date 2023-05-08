unbiased_Quantile_normal <- function(y,B,alpha){

   options(warn=-1)
   #Input:
   #B - anzahl der bootstaps
   #y - beobachtete Stichprobe
   #alpha - theoretische quantil
   #
   #Output: alpha_pu = probability unbiased alpha

   #1.Schritt: Schätzen des parametervectors via MLE
   out <- optim(c(0.2,0.2), normal.lik, y= y, method="L-BFGS-B",lower=c(-Inf,.001), control=list(maxit=10000))
   theta <- out$par

   #2.Schritt: Ziehe B samples um umfang der Stichprobe mit der plug-in Vtl

   #erzeuge eine leere Matrix (Bx2) für schätzwerte
   est_matrix <-  matrix(0, nrow=B, ncol=2)
   n <- length(y)
   for(i in 1:B){
     y_hat <- rnorm(n, mean= theta[1], sd=theta[2])
     #berechne neue Schätzwerte
     out <- optim(c(0.2,0.2), normal.lik, y= y_hat, method="L-BFGS-B",lower=c(-Inf,.001), control=list(maxit=10000))
     est_matrix[i,] <- out$par
   }

   #3.Schritt: Minimiere den Abstand

   #out1 <- optim(alpha , min_quant_norm, method="L-BFGS-B",lower=0.01,upper=.9, control=list(maxit=10000), est_m=est_matrix, alpha=alpha, theta=theta)
   out1 <- optim(alpha , min_quant_norm, method="BFGS", control=list(maxit=10000), est_m=est_matrix, alpha=alpha, theta=theta)

   #output <- list(x = out1$par, y = out1$convergence,  w= out1$counts, z = out1$message)
   output <- list(x = out1$par, y = out1$convergence)

   return(output)
 }


 ####################################################################################
 normal.lik <- function(theta, y){
   mu <- theta[1]
   sigma <- theta[2]
   n <- nrow(y)
   logl<- sum(dnorm(y, mu, sigma, log = TRUE))
   #maximum likelihood, dazu minimieren wir minus die funktion
   return(-logl)
 }
 ####################################################################################
 min_quant_norm <- function(gamma, est_m,alpha,theta){


   inv_Vtl   <- qnorm(gamma, mean = est_m[,1], sd = est_m[,2])
   Vtl       <- pnorm(inv_Vtl, mean=theta[1], sd=theta[2])
   abs_value <- abs(sum(Vtl)/length(Vtl) - alpha)
   return(abs_value)
 }