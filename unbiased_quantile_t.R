unbiased_Quantile_t<- function(y,B,alpha){

options(warn=-1)

   #Input:
   #B - anzahl der bootstaps
   #y - beobachtete Stichprobe
   #alpha - theoretische quantil
   #
   #Output: alpha_pu = probability unbiased quantile

   #1.Schritt: Schätzen des parametervectors via MLE

 
   tmp <- optim(c(0.2,0.2,1), t.lik, y=y, method="BFGS")
   tmp <- tmp$par
   theta <- c(tmp[3], tmp[1], tmp[2])
   #2.Schritt: Ziehe B samples um umfang der Stichprobe mit der plug-in Vtl

   #erzeuge eine leere Matrix (Bx2) für schätzwerte
   est_matrix <-  matrix(0, nrow=B, ncol=3)
   n <- length(y)
   for(i in 1:B){
     y_hat <- theta[[2]] + theta[[3]]*rt(n, df=theta[[1]])
     #berechne neue Schätzwerte(theta*)
    
     out <- optim(c(0.2,0.2,1), t.lik, y=y_hat, method="BFGS")
     out <- out$par
     out_order <- c(out[3], out[1], out[2])
     est_matrix[i,] <- out_order
    
   }

   #3.Schritt: Minimiere den Abstand

   out1 <- optim(alpha , min_quant_t, method="BFGS", control=list(maxit=10000),est_m=est_matrix, alpha=alpha, theta=theta)

   output <- list(x = out1$par, y = out1$convergence)

   return(output)
 }
 ######################################################################
 ####################################################################################
 min_quant_t <- function(gamma, est_m,alpha,theta){
   abs_value <- 0;
   for(i in 1:dim(est_m)[1]){
   inv_Vtl   <- est_m[i,2]+ est_m[i,3]*qt(gamma, df=est_m[i,1])
   Vtl       <- pt((inv_Vtl-theta[[2]])/theta[[3]],df=theta[[1]])
   abs_value <- abs_value + Vtl
   }

   abs_value <- abs(1/dim(est_m)[1]*(abs_value)- alpha)
   return(abs_value)
 }

#############################################################################################
t.lik <-function(y,par){
   if(par[2]>0 & par[3]>0) return(-sum(log(dt((y-par[1])/par[2],df=par[3])/par[2])))
   else return(Inf)
 }
