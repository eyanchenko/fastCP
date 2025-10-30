# Naive implementation where we recalculate the entire objective function each time.

cppFunction('IntegerVector borgatti_naiveCpp(IntegerMatrix A){

  double n = A.nrow();
  double m = n*(n-1)/2;
  int k = 0;
  
  
  // Initialize vector 
  
  IntegerVector C(n);
  
  
  for(int i = 0; i < n; i++){
    
    if(rand()%10 > 5){
      C(i) = 1;
      k++;
    }else{
      C(i) = 0;
    }
     
  }
  
  int kk = k;
  
  
  double sum_A = 0; 
  double sum_CP = 0.5*k*(k-1)+(n-k)*k; 
  double sum_ACP = 0;
  
  for (int i = 1; i < n; i++){
      for(int j = 0; j < i; j++){
        sum_A   += A(i,j);
        sum_ACP += A(i,j) * (C(i) + C(j) - C(i)*C(j)); 
      }
    }
  
  double obj_last = (sum_ACP - sum_A / m * sum_CP) / 
    (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP - sum_CP/m*sum_CP) );
  
                   
  int ind = 1;
  int n_iter = 0;
  int max_iter = 100;
  
  IntegerVector Ctest = clone(C);
  
  IntegerVector v(n);
  for(int i=0; i<n;i++){v(i) = i;}
  
  double obj_new = 0;
  
  while(ind > 0 & n_iter < max_iter){
  
  ind = 0;
  /* Randomly shuffle node order */
  int N = n;
  for(int a=0;a<n;a++){
    int b = round(rand()%N);
    int c = v(a); v(a)=v(b); v(b)=c;
  }
  
  
  
   for(int i = 0; i < n; i++){
    
     Ctest = clone(C);
     
     Ctest(v(i)) = 1 - Ctest(v(i));
     
    /* Update size of core for test vector, kk */
    
     if(Ctest(v(i))==0){
       kk = k - 1;
     }else{
       kk = k + 1;
     }


     /* (Naively) update value of objective function */
     double sum_CP_test = 0.5*kk*(kk-1) + (n-kk)*kk; 
     double sum_ACP_test = 0;
     
     for (int i = 1; i < n; i++){
      for(int j = 0; j < i; j++){
        sum_ACP_test += A(i,j) * (Ctest(i) + Ctest(j) - Ctest(i)*Ctest(j)); 
      }
     }
  
     obj_new = (sum_ACP_test - sum_A / m * sum_CP_test) / 
      (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP_test - sum_CP_test/m*sum_CP_test) );

    
     if(obj_new > obj_last){
       C(v(i)) = Ctest(v(i));
       sum_ACP = sum_ACP_test;
       ind++;
       obj_last = obj_new;
       k = kk;
      
       
     }
    
    
    }
    
    n_iter++;
    
  }
  
  return( C );
}')


# Proposed implementation
cppFunction('IntegerVector borgattiCpp(IntegerMatrix A){

  double n = A.nrow();
  double m = n*(n-1)/2;
  int k = 0;
  
  
  // Initialize vector 
  
  IntegerVector C(n);
  
  
  for(int i = 0; i < n; i++){
    
    if(rand()%10 > 5){
      C(i) = 1;
      k++;
    }else{
      C(i) = 0;
    }
     
  }
  
  int kk = k;
  
  
  double sum_A = 0; 
  double sum_CP = 0.5*k*(k-1)+(n-k)*k; 
  double sum_ACP = 0;
  
  for (int i = 1; i < n; i++){
      for(int j = 0; j < i; j++){
        sum_A   += A(i,j);
        sum_ACP += A(i,j) * (C(i) + C(j) - C(i)*C(j)); 
      }
    }
  
  double obj_last = (sum_ACP - sum_A / m * sum_CP) / 
    (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP - sum_CP/m*sum_CP) );
  
                   
  int ind = 1;
  int n_iter = 0;
  int max_iter = 100;
  
  IntegerVector Ctest = clone(C);
  
  IntegerVector v(n);
  for(int i=0; i<n;i++){v(i) = i;}
  
  double obj_new = 0;
  
  while(ind > 0 & n_iter < max_iter){
  
  ind = 0;
  /* Randomly shuffle node order */
  int N = n;
  for(int a=0;a<n;a++){
    int b = round(rand()%N);
    int c = v(a); v(a)=v(b); v(b)=c;
  }
  
  
  
   for(int i = 0; i < n; i++){
    
     Ctest = clone(C);
     
     Ctest(v(i)) = 1 - Ctest(v(i));
     
    /* Update size of core for test vector, kk */
    
     if(Ctest(v(i))==0){
       kk = k - 1;
     }else{
       kk = k + 1;
     }

     double sum_CP_test = 0.5*kk*(kk-1) + (n-kk)*kk; 
     double sum_ACP_test = sum_ACP;
    
    /* Update value of objective function */
    
     for(int j = 0; j < n; j++){
        sum_ACP_test += A(v(i),j) * (Ctest(v(i)) + Ctest(j) - Ctest(v(i))*Ctest(j)); 
        sum_ACP_test -= A(v(i),j) * (C(v(i)) + C(j) - C(v(i))*C(j)); 
     }

    
     obj_new = (sum_ACP_test - sum_A / m * sum_CP_test) / 
      (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP_test - sum_CP_test/m*sum_CP_test) );

    
     if(obj_new > obj_last){
       C(v(i)) = Ctest(v(i));
       sum_ACP = sum_ACP_test;
       ind++;
       obj_last = obj_new;
       k = kk;
      
       
     }
    
    
    }
    
    n_iter++;
    
  }
  
  return( C );
}')


### Plot adjacency matrices
plot_matrix = function(in_mat){
  require(ggplot2)
  
  x = seq(1, nrow(in_mat))
  y = seq(1, ncol(in_mat))
  
  df = expand.grid(X=x, Y=y)
  df$Z = c(in_mat)
  
  # doesn't seem to care about those colors at all, but it works
  ggplot(df, aes(X, Y, fill= Z)) + 
    geom_tile()+
    scale_fill_continuous(type = "viridis",
                          limits = c(0,1), 
                          breaks = c(0, 1/4, 2/4, 3/4, 4/4),
                          guide_colourbar(nbin = 100),
                          name = "theta")+
    theme_minimal()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    theme(legend.position = "none")
  
  
}

#' @title Borgatti and Everett core-periphery metric
#' @description This function evaluates the BE CP metric
#' @param A adjacency matrix
#' @param C CP labels (C[i] = 1 if node i in core, and 0 otherwise)
#' @return value of BE metric of A at C
#' @export
obj.fun <- function(A,C){
  
  if(nrow(A)!=ncol(A)){
    stop("A must be a square matrix.")
  }
  
  if(nrow(A)!=length(C)){
    stop("Number of rows/columns in A must be the same as number of elements in C.")
  }
  
  if( sum(sort(unique(C))==c(0,1))!=2 ){
    stop("C must only contain 0's and 1's.")
  }
  
  if(sum(C)==length(C) || sum(C)==0){
    stop("C must contain at least one 1 and one zero.")
  }
  
  n = length(C)
  CP <- C #idealized CP matrix to compare
  
  CP[C==1] = 2
  CP[C==0] = 1
  
  CP = CP%*%t(CP)%%2
  CP =  1 - CP
  diag(CP) <- 0
  
  #don't allow partitions of all core or all periphery
  if( sum(CP)==0 || sum(CP)==n*(n-1) ){
    rho=0
  }else{
    rho = cor(A[upper.tri(A)], CP[upper.tri(CP)])
  }
  
  return(rho)
}

#' @title Model-agnostic core-periphery population parameter
#' @description This function evaluates the CP population parameter
#' @param P data-generating probability matrix
#' @param C CP labels (C[i] = 1 if node i in core, and 0 otherwise)
#' @return value of CP populations parameter of P at C
#' @export
rho.fun <- function(P, C){
  if(nrow(P)!=ncol(P)){
    stop("P must be a square matrix.")
  }
  
  if(nrow(P)!=length(C)){
    stop("Number of rows/columns in P must be the same as number of elements in C.")
  }
  
  if( sum(sort(unique(C))==c(0,1))!=2 ){
    stop("C must only contain 0's and 1's.")
  }
  
  if(sum(C)==length(C) || sum(C)==0){
    stop("C must contain at least one 1 and one zero.")
  }
  
  n = length(C)
  k = sum(C)
  
  CP <- C #idealized CP matrix to compare
  
  CP[C==1] = 2
  CP[C==0] = 1
  
  CP = CP%*%t(CP)%%2
  CP =  1 - CP
  diag(CP) <- 0
  
  p     = mean(P[upper.tri(P)])
  d =  (0.5*k*(k-1)+(n-k)*k) / choose(n,2)
  
  
  rho = (0.5*sum(P*CP) - 0.5*n*(n-1)*p*d) / (0.5*n*(n-1)*sqrt(p*(1-p)*d*(1-d)))
  
  return(rho)
}

generateP <- function(n, p11, p12, p22, prop=0.50){
  
  P = matrix(p22, ncol=n, nrow=n)
  k = round(n*prop)
  
  P[1:k, 1:k] <- p11
  P[1:k, (k+1):n] <- P[(k+1):n, 1:k] <- p12
  diag(P) <- 0
  
  return(P)
}
generatePDCBM <- function(theta, p11, p12, p22, prop=0.1){
  n = length(theta)
  P = matrix(p22, ncol=n, nrow=n)
  k = round(n*prop)
  
  P[1:k, 1:k] <- p11
  P[1:k, (k+1):n] <- P[(k+1):n, 1:k] <- p12
  
  P <- theta%*%t(theta) * P
  
  diag(P) <- 0
  
  P[P > 1] <- 1
  
  return(P)
}
generateA <- function(n, p11, p12, p22, prop=0.50){
  m = as.integer(round(n^2*prop^2, digits=0))
  A11 = matrix(rbinom(m, 1, p11), ncol=as.integer(round(n*prop)))
  A11[lower.tri(A11)] <- 0
  diag(A11) <- 0
  A11 = A11 + t(A11)
  
  m = as.integer(round(n^2*(1-prop)^2))
  A22 = matrix(rbinom(m, 1, p22), ncol=as.integer(round(n*(1-prop))), nrow = as.integer(round(n*(1-prop))))
  A22[lower.tri(A22)] <- 0
  diag(A22) <- 0
  A22 = A22 + t(A22)
  
  m = as.integer(round(n^2*prop*(1-prop)))
  A12 <- matrix(rbinom(m, 1, p12), ncol=as.integer(round(n*(1-prop))), nrow = as.integer(round(n*prop)))
  A21 <- t(A12)
  
  A = cbind(rbind(A11, A21), rbind(A12, A22))
  return(A)
}
generateDCBM <- function(theta, p11, p12, p22, q=0.1){
  #theta: different degrees
  #pij: block-block probabilities
  #q: proportion of nodes in block 1
  
  n <- length(theta)
  
  Omega11 <- matrix(p11, ncol = n*q, nrow = n*q)
  Omega22 <- matrix(p22, ncol = n*(1-q), nrow = n*(1-q))
  Omega12 <- matrix(p12, ncol = n*(1-q), nrow = n*q)
  
  Omega <- rbind(cbind(Omega11, Omega12), cbind(t(Omega12), Omega22) )
  
  P = as.vector(theta %*% t(theta) * Omega)
  
  P[P > 1] <- 1
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  
  A[lower.tri(A, diag=TRUE)] <- 0
  A = A + t(A)
  
  return(A) 
  
}





# Classification accuracy
class_acc <- function(Cpred, Ctrue){
  sum(Cpred == Ctrue) / length(Cpred)
}

#source_python('~/Documents/Research/Srijan/Cp/cp_be.py')

C_be <- function(A){
  unlist(cp_be(A))
}














