library(VCM)

# standardize phenotype (no covariates)
stand_pheno = function(y, P){
  # y: {0, 1}
  z = (y - P)/sqrt(P*(1 - P))
  return(z)
}

# standardize genotype
stand_geno = function(G){
  X = scale(G)
  return(X)
}

# liability threshold conditional on covariates
get_ti = function(t, Z, alpha){
  # Z: covariates (n*c matrix)
  # beta: fixed effect size (c-vector)
  ti = drop(t - Z%*%alpha)
  return(ti)
}

# population prevalence conditional on covariates
get_Ki = function(ti){
  Ki = 1 - pnorm(ti)
  return(Ki)
}

# in-sample prevalence conditional on covariates
get_Pi = function(Ki){
  inclu_prob = Ki + (1 - Ki)*K*(1 - P)/P/(1 - K)
  # under full ascertainment
  Pi = Ki/inclu_prob
  return(Pi)
}

# no covariates
pcgc_reg = function(y, X, K, P, t){
  
  # y: standardized phenotype (mean 0, var 1)
  # X: standardized genotype matrix (column mean 0, var 1)
  # K: population prevalence
  # P: in-sample prevalence
  # t: liability threshold
  
  n = nrow(X)
  p = ncol(X)
  
  # phenotype correlations
  pheno_corr = outer(y, y) # n*n matrix, each entry is the phenotype correlation between two individuals
  
  # genotype correlations
  geno_corr = X%*%t(X)/p # n*n matrix, each entry is the genetic correlation between two individuals
  
  # multiply genotype correlations with the constant c
  phi_t = dnorm(t)
  c = P*(1 - P)*phi_t^2/(K^2*(1 - K)^2)
  geno_corr_c = c*geno_corr
  
  # regress phenotype correlations on genotype correlations
  pheno_corr_vec = pheno_corr[upper.tri(pheno_corr)]
  geno_corr_c_vec = geno_corr_c[upper.tri(geno_corr_c)]
  
  # hlsq = sum(geno_corr_c_vec*pheno_corr_vec)/sum(geno_corr_c_vec^2)
  # resi = pheno_corr_vec - hlsq*geno_corr_c_vec
  # hlsq_se = sqrt(sum(resi^2)/length(resi)/sum(geno_corr_c_vec^2))
  # res = list(hlsq=hlsq, hlsq_se=hlsq_se)
  # return(res)
  
  data = data.frame(y=pheno_corr_vec, X=geno_corr_c_vec)
  fit = lm(y~X, data=data)
  hlsq = summary(fit)$coefficients[2,1] # liability scale heritability
  hlsq_se = summary(fit)$coefficients[2,2] # standard error of hlsq
  
  res = list(hlsq=hlsq, hlsq_se=hlsq_se)
  return(res)
  
}

# with covariates
pcgc_reg_covar = function(y, X, K, P, t, Ki, Pi, ti){
  
  # y: standardized phenotype (mean 0, var 1)
  # X: standardized genotype (col mean 0, var 1)
  # K: population prevalence
  # P: in-sample prevalence
  # t: liability threshold
  # Ki: population prevalence conditional on covariates of individual i (n-vector)
  # Pi: in-sample prevalence conditional on covariates of individual i (n-vector)
  # ti: liability threshold conditional on covariates of individual i (n-vector)
  
  n = nrow(X)
  p = ncol(X)
  
  # phenotype correlations
  pheno_corr = outer(y, y) # n*n matrix, each entry is the phenotype correlation between two individuals
  
  # genotype correlations
  geno_corr = X%*%t(X)/p # n*n matrix, each entry is the genetic correlation between two individuals
  
  # multiply genotype correlations with the constant ci
  phi_ti = dnorm(ti) # n-vector
  c1 = (P-K)/P/(1-K)
  c2 = K*(1-P)/P/(1-K)
  c_vec_numer = phi_ti*(1 + c1*Pi)
  c_vec_denom = sqrt(Pi*(1 - Pi))*(Ki + c2*(1 - Ki))
  c_vec = c_vec_numer/c_vec_denom # n-vector
  c_mat = outer(c_vec, c_vec) # n*n matrix
  geno_corr_c = c_mat*geno_corr
  
  # regress phenotype correlations on genotype correlations
  pheno_corr_vec = pheno_corr[upper.tri(pheno_corr)]
  geno_corr_c_vec = geno_corr_c[upper.tri(geno_corr_c)]
  
  data = data.frame(y=pheno_corr_vec, X=geno_corr_c_vec)
  fit = lm(y~., data=data)
  hlsq = summary(fit)$coefficients[2,1] # liability scale heritability
  hlsq_se = summary(fit)$coefficients[2,2] # standard error of hlsq
  
  res = list(hlsq=hlsq, hlsq_se=hlsq_se)
  return(res)
  
}

# simulate case-control data using dynamic programming
sim_dym = function(n, p, K, P, hlsq, beta){
  
  # n: number of samples
  # p: number of SNPs
  # K: population prevalence
  # P: in-sample prevalence
  # hlsq: true liability scale heritability
  # beta: random effect size
  
  n_case = n*P
  n_control = n*(1 - P)
  t = qnorm(1 - K)
  
  mafs = runif(p, 0.05, 0.5)
  
  n1 = n0 = 0
  G1 = X1 = matrix(NA, n_case, p)
  G0 = X0 = matrix(NA, n_control, p)
  l1 = l0 = c()
  ## simulate cases
  # break the joint probability into a serial product of conditional probabilities
  # P(g1,...,gm|p=1) = P(g1|p=1)\prod_2^m{P(gi|g1,...,gi-1,p=1)}
  # P(gi|g1,...,gi-1,p=1) = P(p=1|g1,...,gi)P(gi|g1,...,gi-1)/P(p=1|g1,...,gi-1)
  # under linkage equilibrium, P(gi|g1,...,gi-1)=P(gi)
  # P(p=1|g1,...,gi) = 1 - pnorm((t - \sum_1^i{uj*gj})/sqrt(1 - \sum_1^i{uj^2}))
  # column-by-column
  # for the first SNP (not sure)
  G1[,1] = rbinom(n_case, 2, mafs[1])
  X1[,1] = (G1[,1] - 2*mafs[1])/sqrt(2*mafs[1]*(1 - mafs[1]))
  # for the second to end SNPs
  for(i in 2:p){
    tmp1_numer = t - drop(X1[,1:i]%*%beta[1:i])
    tmp2_numer = sqrt(1 - sum(beta[1:i]^2))
    prob_numer = 1 - pnorm(tmp1_numer/tmp2_numer)
    tmp1_denom = t - drop(X1[,1:(i-1)]%*%beta[1:(i-1)])
    tmp2_denom = sqrt(1 - sum(beta[1:(i-1)]^2))
    prob_denom = 1 - pnorm(tmp1_denom/tmp2_denom)
    prob_i = prob_numer*mafs[i]/prob_denom
    G1[,i] = rbinom(n_case, 2, prob_i)
    X1[,i] = scale(G1[,i])
  }
  ## simulate controls
  G0[,1] = rbinom(n_control, 2, mafs[1])
  X0[,1] = (G0[,1] - 2*mafs[1])/sqrt(2*mafs[1]*(1 - mafs[1]))
  # for the second to end SNPs
  for(i in 2:p){
    tmp1_numer = t - drop(X0[,1:i]%*%beta[1:i])
    tmp2_numer = sqrt(1 - sum(beta[1:i]^2))
    prob_numer = 1 - pnorm(tmp1_numer/tmp2_numer)
    tmp1_denom = t - drop(X1[,1:(i-1)]%*%beta[1:(i-1)])
    tmp2_denom = sqrt(1 - sum(beta[1:(i-1)]^2))
    prob_denom = pnorm(tmp1_denom/tmp2_denom)
    prob_i = prob_numer*mafs[i]/prob_denom
    G0[,i] = rbinom(n_control, 2, prob_i)
    X0[,i] = scale(G1[,i])
  }
  G = rbind(G1, G0)
  X = rbind(X1, X0)
  
  y = c(rep(1, n_case), rep(0, n_control))
  y_std = stand_pheno(y, P)
  
  sim_data = list(y=y, y_std=y_std, G=G, X=X)
  return(sim_data)
  
}
