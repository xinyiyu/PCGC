library(VCM)
library(ggplot2)
### implement pcgc-regression

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
  
  n = nrow(X)
  p = ncol(X)
  t = qnorm(1 - K)
  
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

### generate genotypes and phenotypes
## use generative model (for scenarios with a small number of SNPs, e.g 100 or 1,000 SNPs)
p = 1000 # number of SNPs
mafs = runif(p, 0.05, 0.5)
hlsq = 0.5 # true liability scale heritability
beta = rnorm(p, 0, sqrt(hlsq/p)) # SNP effect sizes
K = 0.1 # population prevalence
t = qnorm(1 - K) # liability threshold
P = 0.15 # in-sample prevalence

# full ascertainment (if the generated individual is a case, then she will be included in the study)
n = 4000
n_case = n*P
n_control = n*(1 - P)

c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2) # correction constant between hlsq and hosq

nrep = 30
EST_HLSQ = data.frame(matrix(NA, nrep, 3))
names(EST_HLSQ) = c('rep', 'pcgc', 'reml')
EST_HLSQ$rep = seq(nrep)
for(rep in 1:nrep){
  ## generate data
  n1 = n0 = 0
  G1 = X1 = matrix(NA, n_case, p)
  G0 = X0 = matrix(NA, n_control, p)
  l1 = l0 = c()
  # no covariates
  while(n1 < n_case | n0 < n_control){
    # generate a genotype
    g = rbinom(p, 2, mafs)
    # standarized the genotype
    x = (g - 2*mafs)/sqrt(2*mafs*(1 - mafs))
    # liability for this individual
    l = sum(x*beta) + rnorm(1, 0, sqrt(1-hlsq))
    # phenotype
    if(l > t & n1 < n_case){
      n1 = n1 + 1
      G1[n1,] = g
      X1[n1,] = x
      l1[n1] = l
    }else{
      if(n0 < n_control){
        n0 = n0 + 1
        G0[n0,] = g
        X0[n0,] = x
        l0[n0] = l
      }
    }
  }
  G = rbind(G1, G0)
  X = rbind(X1, X0)
  l = c(l1, l0)
  y = c(rep(1, n_case), rep(0, n_control))
  y_std = stand_pheno(y, P)
  
  message(sprintf('Finish generating data: %s-th rep.', rep))
  
  ## pcgc regression
  res_pcgc = pcgc_reg(y_std, X, K, P, t)
  EST_HLSQ[EST_HLSQ$rep==rep,]$pcgc = res_pcgc$hlsq
  message(sprintf('Finish pcgc regression: %s-th rep.', rep))
  
  ## reml
  res_reml = linRegPXEM(X=X, y=y, verbose=T)
  hosq_reml = res_reml$sb2/(res_reml$sb2 + res_reml$se2)
  EST_HLSQ[EST_HLSQ$rep==rep,]$reml = c*hosq_reml
  message(sprintf('Finish reml estimation: %s-th rep.', rep))
  
}

## boxplot
PLOT_EST_HLSQ = data.frame(matrix(NA, 2*nrep, 3))
names(PLOT_EST_HLSQ) = c('rep', 'method', 'est_hlsq')
PLOT_EST_HLSQ$rep = rep(seq(nrep), each=2)
PLOT_EST_HLSQ$method = rep(c('pcgc', 'reml'), nrep)
PLOT_EST_HLSQ[PLOT_EST_HLSQ$method=='pcgc',]$est_hlsq = EST_HLSQ$pcgc
PLOT_EST_HLSQ[PLOT_EST_HLSQ$method=='reml',]$est_hlsq = EST_HLSQ$reml

g_est = ggplot(data=PLOT_EST_HLSQ, aes(x=method, y=est_hlsq, color=method)) +
  geom_boxplot() + geom_hline(yintercept=hlsq, lty='dashed') +
  xlab('Method') + ylab(paste('Estimated ', expression(h[l]^2))) +
  ggtitle('Heritability Estimation: PCGC vs REML') +
  theme(plot.title=element_text(hjust=0.5))
g_est

### test
n1 = n0 = 0
G1 = X1 = matrix(NA, n_case, p)
G0 = X0 = matrix(NA, n_control, p)
l1 = l0 = c()
# no covariates
while(n1 < n_case | n0 < n_control){
  # generate a genotype
  g = rbinom(p, 2, mafs)
  # standarized the genotype
  x = (g - 2*mafs)/sqrt(2*mafs*(1 - mafs))
  # liability for this individual
  l = sum(x*beta) + rnorm(1, 0, sqrt(1-hlsq))
  # phenotype
  if(l > t & n1 < n_case){
    n1 = n1 + 1
    G1[n1,] = g
    X1[n1,] = x
    l1[n1] = l
  }else{
    if(n0 < n_control){
      n0 = n0 + 1
      G0[n0,] = g
      X0[n0,] = x
      l0[n0] = l
    }
  }
}
G = rbind(G1, G0)
X = rbind(X1, X0)
l = c(l1, l0)
y = c(rep(1, n_case), rep(0, n_control))
y_std = stand_pheno(y, P)

## dynamic programming (for scenarios with more SNPs, e.g 10,000 SNPs)

### estimate heritability using pcgc
res_pcgc = pcgc_reg(y_std, X, K, P, t)
res_pcgc

### estimate heritability using reml
res_reml = linRegPXEM(X=X, y=y, verbose=T)
hosq_reml = res_reml$sb2/(res_reml$sb2 + res_reml$se2)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)
hlsq_reml = c*hosq_reml
hlsq_reml


### induced G-E interaction: simulation 1
library(mvtnorm)
K = 0.01
P = 0.5
h2l = 0.9
n = 4000
n_case = n*P
n_control = n - n_case
N = 1000000 # a large pool
GE = rmvnorm(N, c(0,0), diag(c(h2l, 1-h2l)))
# plot(x=GE[,1], y=GE[,2], asp=1)
t = qnorm(1 - K)
idx_all_case = which((GE[,1]+GE[,2])>t) # all affected individuals in the pool
idx_all_control = seq(N)[-idx_all_case] # all unaffected individuals in the pool
length(idx_all_case)
# sample cases and controls
idx_study_case = sample(idx_all_case, n_case) # cases in the study
idx_study_control = sample(idx_all_control, n_control) # controls in the study
GE_study = GE[c(idx_study_case,idx_study_control),]
dim(GE_study)
# correlation of g and e in the study
cor(GE_study[,1], GE_study[,2])


h2ls = c(0.1, 0.3, 0.5, 0.7, 0.9)
Ks = c(0.01, 0.05, 0.1)
Ps = c(0.1, 0.3, 0.5)
nrep = 10

corr_g_e = data.frame(matrix(NA, 5*3*3*nrep, 5))
names(corr_g_e) = c('h2l', 'K', 'P', 'rep', 'corr')
corr_g_e$h2l = rep(h2ls, each=3*3*nrep)
corr_g_e$K = rep(rep(Ks, each=3*nrep), 5)
corr_g_e$P = rep(rep(Ps, each=nrep), 5*3)
corr_g_e$rep = rep(seq(nrep), 5*3*3)
for(h2l in h2ls){
  for(K in Ks){
    for(P in Ps){
      GE = rmvnorm(N, c(0,0), diag(c(h2l, 1-h2l)))
      t = qnorm(1 - K)
      idx_all_case = which((GE[,1]+GE[,2])>t) # all affected individuals in the pool
      idx_all_control = seq(N)[-idx_all_case] # all unaffected individuals in the pool
      n_case = n*P
      n_control = n - n_case
      for(rep in 1:nrep){
        idx_study_case = sample(idx_all_case, n_case) # cases in the study
        idx_study_control = sample(idx_all_control, n_control) # controls in the study
        GE_study = GE[c(idx_study_case,idx_study_control),]
        corr_g_e[corr_g_e$h2l==h2l&corr_g_e$K==K&corr_g_e$P==P&corr_g_e$rep==rep,]$corr = cor(GE_study[,1], GE_study[,2])
      
        # message(sprintf('h2l = %s, K = %s, P = %s, rep = %s', h2l, K, P, rep))
      }
    }
  }
}

mean_corr_g_e = data.frame(matrix(NA, 5*3*3, 4))
names(mean_corr_g_e) = c('h2l', 'K', 'P', 'corr')
mean_corr_g_e$h2l = rep(h2ls, each=3*3)
mean_corr_g_e$K = rep(rep(Ks, each=3), 5)
mean_corr_g_e$P = rep(Ps, 5*3)
for(h2l in h2ls){
  for(K in Ks){
    for(P in Ps){
      mean_corr_g_e[mean_corr_g_e$h2l==h2l&mean_corr_g_e$K==K&mean_corr_g_e$P==P,]$corr =
        mean(corr_g_e[corr_g_e$h2l==h2l&corr_g_e$K==K&corr_g_e$P==P,]$corr)
    }
  }
}

P_corr = ggplot(data=corr_g_e, aes(x=factor(h2l), y=corr, color=factor(h2l))) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(rows=vars(P), cols=vars(K)) +
  xlab(expression(h[l]^2)) + ylab('Corr(g, e)') +
  ggtitle('Induced G-E Interactions') +
  theme(plot.title=element_text(hjust=0.5))
P_corr
ggsave(filename='/Users/yxy/PCGC/boxplot_G_E_interactions_sim1.png', P_corr)

### induced G-E interaction: simulation 2
library(mvtnorm)
K = 0.01
P = 0.3
h2l = 0.5

N = 1000000 # a large pool
GE = rmvnorm(N, c(0,0), diag(c(h2l, 1-h2l)))
# plot(x=GE[,1], y=GE[,2], asp=1)
t = qnorm(1 - K)
idx_all_case = which((GE[,1]+GE[,2])>t) # all affected individuals in the pool
idx_all_control = seq(N)[-idx_all_case] # all unaffected individuals in the pool
length(idx_all_case)

# sample cases and controls
idx_study_case = sample(idx_all_case, n_case) # cases in the study
idx_study_control = sample(idx_all_control, n_control) # controls in the study
GE_study = GE[c(idx_study_case,idx_study_control),]
dim(GE_study)
# correlation of g and e in the study
cor(GE_study[,1], GE_study[,2])


ns = c(4000, 6000, 8000)
nrep = 30

corr_g_e = data.frame(matrix(NA, 3*nrep, 3))
names(corr_g_e) = c('n', 'rep', 'corr')
corr_g_e$n = rep(ns, each=nrep)
corr_g_e$rep = rep(seq(nrep), 3)
for(n in ns){
  GE = rmvnorm(N, c(0,0), diag(c(h2l, 1-h2l)))
  t = qnorm(1 - K)
  idx_all_case = which((GE[,1]+GE[,2])>t) # all affected individuals in the pool
  idx_all_control = seq(N)[-idx_all_case] # all unaffected individuals in the pool
  n_case = n*P
  n_control = n - n_case
  for(rep in 1:nrep){
    idx_study_case = sample(idx_all_case, n_case) # cases in the study
    idx_study_control = sample(idx_all_control, n_control) # controls in the study
    GE_study = GE[c(idx_study_case,idx_study_control),]
    corr_g_e[corr_g_e$n==n&corr_g_e$rep==rep,]$corr = cor(GE_study[,1], GE_study[,2])
    
    # message(sprintf('h2l = %s, K = %s, P = %s, rep = %s', h2l, K, P, rep))
  }
}

P_corr = ggplot(data=corr_g_e, aes(x=factor(n), y=corr, color=factor(n))) +
  geom_boxplot(outlier.shape=NA) +
  xlab('Sample size') + ylab('Corr(g, e)') +
  ggtitle('Induced G-E Interactions') +
  theme(plot.title=element_text(hjust=0.5))
P_corr
ggsave(filename='/Users/yxy/PCGC/boxplot_G_E_interactions_sim2.png', P_corr)

