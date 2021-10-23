library(mvtnorm)
library(ggplot2)

### code for PCGC simulation

pcgc_reg = function(y, X, K, P){
  
  # y: standardized phenotype (mean 0, var 1)
  # X: standardized genotype matrix (column mean 0, var 1)
  # K: population prevalence
  # P: in-sample prevalence
  
  t = qnorm(1 - K)
  
  n = nrow(X)
  
  # phenotype correlations
  pheno_corr = outer(y, y) # n*n matrix, each entry is the phenotype correlation between two individuals
  
  # genotype correlations
  geno_corr = X%*%t(X)/p # n*n matrix, each entry is the genetic correlation between two individuals
  
  # multiply genotype correlations with the constant c
  c = P*(1 - P)*dnorm(t)^2/(K^2*(1 - K)^2)
  geno_corr_c = c*geno_corr
  
  # regress phenotype correlations on genotype correlations
  pheno_corr_vec = pheno_corr[upper.tri(pheno_corr)]
  geno_corr_c_vec = geno_corr_c[upper.tri(geno_corr_c)]
  
  data = data.frame(y=pheno_corr_vec, X=geno_corr_c_vec)
  fit = lm(y~X, data=data)
  h2l = summary(fit)$coefficients[2,1] # liability scale heritability
  h2l_se = summary(fit)$coefficients[2,2] # standard error of h2l
  
  res = list(h2l=h2l, h2l_se=h2l_se)
  return(res)
  
}

# parameters
K = 0.1
P = 0.3
t = qnorm(1 - K)
h2l = 0.5
p = 10000

# sample size
n = 4000
n_case = n*P
n_control = n*(1 - P)

mafs = runif(p, 0.05, 0.5)

## generate data
beta = rnorm(p, 0, sqrt(h2l/p)) # SNP effect sizes

n1 = n0 = 0
G1 = X1 = matrix(NA, n_case, p)
G0 = X0 = matrix(NA, n_control, p)
l1 = l0 = c()

while(n1 < n_case | n0 < n_control){
  # generate a genotype
  g = rbinom(p, 2, mafs)
  # standarized the genotype
  x = (g - 2*mafs)/sqrt(2*mafs*(1 - mafs))
  # liability for this individual
  l = sum(x*beta) + rnorm(1, 0, sqrt(1-h2l))
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
y_std = (y - P)/sqrt(P*(1 - P))

## estimate heritability using pcgc regression
fit_pcgc = pcgc_reg(y=y_std, X=X, K=K, P=P)
fit_pcgc$h2l


### liability plot
# paramters
h2l = 0.5
GE = rmvnorm(1000000, c(0, 0), diag(c(h2l, 1-h2l)))
g = GE[,1]
e = GE[,2]
l = g + e
# parameters
K = 0.01
t = qnorm(1 - K)
P = 0.5
n = 2000 # total sample size in the study
n_case = n*P # number of cases in the study
n_control = n*(1 - P) # number of controls in the study
# case-control sampling
idx_case_popu = which(l > t)
idx_control_popu = which(l < t)
idx_case_study = sample(idx_case_popu, n_case)
idx_control_study = sample(idx_control_popu, n_control)
# liability distribution in the study
l_case_study = l[idx_case_study]
l_control_study = l[idx_control_study]
l_study = c(l_case_study, l_control_study)
hist(main='Liability density under ascertainment', 
     x=l_study, breaks=50, xlab='Liability')



### signal-noise ratio
n = 5000
m = 1000
h2l = 0.5
# effect sizes
beta = rnorm(m, 0, sqrt(h2l/m))
# design matrix
X = matrix(rnorm(n*m), n, m)
# error
e = rnorm(n, 0, sqrt(h2l))
# response
y = X%*%beta + e
# linear regression on each variant
beta_hat = c()
beta_hat_se = c()
for(k in 1:m){
  fit = lm(y~., data=data.frame(y, X[,k]))
  beta_hat[k] = summary(fit)$coefficients[2,1]
  beta_hat_se[k] = summary(fit)$coefficients[2,2]
}
z = beta_hat/beta_hat_se
var(z)

pval = 2*pnorm(abs(z), lower.tail=F)
sum(beta_hat^2)

# beta_hat2 = drop(t(X)%*%solve(X%*%t(X))%*%y)
# beta_hat2 = drop(solve(t(X)%*%X)%*%t(X)%*%y)
# sum(beta_hat2^2)

thres = 0.05/m
idx_sig = which(pval < thres)
length(idx_sig)
sig_beta = beta_hat[idx_sig]
cat(sum(sig_beta^2), sum(beta[idx_sig]^2))

### demo for presentation
# decrease SNR by increasing m or decreasing n (it is better to fix n, so that the threshold do not need to change)
# set a threshold -> calculate sum of beta_hat of these significant SNPs 
# each setting: plot z-score (SNR), number of significant SNPs, variance explained by these SNPs
m = 1000
thres = 0.05/m # Bonferroni correction
ns = c(500, 2000, 10000, 20000)
h2l = 0.5

nrep = 10

RES = data.frame(matrix(NA, length(ns)*nrep, 7))
names(RES) = c('m', 'n', 'rep', 'num_sig', 'sig_h2', 'sig_betasq', 'sum_all')
RES$m = m
RES$n = rep(ns, each=nrep)
RES$rep = rep(seq(nrep), length(ns))
for(n in ns){
  # design matrix
  X = matrix(rnorm(n*m), n, m)
  for(rep in 1:nrep){
    # effect sizes
    beta = rnorm(m, 0, sqrt(h2l/m))
    # error
    e = rnorm(n, 0, sqrt(h2l))
    # response
    y = X%*%beta + e
    # linear regression on each variant
    beta_hat = c()
    beta_hat_se = c()
    for(k in 1:m){
      fit = lm(y~., data=data.frame(y, X[,k]))
      beta_hat[k] = summary(fit)$coefficients[2,1]
      beta_hat_se[k] = summary(fit)$coefficients[2,2]
    }
    z = beta_hat/beta_hat_se
    pval = 2*pnorm(abs(z), lower.tail=F)
    idx_sig = which(pval < thres)
    sig_beta = beta_hat[idx_sig]
    RES[RES$n==n&RES$rep==rep,]$num_sig = length(idx_sig)
    RES[RES$n==n&RES$rep==rep,]$sig_h2 = sum(sig_beta^2)
    RES[RES$n==n&RES$rep==rep,]$sig_betasq = sum(beta[idx_sig]^2)
    RES[RES$n==n&RES$rep==rep,]$sum_all = sum(beta_hat^2)
    message(sprintf('n = %s, rep = %s', n, rep))
  }
}

P_all = ggplot(data=RES, aes(x=factor(n), y=sum_all, color=factor(n))) +
  geom_boxplot() +
  geom_hline(yintercept=h2l, lty='dashed') +
  scale_y_continuous(breaks = sort(c(seq(0, 3, length.out=4), h2l))) +
  xlab('n') + ylab('Estimation') +
  ggtitle('Sum of square of all estimated effect sizes: m = 1000') +
  theme(plot.title=element_text(hjust=0.5), legend.position='none')
P_all
ggsave('/Users/yxy/PCGC/sum_all_m1000.png', P_all, width=6, height=5)

P_est = ggplot(data=RES, aes(x=factor(n), y=sig_h2, color=factor(n))) +
  geom_boxplot() +
  geom_hline(yintercept=h2l, lty='dashed') +
  xlab('n') + ylab('Estimation') +
  ggtitle('Heritability explained by significant SNPs: m = 1000') +
  theme(plot.title=element_text(hjust=0.5), legend.position='none')
P_est

P_sig = ggplot(data=RES, aes(x=factor(n), y=num_sig, color=factor(n))) +
  geom_boxplot() +
  xlab('n') + ylab('Number') +
  ggtitle('Number of significant SNPs: m = 1000') +
  theme(plot.title=element_text(hjust=0.5), legend.position='none')
P_sig
ggsave('/Users/yxy/PCGC/num_sigsnp_m1000.png', P_sig)

P_beta = ggplot(data=RES, aes(x=factor(n), y=sig_betasq, color=factor(n))) +
  geom_boxplot() +
  xlab('n') + ylab('Signal') +
  ggtitle('Signal of significant SNPs: m = 1000') +
  theme(plot.title=element_text(hjust=0.5), legend.position='none')
P_beta

PLOT = data.frame(matrix(NA, length(ns)*nrep*2, 5))
names(PLOT) = c('m', 'n', 'rep', 'method', 'value')
PLOT$m = m
PLOT$n = rep(ns, each=nrep*2)
PLOT$rep = rep(rep(seq(nrep), each=2), length(ns))
PLOT$method = c('sig_betasq', 'sig_h2')
PLOT[PLOT$method=='sig_betasq', 'value'] = RES$sig_betasq
PLOT[PLOT$method=='sig_h2', 'value'] = RES$sig_h2

P_both = ggplot(data=PLOT, aes(x=factor(n), y=value, color=method)) +
  geom_boxplot() +
  geom_hline(yintercept=h2l, lty='dashed') +
  xlab('n') + ylab('Value') +
  ggtitle('Aggregate effect size of significant SNPs: m = 1000') +
  theme(plot.title=element_text(hjust=0.5))
ggsave('/Users/yxy/PCGC/h2_sigsnp_m1000.png', P_both)
P_both

### normal, binomial
m = 1000
mafs = runif(m, 0, 1)
gs = rbinom(m, 2, mafs)
xs = (gs - 2*mafs)/sqrt(2*mafs*(1 - mafs))
hist(xs, breaks=50)

