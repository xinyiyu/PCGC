library(VCM)
library(ggplot2)
library(RhpcBLASctl)
source('/import/home/xyubl/PCGC/code/func.R')

set.seed(1)
# threads = 30
# blas_set_num_threads(threads)
nrep = 30
out_data = '/import/home/xyubl/PCGC/sim_data/'
out_result = '/import/home/xyubl/PCGC/sim_result/'
out_plot = '/import/home/xyubl/PCGC/sim_plot/'

##### no covariates #####
### PCGC vs REML: different hlsq, K and P (Fig 2A in paper)
hlsqs = c(0.1, 0.3, 0.5, 0.7, 0.9)
# Ks = c(0.001, 0.005, 0.01)
Ks = c(0.01, 0.05, 0.1)
Ps = c(0.1, 0.3, 0.5)
p = 10000
n = 4000

mafs = runif(p, 0.05, 0.5)

for(hlsq in 0.9){
  for(K in Ks){
    for(P in 0.1){
      
      n_case = n*P
      n_control = n*(1 - P)
      t = qnorm(1 - K)
      c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2) # correction constant between hlsq and hosq
      
      EST_HLSQ = data.frame(matrix(NA, nrep, 3))
      names(EST_HLSQ) = c('rep', 'pcgc', 'reml')
      EST_HLSQ$rep = seq(nrep)
      
      for(rep in 1:nrep){
        
        beta = rnorm(p, 0, sqrt(hlsq/p)) # SNP effect sizes
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
        
        message(sprintf('Finish generating data: %s-th rep, hlsq = %s, K = %s, P = %s.', rep, hlsq, K, P))
        
        ## pcgc regression
        res_pcgc = pcgc_reg(y_std, X, K, P, t)
        EST_HLSQ[EST_HLSQ$rep==rep,]$pcgc = res_pcgc$hlsq
        message(sprintf('Finish pcgc regression: %s-th rep, hlsq = %s, K = %s, P = %s.', rep, hlsq, K, P))
        
        ## reml
        res_reml = linRegPXEM(X=X, y=y, tol=1e-5, verbose=T)
        hosq_reml = res_reml$sb2/(res_reml$sb2 + res_reml$se2)
        EST_HLSQ[EST_HLSQ$rep==rep,]$reml = c*hosq_reml
        message(sprintf('Finish reml estimation: %s-th rep, hlsq = %s, K = %s, P = %s.', rep, hlsq, K, P))
        
      }
      
      # data not saved
      # save result
      save(EST_HLSQ, file=paste0(out_result, sprintf('no_covar/hlsq_pcgc_vs_reml_hlsq_%s_K_%s_P_%s.rda', hlsq, K, P)))
    
    }
  }
}


