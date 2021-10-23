library(ggplot2)

### analyze result 
out_result = '/import/home/xyubl/PCGC/sim_result/no_covar/'
out_plot = '/import/home/xyubl/PCGC/sim_plot/no_covar/'

# load results (currently available), each file saves a nrep*3 dataframe: EST_HLSQ (rep, pcgc, reml)
hlsqs = c(0.1, 0.3, 0.5, 0.7, 0.9)
Ks = c(0.01, 0.05, 0.1)
Ps = c(0.1, 0.3, 0.5)
nrep = 30

EST = data.frame(matrix(NA, 5*3*3*nrep*2, 6))
names(EST) = c('hlsq', 'K', 'P', 'rep', 'method', 'estimate')
EST$hlsq = rep(hlsqs, each=3*3*nrep*2)
EST$K = rep(rep(Ks, each=3*nrep*2), 5)
EST$P = rep(rep(Ps, each=nrep*2), 5*3)
EST$rep = rep(rep(seq(nrep), each=2), 5*3*3)
EST$method = rep(c('PCGC', 'REML'), 5*3*3*nrep)
for(hlsq in hlsqs){
  for(K in Ks){
    for(P in Ps){
      
      # load result under current setting
      load(file=paste0(out_result, sprintf('hlsq_pcgc_vs_reml_hlsq_%s_K_%s_P_%s.rda', hlsq, K, P)))
      
      EST[EST$hlsq==hlsq&EST$K==K&EST$P==P&EST$method=='PCGC',]$estimate = EST_HLSQ$pcgc
      EST[EST$hlsq==hlsq&EST$K==K&EST$P==P&EST$method=='REML',]$estimate = EST_HLSQ$reml
      
    }
  }
}

# boxplot
g_est = ggplot(data=EST, aes(x=factor(hlsq), y=estimate, color=method)) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(rows=vars(P), cols=vars(K)) +
  xlab(expression(paste('True ', h[l]^2))) + ylab(expression(paste('Estimated ', h[l]^2))) + 
  ggtitle(expression(paste(h[l]^2, ' Estimation: PCGC vs REML'))) +
  theme(plot.title=element_text(hjust=0.5))
ggsave(filename=paste0(out_plot, 'hlsq_est_pcgc_vs_reml.png'), g_est)

