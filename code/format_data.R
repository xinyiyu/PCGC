### format simulated genotype to plink data format
library(genio)
out_data = '/import/home/xyubl/PCGC/sim_data/no_covar/'

## generate genotype matrix: 
# use generative model (for scenarios with a small number of SNPs, e.g 100 or 1,000 SNPs)
p = 1000 # number of SNPs
mafs = runif(p, 0.05, 0.5)
hlsq = 0.5 # true liability scale heritability
beta = rnorm(p, 0, sqrt(hlsq/p)) # SNP effect sizes
K = 0.1 # population prevalence
t = qnorm(1 - K) # liability threshold
P = 0.3 # in-sample prevalence

# full ascertainment (if the generated individual is a case, then she will be included in the study)
n = 2000
n_case = n*P
n_control = n*(1 - P)

c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2) # correction constant between hlsq and hosq

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
y_std = (y - P)/sqrt(P*(1 - P))

## generate data.frame as .bim: chr, id, posg, pos, ref, alt
bim_df = data.frame(matrix(NA, p, 6))
names(bim_df) = c('chr', 'id', 'posg', 'pos', 'ref', 'alt')
bim_df$chr = 1
bim_df$id = sprintf('rs%s', seq(p))
bim_df$posg = seq(p)
bim_df$pos = seq(p)
bim_df$ref = rep(c('A', 'C', 'G', 'T'), p/4)
bim_df$alt = rep(c('T', 'G', 'C', 'A'), p/4)

## generate data.frame as .fam: fam, id, pat, mat, sex, pheno
fam_df = data.frame(matrix(NA, n, 6))
names(fam_df) = c('fam', 'id', 'pat', 'mat', 'sex', 'pheno')
fam_df$fam = sprintf('F%s', seq(n))
fam_df$id = seq(n)
fam_df$pat = 0
fam_df$mat = 0
fam_df$sex = 0
fam_df$pheno = -9

## get plink data format using write_plink()
write_plink(file=paste0(out_data, 'test_plink_data'), X=t(G), bim=bim_df, fam=fam_df)

## generate .phen file
phen_df = data.frame(fam=fam_df$fam, id=fam_df$id, pheno=y_std)
write_phen(file=paste0(out_data, 'test_plink_data.phen'), tib=phen_df)

## using GCTA to calculate heritability
# make .grm file
./gcta64 \
--bfile /import/home/xyubl/PCGC/sim_data/no_covar/test_plink_data \
--make-grm \
--out /import/home/xyubl/PCGC/sim_data/no_covar/test_plink_data

# reml
./gcta64 \
--grm /import/home/xyubl/PCGC/sim_data/no_covar/test_plink_data \
--pheno /import/home/xyubl/PCGC/sim_data/no_covar/test_plink_data.phen \
--reml \
--out /import/home/xyubl/PCGC/sim_data/no_covar/test_plink_data
