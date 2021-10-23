library(data.table)

readGRM <- function(rootname)
{
  bin.file.name <- paste0(rootname, ".grm.bin", sep="")
  n.file.name <- paste0(rootname, ".grm.N.bin", sep="")
  id.file.name <- paste0(rootname, ".grm.id", sep="")
  
  cat("Reading IDs\n")
  id <- read.table(id.file.name, colClasses="character")
  n <- dim(id)[1]
  cat("Reading GRM\n")
  bin.file <- file(bin.file.name, "rb")
  grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(bin.file)
  cat("Reading N\n")
  n.file <- file(n.file.name, "rb")
  N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(n.file)
  
  cat("Creating data frame\n")
  l <- list()
  for(i in 1:n)
  {
    l[[i]] <- 1:i
  }
  col1 <- rep(1:n, 1:n)
  col2 <- unlist(l)
  grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	
  
  ret <- list()
  ret$grm <- grm
  ret$id <- id
  return(ret)
}

transformPheno <- function(phen,fe,K){
  # fit logistic regression 
  model <- glm(phen ~ fe , family="binomial")
  asc.probs <- model$fitted
  if(length(asc.probs) != length(phen)){
    # we have NAs 
    new.asc.probs <- array(NA,length(phen))
    new.asc.probs[which(!is.na(phen))] <- asc.probs
    asc.probs <- new.asc.probs 
  }
  P <- mean(phen,na.rm=T)
  asc.const <- (1-P)/P * K/(1-K)
  non.asc.probs <- asc.const * asc.probs / (1+asc.const * asc.probs - asc.probs)
  # convert probs to thresholds 
  non.asc.thresholds <- qnorm(1-non.asc.probs)
  # correct ascertainment
  phis <- dnorm(non.asc.thresholds)
  transpheno = (phen-asc.probs)/(sqrt(asc.probs*(1-asc.probs)));
  
  const1 <- K/(1-K) * (1-P)/P
  const2 <- (P-K)/(P*(1-K)) 
  
  x <- (1 - (asc.probs) * const2 )
  x <- x / sqrt(asc.probs*(1- asc.probs))
  x <- x * phis 
  x = x / (non.asc.probs*(1-const1)+const1)	
  transpheno = transpheno;
  
  # compute the population-wise variance of the thresholds using the law of total-variance
  non.asc.thresholds[non.asc.thresholds==Inf] = NA;
  var.cases <- var(non.asc.thresholds[phen==1],na.rm=T)
  var.controls <- var(non.asc.thresholds[phen==0],na.rm=T)
  mean.cases <- mean(non.asc.thresholds[phen==1],na.rm=T)
  mean.controls <- mean(non.asc.thresholds[phen==0],na.rm=T)
  total.var <- K*var.cases + (1-K)*var.controls + K*(1-K)*(mean.cases - mean.controls)^2
  
  ret <- list(phen=transpheno,multiplier=x,totalvar=total.var)
  return(ret);
  
}

pcgc_with_grm_covar = function(y, grm, c, Z, lambda, t_var){
  
  # y: standardized phenotype
  # grm: vectorized GRM
  # Z: eigenvectors
  # lambda: eigenvalues
  
  # phenotype correlations
  # standardize y
  pheno_corr = outer(y, y) # n*n matrix, each entry is the phenotype correlation between two individuals
  
  # regress phenotype correlations on genotype correlations
  pheno_corr_vec = pheno_corr[upper.tri(pheno_corr)]
  
  c_mat = outer(c, c)
  c_vec = c_mat[upper.tri(c_mat)]
  
  # remove PCs from grm
  eig_mat = Z%*%diag(lambda)%*%t(Z)
  eig_vec = eig_mat[upper.tri(eig_mat)]
  grm_new = grm - eig_vec
  geno_corr_c_vec = c_vec*grm_new
  
  data = data.frame(y=pheno_corr_vec, X=geno_corr_c_vec)
  fit = lm(y~X, data=data)
  sig2g = summary(fit)$coefficients[2,1] # sigma_g^2
  sig2g_se = summary(fit)$coefficients[2,2] # standard error of sigma_g^2
  
  # liability scale heritability
  h2l = sig2g/(1 + t_var)
  
  res = list(sig2g=sig2g, sig2g_se=sig2g_se, h2l=h2l)
  return(res)
  
}

# load GRM
GRM = readGRM('./PCGC/data/CADqc_withPC')
# get off-diagonal elements of GRM
idx_diag = which(GRM$grm$id1==GRM$grm$id2)
n = length(idx_diag)
idx_off = seq(n*(n+1)/2)[-idx_diag]
summary(GRM$grm$grm[idx_diag])
summary(GRM$grm$grm[idx_off])

# estimate variance-component by reml
system(paste0('/import/home/xyubl/gcta_1.93.2beta/gcta64', ' --grm ',
              wtccc_folder, trait, '/', name, 
              ' --pheno ', wtccc_folder, trait, '/', name, '.phen',
              ' --qcovar ', wtccc_folder, trait, '/', trait, 'pca20.eigenvec',
              ' --reml --out ', out_real, name))

# load .phen file: 2-case, 1-control
y = fread('./PCGC/data/CADqc_withPC.phen')$V3
y = ifelse(y==2, 1, 0)
P = sum(y)/length(y)
K = 0.05
t = qnorm(1 - K)

## considering PCs
# load eigenvectors and eigenvalues
PCs = fread('./PCGC/data/CADpca20.eigenvec')
lambda = fread('./PCGC/data/CADpca20.eigenval')$V1 # eigenval for each individual
PCmat = as.matrix(PCs[,-c(1,2)])
transphen = transformPheno(phen=y, fe=PCmat, K=K)

pcgc_est_pc = pcgc_with_grm_covar(y=transphen$phen, grm=GRM$grm$grm[idx_off], c=transphen$multiplier, Z=PCmat, lambda=lambda[1:20], t_var=transphen$totalvar)
pcgc_est_pc


# BD
K = 0.005 
P = 1812/4707
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.677580
h2l_reml = c*h2o_reml
h2l_reml

# CD
K = 0.001
P = 1675/4570
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.587139
h2l_reml = c*h2o_reml
h2l_reml

# T1D
K = 0.01
P = 1932/4825
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.864441
h2l_reml = c*h2o_reml
h2l_reml

# CAD
K = 0.05
P = 1888/4767
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.506738
h2l_reml = c*h2o_reml
h2l_reml

# T2D
K = 0.06
P = 1870/4770
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.471421
h2l_reml = c*h2o_reml
h2l_reml

# RA
K = 0.005
P = 1812/4709
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.574399
h2l_reml = c*h2o_reml
h2l_reml

# HT
K = 0.06
P = 1892/4788
t = qnorm(1 - K)
c = K^2*(1 - K)^2/(P*(1 - P)*dnorm(t)^2)

h2o_reml = 0.481618
h2l_reml = c*h2o_reml
h2l_reml



# keep samples with PCs
PC = fread('/home/share/WTCCC/HT/HTpca20.eigenvec')
keep = PC[,c(1,2)]
fwrite(keep, '/home/xyubl/PCGC/real_data/HTqc_withPC.txt', sep=' ', col.names=F)
pheno = fread('/home/share/WTCCC/HT/HTqc.phen')
pheno_withPC = merge(keep, pheno, by=c('V1', 'V2'), sort=F)
fwrite(pheno_withPC, file='/home/xyubl/PCGC/real_data/HTqc_withPC.phen', sep=' ', col.names=F)
plink \
--bfile /import/home/share/WTCCC/HT/HTqc \
--keep /import/home/xyubl/PCGC/real_data/HTqc_withPC.txt \
--make-bed \
--out /import/home/xyubl/PCGC/real_data/HTqc_withPC

# make .grm file
./gcta64 \
--bfile /import/home/xyubl/PCGC/real_data/HTqc_withPC \
--make-grm \
--out /import/home/xyubl/PCGC/real_data/HTqc_withPC



