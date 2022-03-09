library(data.table)
library(dplyr)
library(tidyr)
library(car)

pheno= fread(snakemake@input[[1]])

write( paste('snp', 'n', 'beta_h1', 'se_h1', 'pvalue_h1', 'beta_h2', 'se_h2', 'pvalue_h2', 'beta_h3', 'se_h3', 'pvalue_h3', 'beta_h4', 'se_h4', 'pvalue_h4', 'pvalue_maternal', 'pvalue_fetal', 'pval_poe', sep= '\t'), snakemake@output[[1]], append= T)

cols= grep('chr', names(pheno))
cols= unique(gsub('_h.', '', names(pheno)[cols], fixed= F))

results_list= lapply(cols, function(col) {

x= pheno %>% select(names(pheno)[grep(col, names(pheno))], PREG_ID)
names(x)= c('h1', 'h2', 'h3', 'h4', 'PREG_ID')
x= inner_join(pheno, x, by= 'PREG_ID')


m1= glm(miscarriage~ h1 + h2 + h3 + h4 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, x, family= binomial)

print(col)

n= length(resid(m1))
coefs= summary(m1)$coefficients[2:5,]
beta_h1= coefs[1,1]
se_h1= coefs[1,2]
pvalue_h1= coefs[1,4]
beta_h2= coefs[2,1]
se_h2= coefs[2,2]
pvalue_h2= coefs[2,4]
beta_h3= coefs[3,1]
se_h3= coefs[3,2]
pvalue_h3= coefs[3,4]
beta_h4= coefs[4,1]
se_h4= coefs[4,2]
pvalue_h4= coefs[4,4]
snp= col

pval_mat= tryCatch(linearHypothesis(m1, 'h1 - h3 + h2 = 0')[['Pr(>Chisq)']][2], warning= function(w){NA}, error= function(w) {NA})
pval_fet= tryCatch(linearHypothesis(m1, 'h1 - h2 + h3 = 0')[['Pr(>Chisq)']][2], warning= function(w){NA}, error= function(w) {NA})
pval_poe= tryCatch(linearHypothesis(m1, 'h1 - h2 = h3')[['Pr(>Chisq)']][2], warning= function(w){NA}, error= function(w) {NA})


results= paste(snp, n, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, beta_h4, se_h4, pvalue_h4, pval_mat, pval_fet, pval_poe, sep= '\t')
write(results, file= snakemake@output[[1]], append=TRUE)

}

)



