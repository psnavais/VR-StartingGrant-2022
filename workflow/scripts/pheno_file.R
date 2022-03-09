library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]])

names(ids)= c('PREG_ID', 'IID', 'BATCH', 'Role')

ids= filter(ids, BATCH != 'TED')

ids= ids[!duplicated(ids[,c('PREG_ID', 'IID')]), ]

flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

out= readLines(snakemake@input[[4]])


ids= filter(ids, !(IID %in% out), IID %in% flag$IID)

ids= spread(ids, key= Role, value= IID)

mfr= inner_join(mfr, ids, by= c('PREG_ID_315'= 'PREG_ID'))


mfr= filter(mfr, !is.na(SPABORT_12_5), !is.na(SPABORT_23_5))

mfr$SPABORT_12_5= as.numeric(ifelse(mfr$SPABORT_12_5== '4 eller flere', '4', mfr$SPABORT_12_5))


mfr$SPABORT_23_5= as.numeric(ifelse(mfr$SPABORT_23_5== '4 eller flere', '4', mfr$SPABORT_23_5))

mfr$parity= with(mfr, ifelse(PARITET_5== '0 (førstegangsfødende)', 0, ifelse(PARITET_5== '1', 1, ifelse(PARITET_5== '2', 2, ifelse(PARITET_5== '3', 3, ifelse(is.na(PARITET_5), NA, 4)) ))))

mfr$tot_miscarriage= rowSums(mfr[,c('SPABORT_12_5', 'SPABORT_23_5')], na.rm=T)

mfr$tot_miscarriage_log= with(mfr, ifelse(tot_miscarriage> 1, 1, 0))


mfr$MOR_AGE= mfr$FAAR - mfr$MOR_FAAR



mfr= arrange(mfr, desc(parity))
mfr$cohort= mfr$BATCH


mfr$frac_miscarriage= with(mfr, tot_miscarriage / parity)
mfr$frac_miscarriage= with(mfr, ifelse(is.na(frac_miscarriage), 0, frac_miscarriage))

mfr$live_to_misc= with(mfr, ifelse(is.na(tot_miscarriage) | is.na(parity), NA, ifelse(tot_miscarriage>= (parity + 1), 1, 0)))
mfr$live_to_misc12= with(mfr, ifelse(is.na(SPABORT_12_5) | is.na(parity), NA, ifelse(SPABORT_12_5>= (parity + 1), 1, 0)))
mfr$live_to_misc23= with(mfr, ifelse(is.na(SPABORT_23_5) | is.na(parity), NA, ifelse(SPABORT_23_5>= (parity + 1), 1, 0)))

mfr$SPABORT_12_5= with(mfr, ifelse(SPABORT_12_5>0, 1, ifelse(is.na(SPABORT_12_5), NA, 0)))
mfr$SPABORT_23_5= with(mfr, ifelse(SPABORT_23_5>0, 1, ifelse(is.na(SPABORT_23_5), NA, 0)))

mfr$miscarriage= with(mfr, ifelse(SPABORT_12_5>0 | SPABORT_23_5>0, 1, ifelse(is.na(SPABORT_12_5) & is.na(SPABORT_23_5), NA, 0)))

mfr= mfr[order(mfr$miscarriage, decreasing= T),]

moms= mfr[!duplicated(mfr$Mother, incomparables= NA), ]
fets= mfr[!duplicated(mfr$Child, incomparables= NA), ]
dads= mfr[!duplicated(mfr$Father, incomparables= NA), ]

moms= select(moms, Mother, miscarriage, cohort) %>% filter(!is.na(Mother))

fets= select(fets, Child, miscarriage, cohort) %>% filter(!is.na(Child))
dads= select(dads, Father, miscarriage, cohort) %>% filter(!is.na(Father))

names(moms)[1]= 'IID'
names(fets)[1]= 'IID'
names(dads)[1]= 'IID'

fwrite(moms, snakemake@output[[1]], sep= '\t')
fwrite(fets, snakemake@output[[2]], sep= '\t')
fwrite(dads, snakemake@output[[3]], sep= '\t')
