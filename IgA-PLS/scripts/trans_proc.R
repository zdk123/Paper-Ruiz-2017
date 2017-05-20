
### OUT of sample prediction using TransPAT
transY <- sample_data(phylo_list$tpat)
Ytrans <- cbind(as.numeric(transY$IgA),NA)
Ytrans.split <- compPLS:::.Split.variation.one.level(Ytrans,
       as.factor(transY$Treatment), as.numeric(as.factor(transY$sample_name)))
Xtrans <- t(clr(phylo_list$tpat@otu_table@.Data, 2))
Xtrans.split <- compPLS:::.Split.variation.one.level(Xtrans,
       as.factor(transY$Treatment), as.numeric(as.factor(transY$sample_name)))

Ytrans <- scale(Ytrans.split$Xw[,1,drop=FALSE],scale=FALSE, center=FALSE)
Xtrans <- scale(Xtrans.split$Xw,               scale=FALSE, center=TRUE)
