---
title: "IGA-PLS"
author: "Zachary Kurtz"
output:
  pdf_document:
    highlight: tango
    toc: true
---



### Introduction
This Notebook is meant to contain the 16S_OTU data found within the respective manuscript. It contains only the code used to generate the figures found within the main text. This entire dataset is publicly available in QIITA under the ID [10527](https://qiita.ucsd.edu/study/description/10527). More details about the procedures used to generate the data can be found within the **Methods** section of the manuscript. The table found within the data folder has been processed to include consensus set of taxa.

#### 1a. Install and load the necessary libraries

```r
library(stringr)
library(phyloseq)
# devtools::install_github('zdk123/compPLS')
library(compPLS)
library(spls)
```
#### 1b. Import the OTU and sample data

```r
pattern  <- "(.*)_cons\\.biom"
pattern2 <- paste('data', pattern, sep="/")
files    <- list.files('data', pattern, full.names=TRUE)
phylo_list <- lapply(files, import_biom)
names(phylo_list) <- str_match(files, pattern2)[,2]
```

```
## Error in str_match(files, pattern2)[, 2]: subscript out of bounds
```

#### prep data for PLS model ####

```r
patY <- sample_data(phylo_list$expat)
```

```
## Error in access(object, "sam_data", errorIfNULL): sam_data slot is empty.
```

```r
Y <- cbind(as.numeric(patY$IgA),NA)
Y.split <- compPLS:::.Split.variation.two.level(Y,
       as.factor(patY$Treatment), as.factor(patY$Diet),
       as.factor(patY$sample_name))
```

```
## Warning: the vector sample was converted into a numeric vector
```

```r
Y <- scale(Y.split$Xw[,1,drop=FALSE])
X <- t(clr(phylo_list$expat@otu_table@.Data, 2))
```

```
## Error in clr(phylo_list$expat@otu_table@.Data, 2): trying to get slot "otu_table" from an object of a basic class ("NULL") with no slots
```

```r
X.split <- compPLS:::.Split.variation.two.level(X,
       as.factor(patY$Treatment), as.factor(patY$Diet),
       as.factor(patY$sample_name))
```

```
## Warning: the vector sample was converted into a numeric vector
```

```
## Error in is.data.frame(x): object 'X' not found
```

```r
X <- scale(X.split$Xw, scale=FALSE)
```

```
## Error in scale(X.split$Xw, scale = FALSE): object 'X.split' not found
```

#### run model, fit to transfer experiments ####

```r
K <- 2
out.stars <- compPLS:::spls.stars(X, Y, rep.num=100, K=K,
             eta=seq((.799), (.999), length.out=10), ncores=1)
```

```
## eta = 0.799 
## eta = 0.821222222222222 
## eta = 0.843444444444444 
## eta = 0.865666666666667 
## eta = 0.887888888888889 
## eta = 0.910111111111111 
## eta = 0.932333333333333 
## eta = 0.954555555555556 
## eta = 0.976777777777778 
## eta = 0.999
```

```r
p <- ncol(X)
threshs <- seq(0,1,.05)
n <- length(threshs)
out.spls_list <- vector('list', n)
Ypred_list    <- vector('list', n)
rsquared_in   <- vector('numeric', n)
rsquared_out  <- vector('numeric', n)
ntax <- vector('numeric', n)
source('scripts/trans_proc.R')
```

```
## Warning in file(filename, "r", encoding = encoding): cannot open file
## 'scripts/trans_proc.R': No such file or directory
```

```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
```

```r
for (i in 1:n) {
    th      <- threshs[i]
    ind     <- which(out.stars$merge[,1,1] >= th)
    ntax[i] <- length(ind)
    out.spls_list[[i]] <- spls::spls(X[,ind], Y, eta=0,
                          K=ifelse(K<ntax[i], K, ntax[i]))
    out.spls         <- out.spls_list[[i]]
    Ypred_list[[i]]  <- X[,ind] %*% out.spls_list[[i]]$betahat
    rsquared_in[i] <- cor(Y, Ypred_list[[i]])^2
    Ytranspred <- Xtrans[,ind,drop=FALSE] %*%
                  out.spls$betahat[,1,drop=FALSE]
    rsquared_out[i] <- cor(Ytrans[,1], Ytranspred[,1])^2
}
```
