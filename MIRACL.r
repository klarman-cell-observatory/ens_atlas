library(dnormpar2)
library(Matrix)

# -------------------------------------------------------------------------------
# Code for:
# Drokhlyansky, Smillie, van Wittenberghe, et al. Cell, 2020
#
# Description:
# This file contains code for detecting heterotypic multiplets in MIRACL-Seq data
# Code for other analyses found throughout the paper is available here:
# https://github.com/cssmillie/ulcerative_colitis
# For other questions, please email Chris Smillie (cssmillie [AT] gmail.com)
#
# Requirements:
# Before using this code, you must install the "dnormpar2" R package (included):
# > R CMD INSTALL dnormpar2.tar.gz
# This contains a fast method for calculating the doublet statistics
#
# Detailed example:
# This code can be run using the "find_doublets" function, which iteratively
# detects doublets for a series of given odds ratios. By starting with high odds
# ratios, we can initially remove only the doublets we are most confident in. This
# let's us refine our estimates for the other parameters, allowing us to gradually
# ease the odds ratio cutoff, only removing the highest confidence doublets at
# each step. To do this, just run the following code:
# > res = find_doublets(data, ident, sig.use, pd=.25, min_odds=c(8,7,6,5,4,3,2,1))
# The result is a list, which for each iteration, returns important information.
# To get the final doublet estimates, you therefore want the last list element:
# > res.final = res[[length(res)]]
# To get the final matrix of odds ratios for each cell type:
# > res.final$odds
# To get the final list of predicted doublets:
# > res.final$doublets
#
# Quick example:
# > res = find_doublets(data, ident, sig.use, pd=.25, min_odds=c(8,7,6,5,4,3,2,1))
# > doublets = res[[length(res)]]$doublets
#
# Please see "score_doublets" and "find_doublets" for additional information
# -------------------------------------------------------------------------------


score_doublets = function(data, ident, sig, scores=NULL, doublets=NULL, pd=0.25, min_exp=1, min_odds=1){

    # Given normalized expression data and a list of gene signatures, calculate doublet score for each cell
    # -----------------------------------------------------------------------------------------------------
    # Input arguments:
    # - data = log-normalized expression matrix (rows = genes, cols = cell barcodes)
    # - ident = vector of cell type names
    #   example: c('Epithelial', 'Goblet', 'Epithelial', 'Neuron', ...)
    # - sig = named list of gene signatures
    #   example: list(Epithelial=c('Epcam', 'Slc26a3', ...), Goblet=c('Muc2', 'Tff3', ...), Neuron=c('Elavl3', 'Elavl4', 'Chat', ...), ...)
    # - scores = current doublet scores (internal)
    # - doublets = current doublets (internal)
    # - pd = prior probability of a doublet
    # - min_exp = minimum log-normalized expression for a gene signature (ignore signatures with low expression)
    # - min_odds = minimum (doublet/singlet) odds ratio to be considered a doublet
    
    # check input data
    # ----------------
    
    # initialize putative doublets
    if(is.null(doublets)){
        doublets = rep(FALSE, ncol(data))
    }
    doublets[is.na(doublets)] = FALSE
    
    # check data dimensions
    if(ncol(data) != length(ident)){
        quit('error: check dimensions. number of cells in "data" and "ident" do not agree')
    }
    ident = as.factor(ident)

    # check gene signature names
    sig = sig[names(sig) %in% levels(ident)]
    if(length(sig) == 0){
        quit('error: names of cell types and gene signatures must match')
    }
    
    # summarize parameters
    # --------------------
    print(paste('Using gene signatures for:', paste(names(sig), collapse=', ')))
    print(paste('Filtering', sum(doublets), 'doublets'))
    print(paste('Ignoring signatures with log2(TP10K+1) <', min_exp))
    print(paste('Odds ratio cutoff =', min_odds))
    
    # score gene signatures
    # ---------------------
    if(is.null(scores)){
        print('Calculating signature scores')
	scores = sapply(sig, function(a){
	    b = intersect(a, rownames(data))
	    colMeans(data[b,], na.rm=T)
	})
	colnames(scores) = make.names(colnames(scores))
    }
    i = !doublets
    scores.u = data.frame(aggregate(scores[i,], list(make.names(ident[i])), mean), row.names=1)
    scores.s = data.frame(aggregate(scores[i,], list(make.names(ident[i])), sd), row.names=1)

    # calculate doublet statistics
    # ----------------------------
    
    # parameters for null hypothesis (h0)
    h0.u = scores.u[make.names(as.character(ident)),]
    h0.s = scores.s[make.names(as.character(ident)),]
    
    # parameters for alternative hypothesis (h1)
    h1.u = t(matrix(diag(as.matrix(scores.u[colnames(scores.u),])), nrow=ncol(scores.u), ncol=nrow(scores)))
    h1.s = t(matrix(diag(as.matrix(scores.s[colnames(scores.s),])), nrow=ncol(scores.s), ncol=nrow(scores)))
    h1.u = .5*h0.u + .5*h1.u
    h1.s = .5*sqrt(h0.s**2 + h1.s**2)
    
    # probabilities for null hypothesis (h0)
    print('Calculating H0 probabilities')
    p0 = sapply(1:ncol(scores), function(j) {
        print(paste('Signature', colnames(scores)[j]))
	dnormpar2(scores[,j], h0.u[,j], h0.s[,j])
    })
    	
    # probabilities for alternative hypothesis (h1)
    print('Calculating H1 probabilities')
    p1 = sapply(1:ncol(scores), function(j) {
	print(paste('Signature', colnames(scores)[j]))
	dnormpar2(scores[,j], h1.u[,j], h1.s[,j])
    })
    
    # fix names
    rownames(h0.u) = rownames(h0.s) = rownames(h1.u) = rownames(h1.s) = rownames(p0) = rownames(p1) = rownames(scores)
    colnames(h0.u) = colnames(h0.s) = colnames(h1.u) = colnames(h1.s) = colnames(p0) = colnames(p1) = colnames(scores)
    
    # classify putative doublets
    # --------------------------
    
    # calculate odds ratio
    print('Detecting doublets')
    odds = (p1*pd)/(p0*(1-pd))
    odds[scores < min_exp] = 0
    
    # infer doublets
    doublets = apply(odds, 1, max) > min_odds
    
    # return doublet information
    # --------------------------
    # scores = gene signature scores for each cell
    # scores.u = mean gene signature score for each cell type
    # scores.s = standard deviation of gene signature scores for each cell type
    # h0.u, h0.s = mean and standard deviation for null hypothesis (singlet)
    # h1.u, h1.s = mean and standard deviation for alternative hypothesis (doublet)
    # p0, p1 = probability for null and alternative hypotheses
    # odds = odds ratios
    # doublets = list of putative doublets
    list(scores=scores, scores.u=scores.u, scores.s=scores.s, h0.u=h0.u, h0.s=h0.s, h1.u=h1.u, h1.s=h1.s, p0=p0, p1=p1, odds=odds, doublets=doublets)
	
}

find_doublets = function(data, ident, sig, scores=NULL, doublets=NULL, pd=0.25, min_exp=1, min_odds=c(4,2,1,1)){
    
    # Given normalized expression data and a list of gene signatures, iteratively infer putative doublets
    # ---------------------------------------------------------------------------------------------------
    # Input arguments:
    # - data = log-normalized expression matrix (rows = genes, cols = cell barcodes)
    # - ident = vector of cell type names
    #   example: c('Epithelial', 'Goblet', 'Epithelial', 'Neuron', ...)
    # - sig = named list of gene signatures
    #   example: list(Epithelial=c('Epcam', 'Slc26a3', ...), Goblet=c('Muc2', 'Tff3', ...), Neuron=c('Elavl3', 'Elavl4', 'Chat', ...), ...)
    # - scores = current doublet scores (internal)
    # - doublets = current doublets (internal)
    # - pd = prior probability of a doublet
    # - min_exp = minimum log-normalized expression for a gene signature (ignore signatures with low expression)
    # - min_odds = minimum (doublet/singlet) odds ratio to be considered a doublet
    #   note: if min_odds is a vector, then it iteratively calls doublets using different odds ratios
    
    # score gene signatures in each cell
    # ----------------------------------
    ident = as.factor(ident)
    sig = sig[names(sig) %in% levels(ident)]
    sig = sapply(sig, function(a) intersect(a, rownames(data)), simplify=F)
        
    if(is.null(scores)){
        print('Calculating signature scores')
        scores = sapply(sig, function(a){
            b = intersect(a, rownames(data))
            colMeans(data[b,,drop=F], na.rm=T)
        })
        colnames(scores) = make.names(colnames(scores))
    }
    if(is.null(doublets)){
        doublets = rep(FALSE, ncol(data))
    }
    scores = scores[,make.names(names(sig))]

    # iteratively calculate doublet statistics
    # ----------------------------------------
    res = list()
    for(i in 1:length(min_odds)){
        res[[i]] = score_doublets(data=data, ident=ident, sig=sig, scores=scores, doublets=doublets, pd=pd, min_exp=min_exp, min_odds=min_odds[[i]])
    }
    return(res)
}
