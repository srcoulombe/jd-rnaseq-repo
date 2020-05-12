library(data.table)
library(circlize)
library(mixOmics)
library(HTSFilter)
library(tidyr)
library(plyr)
library(genefilter)
library(edgeR)
library(limma)
library(Glimma)
library(DESeq2)

edgeR_DGE_analysis <- function( DGE_obj, design, alpha, coef=NULL, contrasts=NULL,
                                useLRT=FALSE, useQLF=FALSE, useEXACT=TRUE,
                                plotRejCurve=TRUE, filtering.methods=NA, quantiles=NA, 
                                pAdjustMethod="BH", verbose=TRUE) {
  DGE_obj <- calcNormFactors(DGE_obj)
  DGE_obj <- estimateGLMCommonDisp(DGE_obj, design, verbose=verbose)
  DGE_obj <- estimateGLMTrendedDisp(DGE_obj, design, verbose=verbose)
  DGE_obj <- estimateGLMTagwiseDisp(DGE_obj, design)
  contrastwise.edger.exactTest.IF.results <- list()
  for(i in 1:dim(contrasts)[2]) {
    edger.result <- edgeR_DGE(
      DGE_obj, design, contrast=contrasts[,i], 
      useLRT=FALSE, useQLF=FALSE, useEXACT=useEXACT
    )
    contrast.name <- colnames(contrasts)[i]
    print(contrast.name)
    edger.results.copy <- edger.result$edgeR.exactTest
    edger.results.copy$unadjPvalues <- edger.result$edgeR.exactTest$table$PValue
    edger.results.copy$counts <- DGE_obj$counts
    #title <- paste0(
    #  "edgeR(", 
    #  contrast.name, 
    #  ") rejection curve (alpha=",
    #  alpha,
    #  ")\n#rej=",
    #  sum(contrast.result['padj'][[1]] < fdr, na.rm=TRUE),
    #  ", threshold=",
    #  round(metadata(contrast.result)$filterThreshold,2)
    #)
    edger.exactTest.IF.results <- independent_filtering(
      edger.results.copy,
      filtering.methods,
      theta=quantiles,
      fdr=alpha,
      showplots=plotRejCurve,
      fromtool='edger'
    )
    print(edger.exactTest.IF.results)
    contrastwise.edger.exactTest.IF.results[[contrast.name]] <- edger.exactTest.IF.results
  }
  return(contrastwise.edger.exactTest.IF.results)
}

edgeR_DGE <- function(DGE_obj, design, coef=NULL, contrast=NULL,
                      useLRT=FALSE, useQLF=FALSE, useEXACT=TRUE) {
  # Wrapper function to run any of edgeR's three tests: exactTest, glmLRt, and glmQLFTest.
  #
  # PARAMETERS:
  #
  #   DGE_obj:  DGEList object on which to perform the statistical tests.
  #             see DGEList in the edgeR package
  #
  #   design:   design matrix to pass onto edgeR's GLM models.
  #             see ?model.matrix
  #
  #   coef, contrasts:  optional integer/character index vector (coef) and 
  #                     numeric vector/matrix (contrast) to pass onto edgeR's GLM tests.
  #                     default to NULL as per edgeR's defaults.
  #   useLRT, useQLF, useExact:   booleans indicating which tests to run.
  #
  # RETURNS:
  #
  #   list("edgeR.exactTest" = edgeR.y.exact, "edgeR.glm.lrt" = edgeR.y.glm.lrt, "edgeR.glm.qlf" = edgeR.y.glm.qlf)
  #
  #     where each item is either the object resulting from running that test, or NULL if that test was
  #     specified not to be run.
  #
  # TODO:
  #   passing more kwargs to the edgeR functions?
  
  if(!any(useEXACT, useLRT, useQLF)) {
    stop("at least one of `useEXACT`, `useLRT`, or `useQLF` must be TRUE")
  }
  
  if(useLRT){
    edgeR.y.glm.lrt <- glmFit(DGE_obj, design)
    if(is.null(coef)) coef=edgeR.y.glm.lrt$design
    edgeR.y.glm.lrt <- glmLRT(edgeR.y.glm.lrt, coef=coef, contrast=contrast)
  } else {
    edgeR.y.glm.lrt <- NULL
  }
  if(useEXACT) {
    edgeR.y.exact <- exactTest(DGE_obj, pair=names(contrast)[contrast!=0])  
  } else {
    edgeR.y.exact <- NULL
  }
  if(useQLF) {
    edgeR.y.glm.qlf <- glmQLFit(DGE_obj, design)
    if(is.null(coef)) coef=edgeR.y.glm.qlf$design
    edgeR.y.glm.qlf <- glmQLFTest(edgeR.y.glm.qlf, coef=coef, contrast=contrast)
  } else {
    edgeR.y.glm.qlf <- NULL
  }
  
  list_to_return = list( "edgeR.exactTest" = edgeR.y.exact, 
                         "edgeR.glm.lrt" = edgeR.y.glm.lrt, 
                         "edgeR.glm.qlf" = edgeR.y.glm.qlf)
  
  return(list_to_return)
}

DESeq2_DGE_analysis <- function(DESeq2_dataset, alpha, contrasts, 
                                plotRejCurve=TRUE, filtering.methods=NA, quantiles=NA, 
                                pAdjustMethod="BH") {
  #
  #
  # PARAMETERS:
  #
  #   DESeq2_dataset:  DESeqDataSet object on which to perform the statistical tests.
  #                   see `DESeqDataSet`` in the DESeq2 package.
  #
  #   alpha:  statistical significance threshold cutoff used to optimize the independent filtering.
  #           see `results`` from DESeq2 package.
  #
  #   pAdjustMethod: string indicating method to use in the p-value adjustment method.
  #                   see ?p.adjust for valid values.
  #                   defaults to "BH".
  #
  # RETURNS:
  #   
  #   result of DESeq2 analysis
  #
  # TODO:
  #   passing more kwargs to the DESeq2 functions
  
  DESeq2_dataset <- estimateSizeFactors(DESeq2_dataset)
  DESeq2.out <- DESeq(DESeq2_dataset)

  if(is.na(filtering.methods)) {
    filtering.methods <- data.frame(
      'mean'= rowMeans(counts(DESeq2_dataset, normalized=TRUE))
    )
  }

  contrastwise.standard.DESeq2.results <- list()
  contrastwise.IF.DESeq2.results <- list()
  for (contrast.name in colnames(contrasts)){
    contrast.conditions <- strsplit(contrast.name, "-")[[1]]
    contrast.result <- results(
      DESeq2.out, 
      contrast=c("condition", contrast.conditions[1], contrast.conditions[2]),
      alpha=alpha,
      pAdjustMethod=pAdjustMethod
    )
    contrastwise.standard.DESeq2.results[[contrast.name]] <- contrast.result
    
    if(plotRejCurve) {
      plottitle <- paste0(
        "DESeq2(", 
        contrast.name, 
        ") rejection curve (alpha=",
        alpha,
        ")\n#rej=",
        sum(contrast.result['padj'][[1]] < fdr, na.rm=TRUE),
        ", threshold=",
        round(metadata(contrast.result)$filterThreshold,2)
      )
     
      plot(
        metadata(contrast.result)$filterNumRej, 
        type="b", 
        ylab="number of rejections",
        xlab="quantiles of filter", 
        main=plottitle
      )
      lines(
        metadata(contrast.result)$lo.fit, 
        col="red"
      )
      abline(v=metadata(contrast.result)$filterTheta)
    }
  

    if(!any(is.na(filtering.methods), is.na(quantiles))){
      title <- paste0(
        "DESeq2(", 
        contrast.name, 
        ") rejection curve (alpha=",
        alpha,
        ")\n#rej=",
        sum(contrast.result['padj'][[1]] < fdr, na.rm=TRUE),
        ", threshold=",
        round(metadata(contrast.result)$filterThreshold,2)
      )
      contrast.result.copy <- contrast.result
      contrast.result.copy$unadjPvalues <- contrast.result$pvalue
      IF.DESeq2.results <- independent_filtering(
        contrast.result.copy,
        filtering.methods,
        theta=quantiles,
        title=title,
        fdr=fdr,
        showplots=plotRejCurve,
        fromtool="deseq"
      )
      contrastwise.IF.DESeq2.results[[contrast.name]] <- IF.DESeq2.results
    }
  }
  return(list(contrastwise.standard.DESeq2.results, contrastwise.IF.DESeq2.results))
}

design.pairs <- function(levels) {
  # credits to G. Smyth
  # https://support.bioconductor.org/p/9228/#9254
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    return(design)
}

independent_filtering <- function(dge_obj.with.pvalues, 
                                  filtering.methods.dataframe,
                                  theta=NA, fdr=0.05, 
                                  showplots=TRUE, p.adjust_method="BH", 
                                  title="Filtering Methods' Rejection Curves",
                                  fromtool=NA) {
  #
  # Wrapper function to recreate DESeq2's independent filtering code.
  #
  # PARAMETERS:
  #
  #   dge_obj.with.pvalues: DGEList object along with a unadjPvalues field.
  #
  #   filtering.methods.dataframe:  dataframe containing the methods to assess in the independent filtering.
  #                                 e.g.:
  #                                    filtering.methods <- data.frame(
  #                                      'mean' = rowMeans(dge_obj$counts),
  #                                      'min' = rowMin(dge_obj$counts),
  #                                      'max' = rowMax(dge_obj$counts),
  #                                      'median' = rowMedians(dge_obj$counts),
  #                                      'secondlargest' = apply(
  #                                        dge_obj$counts,
  #                                        1,
  #                                        function(row) sort(row,partial=length(row)-1)[length(row)-1]
  #                                      )
  #                                    )
  #
  #   theta: sequence ranging from [0.0, 1.0) indicating which quantiles to consider dropping.
  #                   defaults to [0.4, 0.95] by increments of 0.01
  #
  #   fdr:  value for the fdr threshold.
  #         defaults to 0.05. 
  #   
  #   showplots:  boolean indicating whether to display plots or not.
  #               defaults to TRUE.
  #
  #   p.adjust_method:  multiple-hypothesis p-value adjustment method.
  #                     gets passed onto filtered_R, see `genefilter` package.
  #                     defaults to "BH".
  #
  # RETURNS:
  #   
  #   list of "rejections" and "best.numRej.per.Method"
  #
  # REFERENCE:
  #   https://bioconductor.org/packages/release/bioc/vignettes/genefilter/inst/doc/independent_filtering.pdf
  
  if(is.na(theta)) {
    theta <- seq(from=0.4, to=0.95, by=0.01)
  }

  stopifnot(!is.na(fromtool) | !((tolower(fromtool) != "deseq") & (tolower(fromtool) != "edger")) )

  if(is.na(filtering.methods.dataframe)) {
    
    if(tolower(fromtool) == "deseq") {
      filtering.methods.dataframe <- data.frame(
        'mean'= rowMeans(counts(dge_obj.with.pvalues, normalized=TRUE))
      )
    }

    else {
      filtering.methods.dataframe <- data.frame(
        'mean'= rowMeans(dge_obj.with.pvalues$counts)
      )
    }    
  }

  rejections <- sapply(
    filtering.methods.dataframe,
    function(f) filtered_R(
      alpha=fdr, filter=f, 
      test=dge_obj.with.pvalues$unadjPvalues, 
      theta=theta, method=p.adjust_method
      )
    )

  best.numRej.per.method <- apply(rejections, 2, max, na.rm = TRUE)
  relevant.quantile.per.method <- numeric(length(filtering.methods.dataframe))
  relevant.threshold.per.method <- numeric(length(filtering.methods.dataframe))
  for(i in 1:length(filtering.methods.dataframe)) {
    # i-th element of `u1` squared into `i`-th position of `usq`
    q = theta[ which(rejections[,i] == max(rejections[,i])) ][1]
    # the extra [1] is to get the lowest quantile that yielded the max number of rejections
    # to filter out fewer genes
    
    relevant.quantile.per.method[i] = q
    relevant.column <- filtering.methods.dataframe[,colnames(filtering.methods.dataframe)[i]]
    relevant.threshold.per.method[i] = quantile(relevant.column,q)
    
  }
  if(showplots==TRUE) {
    colors=brewer.pal(ncol(filtering.methods), "Set1")
    matplot(theta, rejections, type='l', lty=1, col=colors, 
            lwd=2, xlab=expression(theta), ylab="Rejections", 
            main=title)
    legend("bottomleft", legend=colnames(filtering.methods), fill=colors)
  }
  #print(rejections)
  #print(relevant.quantile.per.method)
  #print(relevant.threshold.per.method)
  return.list <- list(
    "rejections" = rejections, 
    "best.numRej.per.method" = best.numRej.per.method,
    "relevant.quantile.per.method" = relevant.quantile.per.method,
    "relevant.threshold.per.method" = relevant.threshold.per.method
  )
  return(return.list)

}
