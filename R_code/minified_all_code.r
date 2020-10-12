## start of installation block ##

prepareREnvironment <- function(
        CRANDependencies = c(
            "ggplot2",
            "ggrepel",
            "tidyr",
            "dplyr",
            "reshape2",
            "RColorBrewer"
        ),
        BiocDependencies = c(
            "genefilter",
            "DESeq2",
            "edgeR",
            "pheatmap",
            "regioneR",
            "karyoploteR",
            "fgsea",
            "reactome.db",
            "ToPASeq",
            "org.Hs.eg.db" # see https://bioinformatics.stackexchange.com/questions/5229/converting-gene-symbol-to-ensembl-id-in-r
        )
    ) {

    # doc

    package.check <- lapply(
        CRANDependencies,
        FUN = function(x) {
            if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE, quiet = TRUE)
            library(x, character.only = TRUE)
            }
        }
    )
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
        # version 3.9 to allow compatibility with R 3.6 on Google Colab
        install.packages("BiocManager", version = 3.9)
    
    package.check <- lapply(
        BiocDependencies,
        FUN = function(x) {
            if (!require(x, character.only = TRUE)) {
            BiocManager::install(x)
            library(x, character.only = TRUE)
            }
        }
    )
    return(TRUE)
}

rawCountsMatrix_to_dataframe <- function( rawCountsMatrix_filepath,
                                          sep = ",", 
                                          header = TRUE, 
                                          row.names = 1,
                                          keep_columns = c(),
                                          drop_columns = c(),
                                          encoding = 'utf-8',
                                          make_histogram = FALSE, 
                                          make_boxplot = FALSE,
                                          remove_zero_rows = TRUE,
                                          conditions = function(colnames) sub("\\_.*", "", colnames),
                                          verbose = FALSE) {
  
  # Loads the contents of the `sep`-separated file specified by `rawCountsMatrix_filepath`,
  # keeps the columns specified by `keep_columns` and drops columns specified by `drop_columns`,
  # creates a vector of condition labels inferred from the retained columns' names, and
  # returns the filtered dataframe, the vector of condition labels
  #
  # PARAMETERS:
  #
  #   rawCountsMatrix_filepath: string path to the file containing raw counts.
  #
  #   sep, header, row.names:   (optional) parameters to pass onto read.table.
  #                             default to sep = ",", header = TRUE, row.names = 1.
  #
  #   keep_columns, drop_columns:   (optional) vectors of column names (strings) to keep/drop.
  #                                 cannot contain any shared strings.
  #                                 default to c(), c().
  #
  #   encoding: (optional) string used to specify the encoding of the file specified by `rawCountsMatrix_filepath`.
  #
  #   make_histogram, make_boxplot: (optional) booleans indicates whether to produce exploratory plots or not.
  #                 default to TRUE, TRUE.
  #
  #   verbose: (optional) boolean indicator of verbosity. 
  #            defaults to TRUE.
  #
  # RETURNS:
  #
  #   list("raw.data" = raw.data, "conditions" = conditions, "ENSEMBL_ID_to_Symbol" = ENSEMBL_ID_to_Symbol)
  #
  #     where:  raw.data is the raw counts dataframe after column selection,
  #
  #             conditions is the vector of conditions inferred from the selected columns' names, and
  #
  #             ENSEMBL_ID_to_Symbol is either NULL (if make_ensembl_to_symbol=FALSE) 
  #             or a [ENSEMBL ID, Gene Symbol] dataframe (if make_ensembl_to_symbol=TRUE)
  #
  # TODO:   make conditions a kwarg?
  #         replace gene symbol functionality with a premade package
  
  if(length(intersect(keep_columns, drop_columns)) > 0) {
    stop("`keep_columns` and `drop_columns` need to be mutually exclusive")
  }
  
  if( (length(keep_columns) == 0) & (length(drop_columns) == 0) ) {
    keep_columns <- "all"
    print("since `keep_columns` and `drop_columns` were empty vectors, all columns will be kept.")
  }
  
  if(verbose) {
    print("Loading raw counts matrix into a dataframe")
  }
  
  # load data into dataframe
  raw.data <- data.frame(
    read.csv(
      rawCountsMatrix_filepath, 
      sep=sep, 
      header=header, 
      row.names=row.names,
      encoding=encoding
    )
  )
  
  if(verbose) {  
    # print metadata
    print("Dimensions of raw count dataframe before any filtering:")
    print(dim(raw.data))
    print("Columns of raw count dataframe before any filtering:")
    print(colnames(raw.data))
    
    # print a preview
    print("Head of raw count dataframe before any filtering:")
    print(head(raw.data))
    print("Tail of raw count dataframe before any filtering:")
    print(tail(raw.data))
  }

  # positive column selection
  if(keep_columns != "all"){
    if(verbose) {
      print("Keeping the following columns in raw counts matrix:")
      print(keep_columns)
    }
    raw.data <- raw.data[,which(names(raw.data) %in% keep_columns)]
  }  

  # negative column selection
  if(length(drop_columns) > 0) {
    if(verbose) {
      print("Dropping the following columns in raw counts matrix:")
      print(drop_columns)
    }
    raw.data <- raw.data[,-which(names(raw.data) %in% drop_columns)]
  }
  
  if(remove_zero_rows) {
    zero_rows <- rowSums(raw.data) == 0
    if(verbose) {
        print(paste("Removing ", sum(zero_rows), " genes that had 0 raw reads in all kept samples."))
    }
    raw.data <- raw.data[!zero_rows, ]
    if(verbose) {
        print(paste("Kept ", dim(raw.data)[1], " genes that had > 1 raw reads in any kept sample."))
    }
  }

  if( (typeof(conditions) == "closure") | (class(conditions) == "closure")){
        conditions <- conditions(colnames(raw.data))
  } else if(is.vector(conditions)) {
        same_length <- len(conditions) == len(colnames(raw.data))
        if(!same_length) {
            stop(
                paste(
                    "`conditions` vector is ",
                    length(conditions),
                    " long, but there are ",
                    length(colnames(raw.data)),
                    " columns that need a condition assignment"
                )
            )
        } 
  } else {
        stop("`conditions` must be a vector or a function to apply to the columns of the raw data data.frame")
  }

  if(verbose){
    print("Conditions:")
    print(conditions)
  }
  
  if(make_histogram | make_boxplot) {
    plotLog1PReadCountsDistribution(
      raw.data, 
      conditions = conditions,
      histogram  = make_histogram, 
      boxplot = make_boxplot,
      verbose = verbose
    )
  }
  
  return(
    list(
    "raw.data" = raw.data, 
    "conditions" = conditions
    )
  )

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

# MDS plot
#names(colors)=c(unique(conditions))
#plotMDS(raw.data, labels=colnames(raw.data),   col = colors[conditions])

plotLog1PReadCountsDistribution <- function(raw_reads_dataframe, 
                                            conditions = NULL,
                                            histogram = TRUE, 
                                            boxplot = TRUE,
                                            title = "", 
                                            verbose = FALSE){
  # Function used to plot histogram and box plots from a dataframe.
  #                                              
  # PARAMETERS:
  #
  #   raw_reads_dataframe:  data.frame object with samples as columns and genes as rows.
  #
  #   conditions: (optional) iterable of strings assigning each sample 
  #               (column in `raw_reads_dataframe`) to a condition.
  #               defaults to NA.
  #               e.g. if `raw_reads_dataframe`'s columns represented 
  #               wt_sample1, wt_sample2, cancer_sample1, cancer_sample2,
  #               then `conditions` could be set to = c("wt", "wt", "cancer", "cancer")
  #
  #   histogram, boxplot: (optional) booleans indicating whether to plot a specific graph.
  #                       both default to TRUE.
  #   title:  (optional) string to use as title.
  #           defaults to "" (empty string).
  #
  #   verbose:  (optional) boolean indicator of verbosity.
  #             defaults to FALSE.
  #
  # RETURNS:
  #   
  #   Nothing.


  if(is.null(conditions)) {
    conditions <- sub(
      "\\_.*", "",
      colnames(raw_reads_dataframe)
    )
    if(verbose){
      print("Inferred following conditions from filtered dataframe:\n")
      print(conditions)
    }
  } else {
    stopifnot(length(conditions) == dim(raw_reads_dataframe)[2])
  }

  colors = c(1:length(unique(conditions)))
  
  if(histogram) {
    # is as.data.frame needed?
    readcounts.histogram <- ggplot(gather(as.data.frame(log(cpm(raw_reads_dataframe)))), aes(value)) + 
      geom_histogram(binwidth = 1) + 
      facet_wrap(~key) +
      labs(x = "log(cpm(raw reads mapped to genes))")
    print(readcounts.histogram + ggtitle(title))

    colors = brewer.pal(length(colnames(raw_reads_dataframe)), "Set1")
    plot(
        density(log(cpm(raw_reads_dataframe[,colnames(raw_reads_dataframe)[1]]))),
        col=colors[1], 
        lwd=2, 
        las=2, 
        main="", 
        xlab=""
    )

    title(main="raw data", xlab="log(cpm(raw reads mapped to genes))")
    for (i in 2:dim(raw_reads_dataframe)[2]){
        den <- density(log(cpm(raw_reads_dataframe[,colnames(raw_reads_dataframe)[i]])))
        lines(den$x, den$y, col=colors[i], lwd=2)
    }
    legend("topright", colnames(raw_reads_dataframe), text.col=colors, bty="n")
  }
  
  if(boxplot) {
    # is as.data.frame needed?
    readcounts.boxplot <- ggplot(
      data = reshape2::melt(as.data.frame(log(cpm(raw_reads_dataframe)))), aes(x=variable, y=value)) + 
      geom_boxplot(aes(fill=variable)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs(y = "log(cpm(raw reads mapped to genes))", x = "sample")
    print(readcounts.boxplot + ggtitle(title))
  }

  nonzero_genes <- paste(
    colnames(raw_reads_dataframe),
    apply(raw_reads_dataframe, 2, function(x) sum((x > 0)))
  )

  par(mar=c(12,6,4,1)+.1)
  barplot(
    apply(raw_reads_dataframe, 2, sum), 
    # names.arg = colnames(raw_reads_dataframe),
    names.arg = paste(nonzero_genes, "\ngenes with >= 1 raw read(s)"),
    ylab = "Raw read depth\n(library size)",
    las = 2
  )
  

}

plotPheatMap <- function(data_obj, fromtool = NA, title = "", deseq2_norm = FALSE, ...){
  # Wrapper around `pheatmap` to draw clustered heatmaps.
  #
  # PARAMETERS:
  #
  #   data_obj: one of edgeR's DGEList or DESeq2's DESeqDataSet.
  #
  #   fromtool: string among 'edger', 'deseq' indicating which tool generated `data_obj`
  #
  #   title: (optional) string to use as plot title.
  #
  # RETURNS:
  #   
  #   Nothing.

  stopifnot(!is.na(fromtool))

  if (tolower(fromtool) == "edger") {
    sample.counts <- cpm(data_obj) # cpm(data_obj$counts)
    sample.names <- colnames(data_obj)
  } else if (tolower(fromtool) == "deseq") {
    sample.counts <- counts(data_obj, normalized = deseq2_norm)
    sample.names <- colnames(data_obj)
  } else {
    print("plotPheatMap's `fromtool` keyword argument should be 'edger' or 'deseq'")
    return(NA)
  }

  sample.dists <- cor(sample.counts, method = c("pearson"))
  # round to integer
  sample.dists.vals <- round(sample.dists, digits=3)
  # remove NAs
  sample.dists.vals[ is.na(sample.dists.vals) ] <- ""
  
  rownames(sample.dists) <- sample.names
  colnames(sample.dists) <- sample.names
  colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
  
  pm <- pheatmap(sample.dists,
                 col=colors,
                 main=paste("PheatMap\n(sample correlation)", paste("(",title,")", sep=''), sep='\n'),
                 display_numbers = sample.dists.vals,
                 ...
  )
  print(pm)
}

plotRejectionCurve <- function(filtering.methods, rejections, theta){
  # 
  colors = brewer.pal(ncol(filtering.methods), "Set1")
  matplot(theta, rejections, type='l', lty=1, col=colors, 
          lwd=2, xlab=expression(theta), ylab="Rejections", 
          main=title)
  legend("bottomleft", legend=colnames(filtering.methods), fill=colors)
}

plotChromoMap <- function(DEG_gene_names, reference_file="UCSC_hg38.ncbiRefSeq.tsv") {
  gene.info <- read.table(
    reference_file, 
    stringsAsFactors = FALSE, 
    sep="\t"
  )
  gene.info <- gene.info[!duplicated(gene.info$V6),]
  condensed.gene.info <- data.frame(gene.info$V2, gene.info$V4, gene.info$V5, gene.info$V6)
  names(condensed.gene.info) <- c("chr", "start","end","name")
  genomewide.ranges.for.density <- toGRanges(condensed.gene.info)
  
  chromomap.data <- gene.info[ gene.info$V6 %in% DEG_gene_names, ]
  genes.notin.chromomap <- DEG_gene_names[ !(DEG_gene_names %in% chromomap.data$V6) ]
  
  print("Could not include the following genes in the ChromoMap:")
  print(genes.notin.chromomap)
  
  formatted.chromomap.data <- data.frame(chromomap.data$V2,chromomap.data$V4,chromomap.data$V5,chromomap.data$V6)
  names(formatted.chromomap.data) <- c("chr", "start","end","name")
  
  kp <- plotKaryotype(genome="hg38", plot.type=6)
  kpDataBackground(kp)
  gr <- toGRanges(formatted.chromomap.data)
  kpPlotDensity(kp, data=genomewide.ranges.for.density)
  kpPlotMarkers(kp, data=gr, 
                labels=gr$name, text.orientation = "horizontal",
                r1=0.5, cex=0.6, adjust.label.position = TRUE)
}

make.legend <- function(named_num, quantile_df, threshold_df) {
  stopifnot( all(names(named_num) == colnames(quantile_df)) & all(colnames(threshold_df) == names(named_num)) )
  out <- c()
  for(i in 1:length(named_num)) {
    out <- c(
      out, 
      paste(
        names(named_num)[i], ": ", 
        named_num[i], 
        " DEGs; q: ", 
        quantile_df[[names(named_num)[i]]],
        ", read threshold: ",
        round(threshold_df[[names(named_num)[i]]], 3)
      )
    )
  }
  return(out)
}

plotDESeq2RejectionCurve <- function(DESeq2_result_obj, title, ...) {
    #
    #stopifnot(is(DESeq2_result_obj, result))
    plot(
        metadata(DESeq2_result_obj)$filterNumRej,
        type = "b", 
        ylab = "number of rejections (DEGs)",
        xlab = "quantiles of filter", 
        main = title
    )
    lines(
        metadata(DESeq2_result_obj)$lo.fit, 
        col="red"
    )
    abline(v=metadata(DESeq2_result_obj)$filterTheta)   
}

convert.names <- function(df.with.rownames,
                          inputType,
                          outputType,
                          regex_to_replace = ".[0-9]+$") {
  # 
  inputType_ids <- sub(regex_to_replace, "", rownames(df.with.rownames))
  out <- mapIds(
    org.Hs.eg.db, 
    keys = inputType_ids,
    column = outputType,
    keytype = inputType,
    multiVals = "first"
  )
  return(out)
}

DESeq2.pairwise.explo <- function(DESeq2_dataset,
                            alpha = 0.1, 
                            test = "Wald", 
                            showplots = TRUE, 
                            independent.filtering.statistics = NA,
                            quantiles = seq(from = 0.0, to = 0.99, by = 0.025),
                            chosen_filter = "mean",
                            chosen_filter_threshold = NULL,
                            pAdjustMethod = "BH",
                            lfcThreshold = 0,
                            verbose = FALSE, 
                            genes_of_interest = NA,
                            ...) {

  # docs
  
  # argument validation and sanity checks
  stopifnot(is(DESeq2_dataset, "DESeqDataSet"))
  stopifnot( (0 <= alpha) & (alpha <= 1) )
  test <- match.arg(test, choices = c("Wald", "LRT"))
  
  # prepare for independent.filtering.statistics
  DESeq2_dataset <- DESeq(
    DESeq2_dataset, 
    test = test, 
    ...
  )
  
  # continue with argument validation and sanity-checks
  if(is.na(independent.filtering.statistics)) {
    independent.filtering.statistics <- data.frame(
      'mean' = rowMeans(counts(DESeq2_dataset, normalized = TRUE)),
      'median' = rowMedians(counts(DESeq2_dataset, normalized = TRUE))
    )
  } else {
    stopifnot(class(independent.filtering.statistics) == "data.frame")
  }
  quantiles <- sort(quantiles)
  stopifnot( (min(quantiles) >= 0 ) & (max(quantiles) <= 1.0 ) )
  stopifnot(chosen_filter %in% colnames(independent.filtering.statistics))
  if(!is.null(chosen_filter_threshold)) {
    stopifnot(
      ( class(chosen_filter_threshold) == "function" ) |
      ( typeof(chosen_filter_threshold) == "closure")
    )
  }
  if(!is.na(genes_of_interest)) {
    missing <- !(genes_of_interest %in% names(DESeq2_dataset))
    if(any(missing)) {
      stop(
        paste0(
          "the following genes of interest were not found in the dataset:\n",
          genes_of_interest[missing]
        )
      )
    } else {
      goi_survival <- data.frame(
        matrix(
          nrow = length(genes_of_interest),
          ncol = dim(independent.filtering.statistics)[2]
        )
      )
      rownames(goi_survival) <- genes_of_interest
      colnames(goi_survival) <- colnames(independent.filtering.statistics)
      for(i in 1:dim(independent.filtering.statistics)[2]){
        statistic <- colnames(independent.filtering.statistics)[i]
        filter.statistic.column <- independent.filtering.statistics[,i]
        expression_percentile <- ecdf(sort(filter.statistic.column))
        goi_filter_statistic_survival <- expression_percentile(
          independent.filtering.statistics[genes_of_interest,i]
        )

        goi_survival[genes_of_interest,statistic] = goi_filter_statistic_survival
        
      }
      print(goi_survival)
    }
  }
  
  # run standard DESeq2 analysis
  DESeq2.default.result <- results(
    DESeq2_dataset,
    alpha = alpha,
    theta = quantiles, 
    pAdjustMethod = pAdjustMethod,
    lfcThreshold = lfcThreshold,
    independentFiltering = TRUE
  )
  if(showplots){
    plotDESeq2RejectionCurve(
      DESeq2.default.result,
      paste0(
          "Rejection curve for default DESeq2",
          "\n(alpha (a.k.a. FDR) = ",
          alpha,
          ", test = ", 
          test,
          ")\n",
          sum(DESeq2.default.result['padj'][[1]] < alpha, na.rm=TRUE),
          " DEGs when ignoring genes with mean normalized expression < ",
          round(metadata(DESeq2.default.result)$filterThreshold,2)
        )
    )
  }
  
  # run customized DESeq2 analysis if specified
  if((class(chosen_filter_threshold) == "function") | (typeof(chosen_filter_threshold) == "closure")) {
    filter <- chosen_filter_threshold( independent.filtering.statistics[[chosen_filter]] )
  } else if(is.null(chosen_filter_threshold)) {
    filter <- independent.filtering.statistics[[chosen_filter]]
  } else {
    stop("chosen_filter_threshold should be a function.")
  }
  
  DESeq2.specified.results <- results(
    DESeq2_dataset, 
    alpha = alpha,
    theta = quantiles, 
    pAdjustMethod = pAdjustMethod,
    lfcThreshold = lfcThreshold,
    independentFiltering = TRUE,
    filter = filter
  )
  if(showplots){
    plotDESeq2RejectionCurve(
      DESeq2.specified.results,
      paste0(
        "Rejection curve for DESeq2 with specified IF parameters",
        "\n(alpha (a.k.a. FDR) = ",
        alpha,
        ", test = ", 
        test,
        ")\n",
        sum(DESeq2.specified.results['padj'][[1]] < alpha, na.rm=TRUE),
        " DEGs"
      )
    )

    #un.NA.pvalues <- DESeq2.specified.results[['pvalue']]
    #un.NA.pvalues[is.na(un.NA.pvalues)] <- 1.0
    #un.NA.padj <- DESeq2.specified.results[['padj']]
    #un.NA.padj[is.na(un.NA.padj)] <- 1.0

    #DESeq2.specified.results.df <- data.frame(
    #  chosen.filter.stat = filter, # y
    #  log2FoldChange = DESeq2.specified.results[['log2FoldChange']], # x
    #  alpha = 1.-un.NA.pvalues, # alpha
    #  un.NA.padj = un.NA.padj
    #)
    
    #DESeq2.specified.results.df$label <- 1:nrow(DESeq2.specified.results.df)
    #DESeq2.specified.results.df$label[ 
    #  DESeq2.specified.results.df[['chosen.filter.stat']] < 0.25*max(DESeq2.specified.results.df[['chosen.filter.stat']]) 
    #] <- NA
    
    #v2plot <- ggplot(DESeq2.specified.results.df, 
    #  aes(x=log2FoldChange, chosen.filter.stat, label = label)) +
    #  geom_point(aes(alpha = alpha)) +
    #  theme_minimal() + 
    #  geom_text_repel() 
    #print(v2plot)

    #DESeq2.specified.results.df$label <- 1:nrow(DESeq2.specified.results.df)
    #DESeq2.specified.results.df$label[ abs(DESeq2.specified.results.df[['log2FoldChange']]) < 0.6 ] <- NA
    #vplot <- ggplot(data=DESeq2.specified.results.df, 
    #  aes(x=log2FoldChange, y=-log10(un.NA.padj), label = label)) +
    #    geom_point() + 
    #    theme_minimal() +
    #   geom_text_repel() +
    #    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    #    geom_hline(yintercept=-log10(alpha), col="red")
    #print(vplot)
  }
  
  # independent filtering section

  rejections <- sapply(
    independent.filtering.statistics,
    function(filter_statistic) filtered_R(
      alpha = alpha, 
      filter = filter_statistic,
      test = DESeq2.default.result$pvalue,
      theta = quantiles,
      method = pAdjustMethod
    )
  )
  
  best.numRej.per.method <- apply(rejections, 2, max, na.rm = TRUE)
  relevant.quantile.per.method <- list()
  relevant.threshold.per.method <- list()
  best.results.per.method <- list()
  
  for(i in 1:length(independent.filtering.statistics)) {
    filter_name <- colnames(independent.filtering.statistics)[i]
    # the extra [1] is to get the lowest quantile that yielded the max number of rejections
    # to filter out fewer genes
    q <- quantiles[ which(rejections[,i] == max(rejections[,i])) ][1]
    relevant.quantile.per.method[filter_name] = q
    relevant.column <- independent.filtering.statistics[,filter_name]
    relevant.threshold.per.method[filter_name] = quantile(relevant.column,q)
    best.results.per.method[filter_name] = results(
      DESeq2_dataset,
      alpha = alpha,
      theta = quantiles, 
      pAdjustMethod = pAdjustMethod,
      lfcThreshold = lfcThreshold,
      independentFiltering = TRUE,
      filter = relevant.column
    )
  }
  
  # plot rejection curves for statistics in independent.filtering.statistics
  if(showplots) {
    colors=brewer.pal(ncol(independent.filtering.statistics), "Set1")
    matplot(
      quantiles, 
      rejections, 
      type = 'l', 
      lty = 1, 
      col = colors, 
      lwd = 2, 
      xlab = paste0(expression(theta), " (filter's percentile threshold)"),
      ylab = "Rejections (number of DEGs)", 
      main = paste0(
        "Rejection curves for DESeq2",
        "\nusing different statistic filters"
      )
    )
    legend(
      "bottomleft",
      legend = make.legend(best.numRej.per.method, relevant.quantile.per.method, relevant.threshold.per.method),
      fill = colors
    )
  }
  
  return(
    list(
      data = DESeq2_dataset,
      DESeq2.default.result = DESeq2.default.result,
      DESeq2.specified.results = DESeq2.specified.results,
      best.numRej.per.method = best.numRej.per.method,
      relevant.quantile.per.method = relevant.quantile.per.method,
      relevant.threshold.per.method = relevant.threshold.per.method,
      best.results.per.method = best.results.per.method
    )
  )
}

edgeR.pairwise.explo <- function( edgeR_dataset,
                                  design.matrix,
                                  alpha = 0.1, 
                                  test = "exactTest", 
                                  showplots = TRUE, 
                                  independent.filtering.statistics = NA,
                                  quantiles = seq(from = 0.0, to = 1.0, by = 0.025),
                                  chosen_filter = "mean",
                                  chosen_filter_threshold = NULL,
                                  pAdjustMethod = "BH",
                                  verbose = FALSE, 
                                  genes_of_interest = NA,
                                  ...) {
  # docs
  
  # argument validation and sanity checks
  stopifnot(is(edgeR_dataset, "DGEList"))
  stopifnot( (0 <= alpha) & (alpha <= 1) )
  test <- match.arg(test, choices = c("exactTest", "LRT", "QLF"))
  
  if(is.na(independent.filtering.statistics)) {
    independent.filtering.statistics <- data.frame(
      'mean' = rowMeans(cpm(edgeR_dataset$counts)),
      'median' = rowMedians(cpm(edgeR_dataset$counts))
    )
  } else {
    stopifnot(class(independent.filtering.statistics) == "data.frame")
  }

  quantiles <- sort(quantiles)
  stopifnot( (min(quantiles) >= 0 ) & (max(quantiles) <= 1.0 ) )

  sorted.independent.filtering.statistics <- as.data.frame(sapply(independent.filtering.statistics, sort))
  quantiled.independent.filtering.statistics <- sorted.independent.filtering.statistics[
    dim(sorted.independent.filtering.statistics)[1] * quantiles, 
  ]
  
  stopifnot(chosen_filter %in% colnames(independent.filtering.statistics))
  if((class(chosen_filter_threshold) == "function") | (typeof(chosen_filter_threshold) == "closure")) {
    filter <- chosen_filter_threshold( independent.filtering.statistics[[chosen_filter]] )
  } else if(is.null(chosen_filter_threshold)) {
    filter <- independent.filtering.statistics[[chosen_filter]]
  } else {
    stop("chosen_filter_threshold should be a function.")
  }
  if(!is.na(genes_of_interest)) {
    missing <- !(genes_of_interest %in% names(DESeq2_dataset))
    if(any(missing)) {
      stop(
        paste0(
          "the following genes of interest were not found in the dataset:\n",
          genes_of_interest[missing]
        )
      )
    } else {
      goi_survival <- data.frame(
        matrix(
          nrow = length(genes_of_interest),
          ncol = dim(independent.filtering.statistics)[2]
        )
      )
      rownames(goi_survival) <- genes_of_interest
      colnames(goi_survival) <- colnames(independent.filtering.statistics)
      for(i in 1:dim(independent.filtering.statistics)[2]){
        statistic <- colnames(independent.filtering.statistics)[i]
        filter.statistic.column <- independent.filtering.statistics[,i]
        expression_percentile <- ecdf(sort(filter.statistic.column))
        goi_filter_statistic_survival <- expression_percentile(
          independent.filtering.statistics[genes_of_interest,i]
        )

        goi_survival[genes_of_interest,statistic] = goi_filter_statistic_survival
        
      }
      print(goi_survival)
    }
  }
  
  # run edgeR analysis with independent filtering
  best.numRej.per.method <- list()
  relevant.quantile.per.method <- list()
  relevant.threshold.per.method <- list()
  best.results.per.method <- list()
  numRej.per.method.and.quantile <- data.frame(
    matrix(0L, nrow = length(quantiles), ncol = dim(independent.filtering.statistics)[2])
  )
  colnames(numRej.per.method.and.quantile) <- colnames(independent.filtering.statistics)

  for(m in 1:dim(independent.filtering.statistics)[2]) {
    filter_name <- colnames(independent.filtering.statistics)[m]
    best.numRej.so.far <- 0
    relevant.quantile.so.far <- NA
    relevant.threshold.so.far <- NA
    best.result.so.far <- NA
    
    relevant.statistic.column <- independent.filtering.statistics[,m]
    print(paste("Running edgeR with the following filter statistic: ", filter_name))
    progress.bar <- txtProgressBar(min = 0, max = length(quantiles), style = 3)
    for(q in 1:length(quantiles)) {
      setTxtProgressBar(progress.bar, q)
      keep <- relevant.statistic.column > quantiled.independent.filtering.statistics[q,m]
      if(sum(keep) == 0){
        break
      }
      if(verbose) {
        mess <- paste(filter_name, " ", quantiles[q], " ", q, " ", as.matrix(quantiled.independent.filtering.statistics)[q,m], " ", sum(keep))
        print(mess)
      }
      filtered_edgeR_dataset <- edgeR_dataset[ keep, , keep.lib.sizes = FALSE ]
      filtered_edgeR_dataset <- estimateDisp(filtered_edgeR_dataset, design = design.matrix)
      
      if(tolower(test) == "exacttest") {
        test.out <- exactTest(filtered_edgeR_dataset)
      } else if(tolower(test) == "lrt") {
        fit <- glmFit(filtered_edgeR_dataset, design.matrix)
        test.out <- glmLRT(fit, coef=coef)
      } else if(tolower(test) == "qlf") {
        fit <- glmQLFit(filtered_edgeR_dataset, design.matrix)
        test.out <- glmQLFTest(fit, coef=coef)
      } else {
        stop(paste("test must be one of 'exactTest', 'LRT', or 'QLF', not ", test))
      }
      
      filtered_edgeR_results <- topTags(
        test.out,
        n = dim(edgeR_dataset)[1],
        adjust.method = pAdjustMethod,
        p.value = alpha
      )
      nRej <- dim(filtered_edgeR_results)[1]
      numRej.per.method.and.quantile[q, filter_name] <- nRej
      if(nRej > best.numRej.so.far) {
        best.numRej.so.far <- nRej
        relevant.quantile.so.far <- quantiles[q]
        relevant.threshold.so.far <- quantiled.independent.filtering.statistics[q,m]
        best.result.so.far <- filtered_edgeR_results
      }
      
      
    }
    close(progress.bar)
    best.numRej.per.method[filter_name] = best.numRej.so.far
    relevant.quantile.per.method[filter_name] = relevant.quantile.so.far
    relevant.threshold.per.method[filter_name] = relevant.threshold.so.far
    best.results.per.method[filter_name] = best.result.so.far

    row.names(numRej.per.method.and.quantile) <- quantiles
  }

  colors = brewer.pal(ncol(numRej.per.method.and.quantile), "Set1")
  matplot(
    numRej.per.method.and.quantile, 
    type = 'l', 
    lty = 1, 
    col = colors, 
    lwd = 2, 
    ylim = c( min(numRej.per.method.and.quantile) - 1, max(numRej.per.method.and.quantile) + 1 ),
    xlab = paste0(expression(theta), " (filter's percentile threshold)"),, 
    ylab = "Rejections", 
    main = "title"
  )
  legend("bottomleft", legend = colnames(numRej.per.method.and.quantile), fill = colors)
  
  return(
    list(
      best.numRej.per.method = best.numRej.per.method,
      relevant.quantile.per.method = relevant.quantile.per.method,
      relevant.threshold.per.method = relevant.threshold.per.method,
      best.results.per.method = best.results.per.method
    )
  )
  
  
}

run_fgsea <- function(  named_vector, 
                        pathway_data = NULL, 
                        fdr = 0.1, 
                        minsize = 10, 
                        maxsize = 100000,
                        nperms = 10000,
                        seed = 42 ) {
    #

    set.seed(seed)

    if(is.null(pathway_data)) {
        pathways <- reactomePathways(names(named_vector))    
        print("defaulting to using reactomePathways since `pathway_data` was NULL")
        print(pathways)
    } else {
        pathways <- gmtPathways(pathway_data)
    }


    fgsea.out <- fgsea::fgsea(
        pathways = pathways, 
        stats = named_vector,
        minSize = minsize,
        maxSize = maxsize,
        nperm = nperms) %>% 
        as.data.frame() %>% 
        dplyr::filter(padj < !!fdr)
    
    overrep <- fgsea.out[ fgsea.out[['NES']] > 0 ]
    underrep <- fgsea.out[ fgsea.out[['NES']] < 0 ]

    fgRes <- fgsea.out[ !is.na(match(fgsea.out$pathway, c( overrep$mainPathways, underrep$mainPathways))), ] %>% 
        arrange(desc(NES))
    fgRes$pathway <- stringr::str_replace(fgRes$pathway, "GO_" , "")
    
    fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
    filtRes <- rbind(head(fgRes, n = 10),
                    tail(fgRes, n = 10 ))
    #g <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    #    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    #    geom_point( size=5, aes( fill = Enrichment),
    #            shape=21, stroke=2) +
    #    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
    #                    "Up-regulated" = "firebrick") ) +
    #    coord_flip() +
    #    labs(x="Pathway", y="Normalized Enrichment Score",
    #        title="GSEA - Biological Processes") + 
    #    theme_minimal()
    #print(g)
    #output = list("Results" = fgRes, "Plot" = g)
    output = list("Results" = fgRes)
    return(output)
}