###############
# DESCRIPTION #
###############


## start of installation block ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.

## First specify the packages of interest from both BiocManager and CRAN
Bioc_dependencies = c(
  "DESeq2", 
  "HTSFilter",
  "genefilter",
  "edgeR",
  "limma",
  "Glimma",
  "pheatmap",
  "biomaRt",
  "regioneR",
  "chromoMap",
  "karyoploteR"
)

CRAN_dependencies = c(
  "data.table",
  "tidyr",
  "plyr",
  "dplyr",
  "reshape2",
  "factoextra",
  "RColorBrewer",
  "ggplot2",
  "ggpubr"
)

## Now load or install & load all
package.check <- lapply(
  Bioc_dependencies,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x)
      library(x, character.only = TRUE)
    }
  }
)

package.check <- lapply(
  CRAN_dependencies,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, quiet = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

### end of installation block ###


## start of file_io.R ##
# external dependencies
library(biomaRt)
library(stringr)
library(dplyr)

rawCountsMatrix_to_dataframe <- function( rawCountsMatrix_filepath,
                                          sep = ",", 
                                          header = TRUE, 
                                          row.names = 1,
                                          keep_columns = c(),
                                          drop_columns = c(),
                                          encoding = 'utf-8',
                                          make_ensembl_to_symbol = FALSE,
                                          gene_symbol_column = NA, 
                                          make_histogram = FALSE, 
                                          make_boxplot = FALSE,
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
  #   make_ensembl_to_symbol:   (optional) if TRUE, will create a [ENSEMBL ID, Gene Symbol] 
  #                             dataframe and include it in the returned list. 
  #                             Assumes that the ENSEMBL IDs are the row names of the
  #                             `raw.data` dataframe.
  #                             if FALSE, the third item in the returned list will be NA.
  #                             defaults to FALSE.
  #
  #   gene_symbol_column: (optional) string name of the column containing gene symbols.
  #                       defaults to NA.
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
  # TODO: make conditions a kwarg?
  
  if(length(intersect(keep_columns, drop_columns)) > 0) {
    stop("`keep_columns` and `drop_columns` aren't mutually exclusive")
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
    print("Dimensions of raw count dataframe:")
    print(dim(raw.data))
    print("Columns of raw count dataframe:")
    print(colnames(raw.data))
    
    # print a preview
    print("Head of raw count dataframe:")
    print(head(raw.data))
    print("Tail of raw count dataframe:")
    print(tail(raw.data))
  }
  
  if(make_ensembl_to_symbol) {
    if(verbose) print("Creating ENSEMBL ID <-> Gene Symbol dataframe")
    # create [ENSEMBLE ID, Gene Symbol] dataframe for future lookups
    stopifnot(!is.na(gene_symbol_column))
    stopifnot(gene_symbol_column %in% colnames(raw.data))
    ENSEMBL_ID_to_Symbol <- data.frame(
      rownames(raw.data), raw.data[gene_symbol_column]
    )
    colnames(ENSEMBL_ID_to_Symbol) <- c("ENSEMBL ID", gene_symbol_column)
    
  } else {
    make_ensembl_to_symbol=FALSE
  }
  
  # positive column selection
  if(all(keep_columns != "all")){
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
  
  conditions <- sub(
    "\\_.*", "",
    colnames(raw.data)
  )
  if(verbose){
    print("Inferred following conditions from filtered count matrix:\n")
    print(conditions)
  }
  
  if(make_histogram | make_boxplot) {
    plotLog1PReadCountsDistribution(
      raw.data, 
      histogram=make_histogram, 
      boxplot=make_boxplot,
      verbose=verbose
    )
  }
  
  if(make_ensembl_to_symbol) {
    return(
      list(
        "raw.data" = raw.data, 
        "conditions" = conditions, 
        "ENSEMBL_ID_to_Symbol" = ENSEMBL_ID_to_Symbol
      )
    )
  } else {
    return(
      list(
        "raw.data" = raw.data, 
        "conditions" = conditions,
        "ENSEMBL_ID_to_Symbol" = NA
      )
    )
  }
}

# to remove
write.gct <- function(gct, filename, check.file.extension=TRUE){
  # Save a representation of `gct` to a file, ensuring the filename has the extension .gct
  # in order to produce the DESeq2-normalized counts file (.gct) to use as input
  # to GSEA.
  # shamelessly ripped off of:
  # https://github.com/genepattern/DESeq2/blob/abac851160ed1f84af2aa06f85b497c4b144656e/src/common.R
  # 
  # PARAMETERS:
  #   
  #   gct: 
  #
  #   filename: string indicating the path of the output file. 
  #             its extension will be replaced with .gct if `check.file.extension`                
  #             is TRUE.
  #
  #   check.file.extension: boolean indicating whether to replace the extension of
  #                         `filename` with `.gct`.
  #                         defaults to TRUE.
  # RETURNS:
  #
  #   filename: string path to the output file.
  
  if(check.file.extension) {
    filename <- check.extension(filename, ".gct") 
  }
  if(is.null(gct$data)) {
    exit("No data given.")
  }
  if(is.null(row.names(gct$data))) {
    exit("No row names given.")
  }
  if(is.null(colnames(gct$data))) {
    exit("No column names given.")
  }
  
  rows <- dim(gct$data)[1]
  columns <- dim(gct$data)[2]
  
  if(rows!=length(row.names(gct$data))) {
    exit("Number of data rows (", rows, ") not equal to number of row names (", length(row.names(gct$data)), ").")
  }
  if(columns!=length(colnames(gct$data))) {
    exit("Number of data columns (", columns , " not equal to number of column names (", length(colnames(gct$data)), ").")
  }
  
  if(!is.null(gct$row.descriptions)) {
    if(length(gct$row.descriptions)!=rows) {
      exit("Number of row descriptions (", length(gct$row.descriptions), ") not equal to number of row names (", rows, ").")
    }
  }
  
  row.descriptions <- gct$row.descriptions
  if(is.null(row.descriptions)) {
    row.descriptions <- ''
  }
  m <- cbind(row.names(gct$data), row.descriptions, gct$data)
  f <- file(filename, "w")
  on.exit(close(f))
  
  cat("#1.2", "\n", file=f, append=TRUE, sep="")
  cat(rows, "\t", columns, "\n", file=f, append=TRUE, sep="")
  cat("Name", "\t", file=f, append=TRUE, sep="")
  cat("Description", file=f, append=TRUE, sep="")
  names <- colnames(gct$data)
  
  for(j in 1:length(names)) {
    cat("\t", names[j], file=f, append=TRUE, sep="")
  }
  
  cat("\n", file=f, append=TRUE, sep="")
  write.table(m, file=f, append=TRUE, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
  return(filename)
}

# to remove
write.factor.to.cls <- function(factor, filename) {
  # Save a representation of `factor` to a file which could be used as input
  # to GSEA.
  # shamelessly ripped off of:
  # https://github.com/genepattern/DESeq2/blob/abac851160ed1f84af2aa06f85b497c4b144656e/src/common.R
  # 
  # PARAMETERS:
  #   
  #   factor: 
  #
  #   filename: string indicating the path of the output file. 
  #
  # RETURNS:
  #
  #   filename: string path to the output file.
  
  file <- file(filename, "w")
  on.exit(close(file))
  codes <- unclass(factor)
  cat(file=file, length(codes), length(levels(factor)), "1\n")
  
  levels <- levels(factor)
  
  cat(file=file, "# ")
  num.levels <- length(levels)
  
  if(num.levels-1 != 0)
  {
    for(i in 1:(num.levels-1))
    {
      cat(file=file, levels[i])
      cat(file=file, " ")
    }
  }
  cat(file=file, levels[num.levels])
  cat(file=file, "\n")
  
  num.samples <- length(codes)
  if(num.samples-1 != 0)
  {
    for(i in 1:(num.samples-1))
    {
      cat(file=file, codes[i]-1)
      cat(file=file, " ")
    }
  }
  
  cat(file=file, codes[num.samples]-1)
  return(filename) 
}

get.conversion.map <- function( raw.counts.data, 
                                biomart="ensembl", 
                                dataset="hsapiens_gene_ensembl",
                                BMattributes=c('ensembl_gene_id','external_gene_name'),
                                BMfilters='ensembl_gene_id',
                                str_replace_pattern=".[0-9]+$",
                                str_replace_replacement="",
                                prefetched_file_name="ensembl_ID_gene_symbol_map_may2020.tsv"){
  #
  #
  #
  print("Using pre-computed conversion file instead...")
    
    x.to.y.map <- as.data.frame(
      read.csv(
        prefetched_file_name,
        sep='\t',
        header=TRUE,
        encoding='utf-8',
        row.names = 1
      )
    )
    message("Loaded map!")
    print(head(x.to.y.map))
    return(x.to.y.map)
  
  # use of biomart is deprecated; getBM returns nothing.
  mart <- useMart(
    biomart=biomart,
    dataset=dataset
  )
  
  x.to.y.map <- tryCatch({
    getBM(
      attributes=BMattributes,
      filters=BMfilters,
      values=str_replace(
        rownames(raw.counts.data$raw.data),
        pattern=str_replace_pattern,
        replacement=str_replace_replacement
      ),
      mart=mart
    )
  }, error=function(exception) {
    message("Caught server-side error:")
    message(exception)
    message("Using pre-computed conversion file instead...")
    
    x.to.y.map <- as.data.frame(
      read.csv(
        prefetched_file_name,
        sep='\t',
        header=TRUE,
        encoding='utf-8',
        row.names = 1
      )
    )
    message("Loaded map!")
    print(head(x.to.y.map))
    return(x.to.y.map)
  })
  print("got map from biomaRt")
  print(head(x.to.y.map))
  rownames(x.to.y.map) <- x.to.y.map$ensembl_gene_id
  print(head(x.to.y.map))
  x.to.y.map<-subset(x.to.y.map, select=('external_gene_name'))
  return(x.to.y.map)
}

save.spreadsheet <- function( DESeq2.results, 
                              edgeR.results, 
                              raw.counts.data, 
                              contrasts, 
                              symbol.to.id.map, 
                              fdr,
                              file_prefix, 
                              mainDir, 
                              subDir, 
                              saveoutput=TRUE, 
                              chosen_filter="mean"){
  #
  #
  #
  print("in save.spreadsheet")
  print(head(symbol.to.id.map))
  if(dir.exists(file.path(mainDir, subDir))) {
    print(
      paste0(
        "Directory ",
        file.path(mainDir, subDir),
        " already exists.",
        "Please change your output directory name"
      )
    )
    return()
  } else {
    dir.create(file.path(mainDir, subDir), recursive = TRUE)
  }
  
  all.dataframes <- list()
  for (contrast.name in colnames(contrasts)) {
    print(contrast.name)
    contrast.conditions <- strsplit(contrast.name, "_")[[1]]
    print(contrast.conditions)
    contrast.conditions <- c(contrast.conditions[2], contrast.conditions[4])
    print(contrast.conditions)
    relevant.columns <- lapply(contrast.conditions, grepl, colnames(raw.counts.data$raw.data))
    relevant.columns.mask <- sapply(transpose(relevant.columns), any) # column-wise OR
    relevant.columns <- colnames(raw.counts.data$raw.data)[relevant.columns.mask] # from indices to names
    
    relevant.raw.data <- raw.counts.data$raw.data[relevant.columns]
    
    relevant.DESeq2.normalized.data <- as.data.frame(counts(DESeq2.results$DGE_obj, normalized=TRUE))[relevant.columns.mask]
    colnames(relevant.DESeq2.normalized.data) <- paste("DESeq_norm", colnames(relevant.DESeq2.normalized.data), sep = "_")
    
    relevant.edger.normalized.data <- as.data.frame(t(
      t(edgeR.results$DGE_obj$pseudo.counts)*(edgeR.results$DGE_obj$samples$norm.factors)
    ))[relevant.columns.mask]
    colnames(relevant.edger.normalized.data) <- paste("edgeR_norm", colnames(relevant.edger.normalized.data), sep = "_")

    relevant.DESeq2.results <- DESeq2.results[[contrast.name]][['IF.results']]
    relevant.DESeq2.metadata <- DESeq2.results[[contrast.name]][['IF.metadata']]
    
    relevant.edger.results <- edgeR.results[[contrast.name]][['IF.results']]
    relevant.edger.metadata <- edgeR.results[[contrast.name]][['IF.metadata']]
    
    all_results <- list(
      as.data.frame(relevant.raw.data),
      as.data.frame(relevant.DESeq2.normalized.data),
      as.data.frame(relevant.edger.normalized.data),
      as.data.frame(relevant.DESeq2.results),
      as.data.frame(relevant.edger.results)
    )

    for(i in 1:length(all_results)){
      all_results[[i]]$ROWNAMES  <- rownames(all_results[[i]])
    }
    
    all_results <- join_all(all_results, by="ROWNAMES", type="full")
    
   
    # make a column that will (eventually) store Gene Symbols
    all_results["ID1"] <- all_results["ROWNAMES"]
    all_results["ID2"] <- all_results["ROWNAMES"]
    
    # re-order columns s.t. ID and ROWNAMES is the first two columns and the other columns' order
    # doesn't change
    all_results <- all_results %>% dplyr::select(ID1, ID2, ROWNAMES, everything())
    
    # trim the Ensembl IDs in the ROWNAMES column to remove the .[0-9]+$
    # as this would prevent conversion to gene symbols:
    # see : https://www.biostars.org/p/302441/
    all_results$ID2 <- str_replace(all_results$ROWNAMES, pattern="\.[0-9]+$", replacement = "")
    print("all_results_ID2")
    print(all_results$ID2[1:6])

    all_results$ROWNAMES <- symbol.to.id.map[c(all_results$ID2),]
    print("ROWNAMES")
    print(all_results$ROWNAMES[1:6])
    print(typeof(all_results))
    print(class(all_results))
    print(length(all_results))
    print(head(all_results))
    colnames(all_results)[1] <- "Ensembl ID"
    colnames(all_results)[2] <- "Version-less Ensembl ID"
    colnames(all_results)[3] <- "external_gene_name"
    
    all.dataframes[[contrast.name]] <- all_results
    
    if(saveoutput) {
      # needs sanity-checks!
      output.file.name<-paste0(
        "DGEA_results_",
        file_prefix,
        "_",
        contrast.name,
        "_fdr_",
        fdr,
        ".tsv"
      )
      
      output.file.path <- file.path(
        mainDir,
        subDir,
        output.file.name
      )
      
      if(file.exists(output.file.path)) {
        print(
          paste0(
            "File\n",
            output.file.path,
            "\nalready exists.\n",
            "Please change your output directory name"
          )
        )
        return()
      } else {
        print(paste0("Saving in: ", output.file.path)  )
      }
      
      
      write.table(
        all_results, 
        file=output.file.path, 
        quote=FALSE, 
        sep='\t', 
        row.names=FALSE
      )
      
      output.metadata<-paste0(
        "#metadata",
        "\n#alpha/fdr: ", fdr, 
        "\n#chosen filter: ", chosen_filter, 
        "\n#input file: ", counts.matrix.filepath, 
        "\n#date: ", date(),
        "\n#DESeq2 version: ", packageVersion("DESeq2"),
        "\n#edgeR version: ", packageVersion("edgeR"),
        "\n#biomaRt version: ", packageVersion("biomaRt"),
        sep=""
      )
      
      write(
        output.metadata[1], 
        output.file.path, 
        append=TRUE
      )
      
      print(
        paste0(
          "saved results for ", 
          contrast.name, 
          " in:", 
          output.file.path
        )
      )
    }
  }
  if(saveoutput){
    return(file.path(mainDir,subDir))
  }
}

### end of file_io.R ###


## start of dge_analysis.R ##
library(data.table)
library(HTSFilter)
library(tidyr)
library(plyr)
library(genefilter)
library(edgeR)
library(limma)
library(Glimma)
library(DESeq2)

edgeR_DGE_analysis <- function( DGE_obj, design, alpha, 
                                coef = NULL, 
                                contrasts = NULL,
                                useLRT = FALSE, 
                                useQLF = FALSE, 
                                useEXACT = TRUE,
                                plotRejCurve = TRUE, 
                                filtering.methods = NA, 
                                chosen_filter = NA,
                                quantiles = seq(from=0.0, to=0.99, by=0.025), 
                                pAdjustMethod = "BH", 
                                verbose = TRUE) {
  
  if(!any(useEXACT, useLRT, useQLF)) {
    stop("at least one of `useEXACT`, `useLRT`, or `useQLF` must be TRUE")
  }

  stopifnot(length(quantiles)>0)

  quantiles <- sort(quantiles)

  if(is.na(filtering.methods)) {
    filtering.methods.dataframe <- data.frame(
      'mean'=rowMeans(DGE_obj$counts),
      'min'= rowMin(DGE_obj$counts),
      'max'= rowMax(DGE_obj$counts),
      'median'= rowMedians(DGE_obj$counts),
      'secondlargest'= apply(DGE_obj$counts, 1,  function(row) sort(row, partial=length(row)-1)[length(row)-1])
    )
  } else {
    filtering.methods.dataframe <- filtering.methods
  }

  DGE_obj <- calcNormFactors(DGE_obj)
  DGE_obj <- estimateCommonDisp(DGE_obj)
  DGE_obj <- estimateTagwiseDisp(DGE_obj)
  DGE_obj <- estimateGLMCommonDisp(DGE_obj, design)
  DGE_obj <- estimateGLMTrendedDisp(DGE_obj, design)
  DGE_obj <- estimateGLMTagwiseDisp(DGE_obj, design)

  contrastwise.output.list <- list('DGE_obj'=DGE_obj)

  for(i in 1:dim(contrasts)[2]) {
    edger.result <- edgeR_DGE(
      DGE_obj, design, contrast=contrasts[,i], 
      useLRT = useLRT, useQLF = useQLF, useEXACT = useEXACT
    )
    contrast.name <- colnames(contrasts)[i]
    if(verbose){ 
      print(paste("Ran initial edgeR analysis with contrast = ", contrast.name)) 
    }
    
    # save a copy of the edger results obtained without any filtering
    contrastwise.output.list[[contrast.name]] <- list(
      "naive.edgeR.results"=edger.result
    )
    
    # preparation for independent filtering
    edger.results.copy <- edger.result$edgeR.exactTest
    edger.results.copy$unadjPvalues <- edger.result$edgeR.exactTest$table$PValue
    edger.results.copy$counts <- DGE_obj$counts
    title <- paste0(
      "edgeR(", 
      contrast.name, 
      ") rejection curve (alpha; a.k.a. FDR)=",
      alpha,
      ")"
    )
    
    edger.exactTest.IF.results <- independent_filtering(
      edger.results.copy,
      filtering.methods.dataframe,
      theta=quantiles,
      fdr=alpha,
      showplots=plotRejCurve,
      fromtool='edger',
      title=title
    )

    # save a copy of the independent filtering results
    contrastwise.output.list[[contrast.name]][["IF.metadata"]] <- edger.exactTest.IF.results
    
    
    if(!is.na(chosen_filter)){
      # get indices of genes whose raw read counts are greater than the chosen filter
      edgeR.exactTest.genes.passes_filter <- edger.exactTest.IF.results$filtering.methods[chosen_filter] > edger.exactTest.IF.results$relevant.threshold.per.method[chosen_filter]
      
      # run topTags only on genes whose raw read counts are greater than the chosen filter
      edgeR.exactTest.topTags <- topTags(
        edger.result$edgeR.exactTest[ edgeR.exactTest.genes.passes_filter, ], 
        n = sum(edgeR.exactTest.genes.passes_filter),
        adjust.method="BH"
      )
      
      colnames(edgeR.exactTest.topTags$table) <- paste(
        "edgeR.exactTest", 
        colnames(edgeR.exactTest.topTags$table), 
        sep = "_"
      )
      
      edgeR.exactTest.topTags=as.data.frame(
        edgeR.exactTest.topTags, 
        row.names=rownames(edgeR.exactTest.topTags$table)
      )
      
      contrastwise.output.list[[contrast.name]][["IF.results"]] <- edgeR.exactTest.topTags  
    }
    
  }
  return(contrastwise.output.list)
}

edgeR_DGE <- function(DGE_obj, 
                      design, 
                      coef=NULL, 
                      contrast=NULL,
                      useLRT=FALSE, 
                      useQLF=FALSE, 
                      useEXACT=TRUE) {
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

DESeq2_DGE_analysis <- function(DESeq2_dataset, 
                                alpha, 
                                contrasts, 
                                showplots = TRUE, 
                                filtering.methods = NA, 
                                quantiles = seq(from=0.0, to=0.99, by=0.025), 
                                chosen_filter = "mean", 
                                pAdjustMethod = "BH", 
                                lfcThreshold = 0,
                                verbose = FALSE) {
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
  
  #DESeq2_dataset <- estimateSizeFactors(DESeq2_dataset)
  DESeq2_dataset <- DESeq(DESeq2_dataset)
  DESeq2.out <- DESeq2_dataset

  means = rowMeans(counts(DESeq2_dataset, normalized=TRUE))
  medians = rowMedians(counts(DESeq2_dataset, normalized=TRUE))
  
  if(is.na(filtering.methods)) { 
    filtering.methods.dataframe <- data.frame(
      mean=means,
      #'min'= rowMin(counts(DESeq2_dataset, normalized=TRUE)),
      #'max'= rowMax(counts(DESeq2_dataset, normalized=TRUE)),
      median=medians
      #'secondlargest'= apply(counts(DESeq2_dataset, normalized=TRUE), 1,  function(row) sort(row, partial=length(row)-1)[length(row)-1])
    )
  } else {
    filtering.methods.dataframe <- filtering.methods
  }

  stopifnot(chosen_filter %in% colnames(filtering.methods.dataframe))

  contrastwise.output.list <- list()
  contrastwise.output.list[["DGE_obj"]] <- DESeq2.out
  
  # for every specified contrast, 
  # run IF on DESeq results
  for (contrast.name in colnames(contrasts)){
    contrast.conditions <- strsplit(contrast.name, "_")[[1]]
    
    if(verbose){
      print(contrast.name)
    }
    
    # get and save IF'd results using the chosen filter
    deseq2.result <- results(
      DESeq2.out,
      contrast = c(contrast.conditions[1], contrast.conditions[2], contrast.conditions[4]),
      #contrast = c("condition", contrast.conditions[1], contrast.conditions[2]),
      alpha = alpha,
      theta = quantiles, 
      pAdjustMethod = pAdjustMethod,
      lfcThreshold = lfcThreshold,
      independentFiltering = TRUE,
      filter = filtering.methods.dataframe[[chosen_filter]]
    )

    if(verbose) {
      print(summary(deseq2.result))
    }

    contrastwise.output.list[[contrast.name]] <- list(
      "DESeq2.results" = deseq2.result
    )

    title <- paste0(
      "DESeq2(", 
      contrast.name, 
      ") rejection curve (alpha (a.k.a. FDR) = ",
      alpha,
      ")\n# rejections (# DEGs) = ",
      sum(deseq2.result['padj'][[1]] < alpha, na.rm=TRUE),
      ", threshold = ",
      round(metadata(deseq2.result)$filterThreshold,2)
    )

    # plot the rejection curve for the chosen filter
    if(showplots) {      
      plot(
        metadata(deseq2.result)$filterNumRej, 
        type="b", 
        ylab="number of rejections",
        xlab="quantiles of filter", 
        main=title
      )
      lines(
        metadata(deseq2.result)$lo.fit, 
        col="red"
      )
      abline(v=metadata(deseq2.result)$filterTheta)
      print(resultsNames(DESeq2_dataset))
      print(contrast.name)
      #plotMA(
      #  lfcShrink(DESeq2_dataset, coef=contrast.name), 
      #  ylim=c(-2,2)
      #)
    }

    contrast.result.copy <- deseq2.result
    contrast.result.copy$unadjPvalues <- deseq2.result$pvalue
    IF.DESeq2.results <- independent_filtering(
      contrast.result.copy,
      filtering.methods.dataframe,
      theta = quantiles,
      title = title,
      fdr = alpha,
      showplots = TRUE, # still need to work on this functionality
      fromtool = "deseq"
    )
    #contrastwise.IF.DESeq2.results[[contrast.name]] <- IF.DESeq2.results
    # save a copy of the independent filtering results
    contrastwise.output.list[[contrast.name]][["DESeq2.IF.results"]] <- IF.DESeq2.results

    ##
    #if(!is.na(chosen_filter)){
    #  contrastwise.output.list[[contrast.name]][["IF.results"]] <- results(
    #    DESeq2.out, 
    #    contrast=c("condition", contrast.conditions[1], contrast.conditions[2]),
    #    alpha=alpha,
    #    pAdjustMethod=pAdjustMethod,
    #    filter=filtering.methods.dataframe[[chosen_filter]]
    #  )
    #}    
  }
  return(contrastwise.output.list)
}

design.pairs <- function(levels) {
  levels <- sort(levels)
  print(levels)
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
      colnames(design)[k] <- paste("condition","_",levels[i],"_vs_",levels[j],sep="")
    }
  return(design)
}

independent_filtering <- function(dge_obj.with.pvalues, 
                                  filtering.methods.dataframe = NULL,
                                  theta = seq(from=0.0, to=0.99, by=0.025), 
                                  fdr = 0.1, 
                                  showplots = TRUE, 
                                  p.adjust_method = "BH", 
                                  title = "Filtering Methods' Rejection Curves",
                                  fromtool = NA,
                                  genes_of_interest = c()) {
  #
  # Wrapper function to recreate DESeq2's independent filtering code.
  #
  # PARAMETERS:
  #
  #   dge_obj.with.pvalues: ... object along with a unadjPvalues field.
  #
  #   filtering.methods.dataframe:  dataframe containing the methods to assess in the independent filtering.
  #                                 defaults to:
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
  #   theta: (optional) sequence ranging from [0.0, 1.0) indicating which quantiles to consider dropping.
  #          defaults to seq(from=0.0, to=0.99, by=0.025).
  #
  #   fdr:  value for the fdr threshold.
  #         defaults to 0.1. 
  #   
  #   showplots:  boolean indicating whether to display plots or not.
  #               defaults to TRUE.
  #
  #   p.adjust_method:  (optional) string specifying the multiple-hypothesis 
  #                     p-value adjustment method.
  #                     gets passed onto filtered_R, see `genefilter` package.
  #                     defaults to "BH".
  #
  #   title: (optional) string used as title to the rejection curve plot.
  #          defaults to 'Filtering Methods' Rejection Curves'
  #
  #   fromtool: string among 'deseq', 'edger' indicating which tool generated `dge_obj.with.pvalues` .
  #
  # RETURNS:
  #   
  #   list of "rejections" and "best.numRej.per.Method"
  #
  # REFERENCE:
  #   https://bioconductor.org/packages/release/bioc/vignettes/genefilter/inst/doc/independent_filtering.pdf
    
  stopifnot(!is.na(fromtool) | !((tolower(fromtool) != "deseq") & (tolower(fromtool) != "edger")) )
  if(is.null(filtering.methods.dataframe)) {
    if(tolower(fromtool) == "deseq") {
      filtering.methods.dataframe <- data.frame(
        'mean'= rowMeans(counts(dge_obj.with.pvalues, normalized=TRUE)),
        'min'= rowMin(counts(dge_obj.with.pvalues, normalized=TRUE)),
        'max'= rowMax(counts(dge_obj.with.pvalues, normalized=TRUE)),
        'median'= rowMedians(counts(dge_obj.with.pvalues, normalized=TRUE)),
        'secondlargest'= apply(counts(dge_obj.with.pvalues, normalized=TRUE), 1,  function(row) sort(row, partial=length(row)-1)[length(row)-1])
      )
    }

    else { # are these normalized?
      filtering.methods.dataframe <- data.frame(
        'mean' = rowMeans(dge_obj.with.pvalues$counts),
        'min' = rowMin(dge_obj.with.pvalues$counts),
        'max' = rowMax(dge_obj.with.pvalues$counts),
        'median' = rowMedians(dge_obj.with.pvalues$counts),
        'secondlargest' = apply(dge_obj.with.pvalues$counts, 1,  function(row) sort(row, partial=length(row)-1)[length(row)-1])
      )
    }    
  }
  
  # check that all genes of interest are present
  if(!is.null(genes_of_interest)){  
    absent <- genes_of_interest[ !(genes_of_interest %in% rownames(filtering.methods.dataframe)) ]
    if(length(absent) > 0) {
      stop(
        paste0(
          "the following genes were not found in the row names:", 
          paste0(absent, sep='\n'),
          sep='\n'
        )
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

  genes_of_interest_survival <- list()
  if(!is.null(genes_of_interest)) {
    goi_rows <- which(rownames(filtering.methods.dataframe) %in% goi)
    for(i in 1:length(filtering.methods.dataframe)) {
      statistic_name <- colnames(filtering.methods.dataframe)[i]
      statistic_values <- filtering.methods.dataframe[,i]
      ecdf_fun <- ecdf(sort(statistic_values))
      goi_values <- statistic_values[goi_rows]
      goi_statistic_percentiles <- ecdf_fun(goi_values)
      goi_df <- cbind(goi_values, goi_statistic_percentiles)
      rownames(goi_df) <- genes_of_interest 
      colnames(goi_df) <- c(paste0(statistic_name, " statistic filter value"), "percentile value")
      print(goi_df)
      genes_of_interest_survival[[statistic_name]] <- goi_df
    }
  }
  
  best.numRej.per.method <- apply(rejections, 2, max, na.rm = TRUE)
  relevant.quantile.per.method <- list()
  relevant.threshold.per.method <- list()

  for(i in 1:length(filtering.methods.dataframe)) {
    q = theta[ which(rejections[,i] == max(rejections[,i])) ][1]
    # the extra [1] is to get the lowest quantile that yielded the max number of rejections
    # to filter out fewer genes
    filter_name <- colnames(filtering.methods.dataframe)[i]
    relevant.quantile.per.method[filter_name] = q
    relevant.column <- filtering.methods.dataframe[,filter_name]
    relevant.threshold.per.method[filter_name] = quantile(relevant.column,q)
  }

  print("relevant.quantile.per.method")
  print(relevant.quantile.per.method)

  if(showplots==TRUE) {
    colors=brewer.pal(ncol(filtering.methods.dataframe), "Set1")
    if(!is.null(genes_of_interest)) {
       par(mfrow=c(1+NCOL(filtering.methods.dataframe),1))
    }
    matplot(
      theta, 
      rejections, 
      type='l', 
      lty=1, 
      col=colors, 
      lwd=2, 
      xlab=paste0(expression(theta), " (filter's percentile threshold)"),
      ylab="Rejections (number of DEGs)", 
      main=title
    )
    
    for(i in 1:length(relevant.quantile.per.method)) {
      abline(
        v = relevant.quantile.per.method[colnames(relevant.quantile.per.method)[i]],
        col = colors[i]
      )
    }
    
    legend(
      "bottomleft",
      legend = make.legend(best.numRej.per.method, relevant.quantile.per.method, relevant.threshold.per.method),
      #legend=colnames(filtering.methods.dataframe), 
      fill=colors
    )
    #if(!is.null(genes_of_interest)) {
    #  for(r in 1:NCOL(filtering.methods.dataframe)+1) {
    #    plot(
    #      0,
    #      type="n", 
    #      xlim=c(0.0,1.1), 
    #      ylim=c(0,length(genes_of_interest)+1), 
    #      xlab="X Label", 
    #      ylab="Y Label", 
    #      main="Survival of Genes Of Interest"
    #    )
    #    statistic_name <- colnames(filtering.methods.dataframe)[r]
    #    survival_data <- genes_of_interest_survival[[statistic_name]]
    #    for(i in 1:length(genes_of_interest)){
    #      lines(
    #        0:survival_data[genes_of_interest[i],"percentile value"], 
    #        rep(i, a[i,2]), 
    #        type="l", col=i
    #      )
    #    }
    #  }
    #}
  }

  return(
    list(
    "rejections" = rejections, 
    "best.numRej.per.method" = best.numRej.per.method,
    "relevant.quantile.per.method" = relevant.quantile.per.method,
    "relevant.threshold.per.method" = relevant.threshold.per.method,
    "filtering.methods" = filtering.methods.dataframe,
    "genes_of_interest_survival" = genes_of_interest_survival
  )
  )
  
}

### end of dge_analysis.R ###


## start of visualization.R ##
library(dplyr)
library(tidyr)
library(reshape2)
library(factoextra)
library(RColorBrewer)
library(pheatmap)
library(regioneR)
library(chromoMap)
library(karyoploteR) 

plotLog1PReadCountsDistribution <- function(raw_reads_dataframe, 
                                            conditions = NA,
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


  if(is.na(conditions)) {
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
    readcounts.histogram <- ggplot(gather(as.data.frame(log1p(raw_reads_dataframe))), aes(value)) + 
      geom_histogram(binwidth = 1) + 
      facet_wrap(~key) +
      labs(x = "ln(1 + raw reads mapped to genes)")
    print(readcounts.histogram + ggtitle(title))
  }
  
  if(boxplot) {
    # is as.data.frame needed?
    readcounts.boxplot <- ggplot(
      data = reshape2::melt(as.data.frame(log1p(raw_reads_dataframe))), aes(x=variable, y=value)) + 
      geom_boxplot(aes(fill=variable)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs(y = "ln(1 + raw reads mapped to genes)", x = "sample")
    print(readcounts.boxplot + ggtitle(title))
  }
}

plotPheatMap <- function(data_obj, fromtool = NA, title = "", ...){
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
    sample.counts <- data_obj$counts
    sample.names <- colnames(data_obj)
  } else if (tolower(fromtool) == "deseq") {
    sample.counts <- counts(data_obj)
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

plotScree <- function(prcomp.obj, ...){
  # Wrapper around `fviz_eig` to plot a scree plot.
  #
  # PARAMETERS:
  #
  #   prcomp.obj: object obtained by prcomp.
  #
  # RETURNS
  #
  #   Nothing.
  print(
    fviz_eig(prcomp.obj, ...)
  )
}

plotPC1vsPC2 <- function( prcomp.obj, 
                          normalized.data,
                          repel=TRUE,
                          title="PC1 vs PC2\n(Large markers are centroids)",
                          ...){
  # Wrapper around `fviz_pca_ind` to plot a PC1 vs PC2 plot.
  #
  # PARAMETERS:
  #
  #   prcomp.obj: object obtained by prcomp.
  #
  #   normalized.data:  data.frame with a `condition` column to use when 
  #                     labelling points.
  #
  #   repel:  (optional) boolean indicating whether the points should be 'wiggled'
  #           when plotting to avoid overplotting.
  #           defaults to TRUE.
  #   
  #   title:  (optional) string to use as the plot's title.
  #           defaults to 'PC1 vs PC2\n(Large markers are centroids)'.
  #
  # RETURNS
  #
  #   Nothing.
  print(
    fviz_pca_ind(
      prcomp.obj, 
      #label=label, 
      habillage = normalized.data$condition,
      #geom = geom, 
      repel = repel,
      title = title,
      ...
    )
  )
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
        threshold_df[[names(named_num)[i]]]
      )
    )
  }
  return(out)
}
### end of visualization.R ###

