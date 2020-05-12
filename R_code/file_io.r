# external dependencies
library(biomaRt)
#library(here)

# local dependencies
#print(getwd())
#source(here::here('visualization.R'))

rawCountsMatrix_to_dataframe <- function( rawCountsMatrix_filepath,
                                          sep=",", header=TRUE, row.names=1,
                                          keep_columns=c(),drop_columns=c(),
                                          make_ensembl_to_symbol=TRUE,
                                          make_histogram=FALSE, make_boxplot=FALSE,
                                          verbose=FALSE) {
  
  # Do-it-all function that 
  # 1. loads the contents of the `sep`-separated file specified by 
  # `rawCountsMatrix_filepath` (passing the `sep`,`header`, and `row.names`
  # parameter values to the read.table function), 
  # 2. keeps the columns specified by `keep_columns` (also drops columns 
  # specified by `drop_columns`),
  # 3. creates a vector of condition labels inferred from the retained columns'
  # names, and
  # 4. returns the filtered dataframe, the vector of condition labels
  
  # PARAMETERS:
  #
  #   rawCountsMatrix_filepath: string path to the file containing raw counts.
  #
  #   sep, header, row.names:   optional parameters to pass onto read.table.
  #
  #                             default to sep=",", header=TRUE, row.names=1.
  #
  #   keep_columns, drop_columns:   vectors of column names (strings) to keep/drop.
  #                                 cannot contain any shared strings.
  #  
  #   `make_ensembl_to_symbol`: if TRUE, will create a [ENSEMBL ID, Gene Symbol] 
  #                             dataframe and return it as the third item in the returned
  #                             list. Assumes that the ENSEMBL IDs are the row names of the
  #                             `raw.data` dataframe, and that the dataframe has a "Symbol" column.
  #                           
  #                             if FALSE, the third item in the returned list will be NULL.
  #                             
  #                             default = TRUE.
  #
  #   `makeplots`:  indicates whether to produce exploratory plots or not.
  #
  #                 default = TRUE.
  #
  #   `verbose`:  indicator of verbosity. 
  #
  #               default = TRUE.
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
  # TODO:
  #   nicer graphs with ggplot
  
  if(length(intersect(keep_columns, drop_columns)) > 0) {
    stop("`keep_columns` and `drop_columns` aren't mutually exclusive")
  }
  
  if( (length(keep_columns) == 0) & (length(drop_columns) == 0) ) {
    stop("at least one of `keep_columns` and `drop_columns` must be non-empty")
  }
  
  if(verbose) {
    print("Loading raw counts matrix")
  }
  
  # load data into dataframe
  raw.data <- data.frame(
    read.csv(
      rawCountsMatrix_filepath, 
      sep=sep, 
      header=header, 
      row.names=row.names,
      encoding='utf-8'
    )
  )
  
  print(head(raw.data))
  print(tail(raw.data))
  
  if(verbose){
    print("Dimensions of raw count matrix:")
    print(dim(raw.data))
    print("Columns of raw count matrix:")
    print(colnames(raw.data))
    print("Head of raw count matrix:")
    print(head(raw.data))
  }
  
  if(make_ensembl_to_symbol) {
    if(verbose) print("Creating ENSEMBL ID <-> Gene Symbol dataframe")
    # create [ENSEMBLE ID, Gene Symbol] dataframe for future lookups
    ENSEMBL_ID_to_Symbol <- data.frame(
      rownames(raw.data), raw.data$Symbol
    )
    colnames(ENSEMBL_ID_to_Symbol) <- c("ENSEMBLID", "Symbol")
    
  } else {
    make_ensembl_to_symbol=FALSE
  }
  
  if(verbose) {
    print("Keeping the following columns in raw counts matrix:")
    print(keep_columns)
  }
  # positive column selection
  raw.data <- raw.data[,which(names(raw.data) %in% keep_columns)]
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
        histogram=make_histogram, boxplot=make_boxplot,
        verbose=verbose
    )
  }

  if(make_ensembl_to_symbol) {
    list_to_return <- list("raw.data" = raw.data, "conditions" = conditions, "ENSEMBL_ID_to_Symbol" = ENSEMBL_ID_to_Symbol)  
  } else {
    list_to_return <- list("raw.data" = raw.data, "conditions" = conditions)
  }
  
  return(list_to_return)
}

write.gct <- function(gct, filename, check.file.extension=TRUE){
  #
  # save a GCT result to a file, ensuring the filename has the extension .gct
  #
  # shamelessly ripped off of:
  # https://github.com/genepattern/DESeq2/blob/abac851160ed1f84af2aa06f85b497c4b144656e/src/common.R
  # in order to produce the DESeq2-normalized counts file (.gct format)
  # to use as input to GSEA
  
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

write.factor.to.cls <- function(factor, filename) {
  # writes a factor to a cls file
  #
  # shamelessly ripped off of:
  # https://github.com/genepattern/DESeq2/blob/abac851160ed1f84af2aa06f85b497c4b144656e/src/common.R
  # in order to generate the .cls file needed for GSEA
  # modifications: removed check.file.extension functionality.
  
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

get.conversion.map <- function( raw.counts.data, biomart="ensembl", dataset="hsapiens_gene_ensembl",
                                BMattributes=c('ensembl_gene_id','external_gene_name'),
                                BMfilters='ensembl_gene_id',
                                str_replace_pattern=".[0-9]+$",
                                str_replace_replacement=""){
  mart <- useMart(
    biomart=biomart,
    dataset=dataset
  )
  x.to.y.map <- getBM(
    attributes=BMattributes,
    filters=BMfilters,
    values=str_replace(
      rownames(raw.counts.data$raw.data),
      pattern=str_replace_pattern,
      replacement=str_replace_replacement
    ),
    mart=mart
  )
  return(x.to.y.map)
}