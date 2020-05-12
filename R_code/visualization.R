library(dplyr)
library(tidyr)
library(reshape2)
library(factoextra)
library(RColorBrewer)
library(pheatmap)

plotLog1PReadCountsDistribution <- function(raw_reads_dataframe, 
                                        histogram=TRUE, boxplot=TRUE,
                                        verbose=FALSE){
    #
    conditions <- sub(
        "\\_.*", "",
        colnames(raw_reads_dataframe)
    )
    if(verbose){
        print("Inferred following conditions from filtered count matrix:\n")
        print(conditions)
    }
    colors = c(1:length(unique(conditions)))
    
    if(histogram) {
        readcounts.histogram <- ggplot(gather(as.data.frame(log1p(raw_reads_dataframe))), aes(value)) + 
            geom_histogram(binwidth = 1) + 
            facet_wrap(~key) +
            labs(x = "ln(1 + raw reads mapped to genes)")
        print(readcounts.histogram)
    }

    if(boxplot) {
        readcounts.boxplot <- ggplot(
            data = reshape2::melt(as.data.frame(log1p(raw_reads_dataframe))), aes( x=variable, y=value)) + 
            geom_boxplot(aes(fill=variable)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            labs(y = "ln(1 + raw reads mapped to genes)", x = "sample")
        print(readcounts.boxplot)
    }
}

plotPheatMap <- function(DESeq_normalized_data){
    sample.dists <- as.matrix(dist(t(assay(DESeq_normalized_data))))
    # round to integer
    sample.dists.vals <- round(sample.dists)
    # remove NAs
    sample.dists.vals[ is.na(sample.dists.vals) ] <- ""
    rownames(sample.dists) <- colnames(DESeq_normalized_data)
    colnames(sample.dists) <- colnames(DESeq_normalized_data)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pm <- pheatmap(sample.dists,
        clustering_distance_rows=sample.dists,
        clustering_distance_cols=sample.dists,
        col=colors,
        main="PheatMap",
        display_numbers = sample.dists.vals
    )
    print(pm)
}

plotScree <- function(prcomp.obj){
    p <- fviz_eig(prcomp.obj)
    print(p)
}

plotPC1vsPC2 <- function(prcomp.obj, normalized.data, 
                        label="none", geom=c("point", "text"),
                        repel=TRUE){
    p <- fviz_pca_ind(
        prcomp.obj, label=label, 
        habillage=normalized.data$condition,
        geom=geom, repel=repel
    )
    print(p)
}

plotRejectionCurve <- function(filtering.methods, rejections, theta){
    # 
    colors = brewer.pal(ncol(filtering.methods), "Set1")
    matplot(theta, rejections, type='l', lty=1, col=colors, 
            lwd=2, xlab=expression(theta), ylab="Rejections", 
            main=title)
    legend("bottomleft", legend=colnames(filtering.methods), fill=colors)
}