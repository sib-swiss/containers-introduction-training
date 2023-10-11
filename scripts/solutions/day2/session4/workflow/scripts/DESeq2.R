# Script designed to find Differentially Expressed Genes using DESeq2


# Redirects all R logs to snakemake log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


cat("Loading packages\n")
suppressMessages({
  library(apeglm)
  library(data.table)
  library(DESeq2)
  library(IHW)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  # Parallelization doesn't work well on Windows; if needed, comment the next
  # line and remove the parallel and BPPARAM parameters below
  library(BiocParallel)
})


cat("Getting data from snakemake\n")
read_counts <- unlist(snakemake@input$table)
condition <- "Treatment"
target <- "high"
control <- "low"
threshold <- 2
## To make sure threshold and a posteriori are valid numbers
if (threshold < 1) {
  cat("\tYou can't perform a DE analysis with a threshold < 1.")
  cat(" Please check <DE_contrats.tsv>. The analysis will now stop")
  quit(save = "no", status = 1)
} else if (threshold == 1) {
  aposteriori <- as.numeric(sub(pattern = "aposteriori", replacement = "",
                                x = unlist(snakemake@wildcards$aposteriori)))
  if (aposteriori <= 1) {
    cat("\tYou can't perform a DE analysis with an a posteriori threshold <= 1.")
    cat(" Please check <DE_contrats.tsv>. The analysis will now stop")
    quit(save = "no", status = 1)
  }
}
fit_type <- 'parametric'
filter <- T
alpha <- 0.05
threads <- snakemake@threads
DEG_table <- unlist(snakemake@output$deg)
plots_pdf <- unlist(snakemake@output$pdf)


# Logs variables imported from snakemake; mainly for debugging purpose
cat("Variables and variable types used for this analysis:\n")
cat("\tread_counts       <--  ", read_counts, ", variable type: ", data.class(read_counts),           "\n", sep = "")
cat("\tcondition         <--  ", condition, ", variable type: ", data.class(condition),               "\n", sep = "")
cat("\ttarget            <--  ", target, ", variable type: ", data.class(target),                     "\n", sep = "")
cat("\tcontrol           <--  ", control, ", variable type: ", data.class(control),                   "\n", sep = "")
cat("\tthreshold         <--  ", threshold, ", variable type: ", data.class(threshold),               "\n", sep = "")
cat("\tfit_type          <--  ", fit_type, ", variable type: ", data.class(fit_type),                 "\n", sep = "")
cat("\tfilter            <--  ", filter, ", variable type: ", data.class(filter),                     "\n", sep = "")
cat("\talpha             <--  ", alpha, ", variable type: ", data.class(alpha),                       "\n", sep = "")
cat("\tDEG_table         <--  ", DEG_table, ", variable type: ", data.class(DEG_table),               "\n", sep = "")
cat("\tplots_pdf         <--  ", plots_pdf, ", variable type: ", data.class(plots_pdf),               "\n", sep = "")


cat("Processing metadata\n")
cat("\tImporting metadata from\n")
metadata <- data.table(Sample = c("lowCO2_sample1", "lowCO2_sample2", "lowCO2_sample3",
                                         "highCO2_sample1", "highCO2_sample2", "highCO2_sample3"),
                              Condition_type = c("Treatment", "Treatment", "Treatment",
                                                 "Treatment", "Treatment", "Treatment"),
                              Condition_value = c("low", "low", "low",
                                                  "high", "high", "high"))

cat("\tSelecting metadata about <", condition, "> condition\n", sep = "")
metadata <- metadata[Condition_type == condition]
cat("\tRenaming <Condition_value> column to <", condition, ">\n", sep = "")
setnames(metadata, "Condition_value", condition)

cat("\tSelecting metadata about <", target, "> and <", control, ">\n", sep = "")
metadata <- metadata[c(target, control), on = condition]

cat("\tRemoving metadata empty rows\n")
if (length(which(rowSums(!is.na(metadata)) == 0)) != 0) {
  metadata <- metadata[- which(rowSums(!is.na(metadata)) == 0)]
}

cat("\tRemoving metadata empty columns\n")
if (length(which(colSums(!is.na(metadata)) == 0)) != 0) {
  metadata[, which(colSums(!is.na(metadata)) == 0) := NULL][, Assembly := NULL]
}

cat("\tReordering metadata rows\n")
metadata <- metadata[order(Sample)]  # Alphabetical order

cat("\tConverting metadata columns into factors\n")
metadata <- metadata[, lapply(.SD, as.factor)]


cat("Processing count data\n")
cat("\tImporting count data from <", read_counts, ">\n", sep = "")
cts <- fread(input = read_counts, header = T, sep = "\t", quote = F,
             stringsAsFactors = F, verbose = F, showProgress = T,
             data.table = T, nThread = threads)

cat("\tFiltering and reordering count data\n")
## Removes unused samples in count table and sorts its columns in the same
## order than sample columns of metadata (required by DESeq2)
cts <- cts[, .SD, .SDcols = c("Geneid", as.character(metadata[, Sample]))]

## Checks whether the reordering worked, message is mainly for debugging purpose
if (all(metadata[, Sample] == colnames(cts)[-1])) {
  cat("\tThe columns of the count table and the rows of the samples information")
  cat(" are identical. The DE analysis can continue\n")
} else {
  cat("\tThe columns of the count table are different from the rows of the")
  cat(" samples information table. This will cause DESeq2 to crash. Please")
  cat(" check the log and the sample names in count table. Analysis will now")
  cat(" stop\n")
  quit(save = "no", status = 1)
}

cat("\tRounding read counts\n")
for (col in 2:ncol(cts)){
  cts[, (col) := round(cts[[col]])]
}

cat("\tReplacing <Unknown/NA> read counts by <0>\n")
cts[is.na(cts)] <- 0  # Avoids a crash in DEseq2

cat("\tRemoving genes with less than <5> reads in at least 3 samples\n")
## Not mandatory, but it avoids to waste time on genes with low read counts
keep <- rowSums(cts >= 5) >= 3
counts_filtered <- cts[keep,]
cat("\t<", nrow(cts) - nrow(counts_filtered), "> genes removed\n", sep = "")

cat("\tSetting rownames in filtered table\n")
# Sets first column as rownames because DESeq2 needs them
counts_filtered <- as.data.frame(counts_filtered)
rownames(counts_filtered) <- counts_filtered$Geneid
counts_filtered$Geneid <- NULL


# Creates a DESeqDataSet
cat("No input of batch effect detected,")
cat(" creating DESeq dataset accordingly\n")
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata,
                              design = formula(paste("~", condition)),
                              tidy = F, ignoreRank = F)

cat("Setting <", control, "> as control condition for DE analysis\n", sep = "")
dds@colData[, condition] <- relevel(dds@colData[, condition], ref = control)

cat("Analysing data with DESeq2\n")
dds <- DESeq(object = dds, test = "Wald", sfType = "poscounts", betaPrior = F,
             fitType = fit_type, modelMatrixType = "standard", parallel = TRUE,
             minReplicatesForReplace = 3, BPPARAM = MulticoreParam(threads))


# Sample quality assessment (QA)
## Initializes a list where plots will be stored before printing
cat("Assessing samples quality:\n")
plot_list <- list()

## Perfoms QA
cat("\tApplying blind variance stabilising transformation\n")
blind_vst <- varianceStabilizingTransformation(dds, blind = T,
                                               fitType = fit_type)

# PCA plot
cat("\tExtracting PCA data\n")
pcaData <- plotPCA(blind_vst, intgroup = condition, returnData = T)
# Calculates the % of variance explained by PC1 and PC2
percentVar <- round(100 * attr(pcaData, "percentVar"))
cat("\t\tPlotting PCA\n")
tmp <- ggplot(data = pcaData, aes(x = PC1, y = PC2,
                                  color = .data[[condition]])) +
          geom_point(size = 3) +
          scale_x_continuous(breaks = pretty_breaks()) +
          scale_y_continuous(breaks = pretty_breaks()) +
          coord_fixed() +
          labs(title = bquote("PCA of"~italic(.("S. cerevisiae"))~"samples"),
             x = paste0("PC1: ", percentVar[1], "% variance"),
             y = paste0("PC2: ", percentVar[2], "% variance")) +
          theme(plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 18,
                                          margin = margin(t = 0, r = 0,
                                                          b = 10, l = 0)),
                plot.margin = margin(t = 10, r = 25, b = 5, l = 10),
                axis.title.y = element_text(margin = margin(t = 0, r = 10,
                                                            b = 0, l = 0)),
                axis.text = element_text(size = 11),
                axis.title = element_text(size = 13),
                panel.background = element_rect(color = "black",
                                                fill = "grey90"),
                panel.grid.major = element_line(color = "black",
                                                linewidth = 0.5,
                                                linetype = 3),
                panel.grid.minor = element_blank())
plot_list <- c(plot_list, list(tmp))
rm(tmp)

# Sample to sample distance heatmap
cat("\tExtracting distance data\n")
sample_dist <- dist(t(assay(blind_vst)))
sample_dist_mat <- as.matrix(sample_dist)
cat("\t\tPlotting distance heatmap\n")
rownames(sample_dist_mat) <- paste(blind_vst$Sample)
colnames(sample_dist_mat) <- paste(blind_vst$Sample)
suppressWarnings(
tmp <- pheatmap(mat = sample_dist_mat,
                color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                clustering_distance_rows = sample_dist,
                clustering_distance_cols = sample_dist,
                main = bquote("Distance between S. cerevisiae samples"),
                border_color = "grey50", treeheight_row = 75,
                treeheight_col = 75, legend = T, show_rownames = T,
                show_colnames = T, fontsize = 13, fontsize_row = 10,
                fontsize_col = 10, angle_col = 315, display_numbers = T,
                number_format = "%.1f", number_color = "red", na_col = "grey",
                fontsize_number = 10, silent = T)
)
plot_list <- c(plot_list, list(tmp))
rm(tmp)

# Count matrix heatmap
cat("\tExtracting count data\n")
# This commands gives the mean of normalised counts across samples for the
# 20 genes with the most reads (on average)
select <- order(rowMeans(counts(dds, normalized = T)), decreasing = T)[1:20]
cat("\t\tPlotting counts heatmap\n")
# Creates a table that will be used for heatmap annotations
annot_col <- as.data.frame(colData(dds)[, condition])
rownames(annot_col) <- colData(dds)[, "Sample"]
colnames(annot_col) <- condition
suppressWarnings(
tmp <- pheatmap(assay(blind_vst)[select, ], legend = T, show_rownames = F,
                cluster_rows = F, cluster_cols = F, border_color = "grey50",
                na_col = "grey", annotation_col = annot_col, silent = T,
                show_colnames = T, fontsize = 13, fontsize_row = 9,
                fontsize_col = 9, angle_col = 315,
                main = bquote("log2 counts of the top 20 genes in S. cerevisiae"))
)
plot_list <- c(plot_list, list(tmp))
rm(tmp)


# Per-gene dispersion estimates and fitted mean-dispersion relationship
# to check dispersion estimation
# ggplot2 version of DESeq2 code. See Love et al., 2014
cat("Assessing DESeq fit:\n")
cat("\tExtracting dispersion estimates\n")
# Extracts strictly positive baseMean
to_keep <- mcols(dds)$baseMean > 0
basemeans <- mcols(dds)$baseMean[to_keep]
# Extracts dispersion estimates
disp <- mcols(dds)$dispGeneEst[to_keep]
# Sets lower bound of y-axis
ymin = 10^floor(log10(min(disp[disp > 0], na.rm = TRUE)) - 0.1)

cat("\tCalculating fitting estimator\n")
est <- signif(median(abs(mcols(dds)$dispGeneEst - mcols(dds)$dispFit)),
              digits = 4)

cat("\tPlotting dispersion estimates\n")
tmp <- ggplot() + geom_point(aes(x = basemeans, y = pmax(disp, ymin),
                                 color = "black"),
                             size = 0.5,
                             shape = ifelse(disp < ymin, 6, 20)) +
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          # Draws a circle over the outliers
          (if (!is.null(dispersions(dds))) {
            geom_point(aes(x = basemeans, y = dispersions(dds)[to_keep],
                           color = "dodgerblue"),
                       size = ifelse(mcols(dds)$dispOutlier[to_keep],
                                     2 * 0.45, 0.45),
                       shape = ifelse(mcols(dds)$dispOutlier[to_keep],
                                      1, 20),
                       stroke = ifelse(mcols(dds)$dispOutlier[to_keep],
                                       2, 1))}) +
          # Draws the fitted model
          (if (!is.null(mcols(dds)$dispFit)) {
            geom_point(aes(x = basemeans, y = mcols(dds)$dispFit[to_keep],
                           color = "red"),
                       size = 0.5, shape = 20)
          }) +
          labs(title = bquote("Dispersion estimates in S. cerevisiae samples"),
               x = "Mean of normalized counts",
               y = "Dispersion") +
          theme(plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 18,
                                          margin = margin(t = 0, r = 0,
                                                          b = 10, l = 0)),
                plot.margin = margin(t = 10, r = 25, b = 5, l = 10),
                axis.title.y = element_text(margin = margin(t = 0, r = 10,
                                                            b = 0, l = 0)),
                axis.text = element_text(size = 11),
                axis.title = element_text(size = 13),
                panel.background = element_rect(color = "black",
                                                fill = "grey95"),
                panel.grid.major = element_line(color = "black",
                                                linewidth = 0.5,
                                                linetype = 3),
                panel.grid.minor = element_blank(), legend.position = "right") +
          scale_color_identity(guide = "legend", name = NULL,
                               labels = c("Gene-est", "Fitted", "Final"),
                               breaks = c("black", "red", "dodgerblue")) +
          annotate(geom = "label", label = paste("Fitting estimator:\n", est),
                   size = 4.25, color = "black", alpha = 1,
                   x = max(basemeans) * 0.3, y = min(pmax(disp, ymin)) * 10)
plot_list <- c(plot_list, list(tmp))
rm(tmp)


# Cook's distances to check if a sample has a more outliers than the others
cat("Extracting Cook's distances\n")
cooks <- stack(as.data.frame(log10(assays(dds)[["cooks"]])))

cat("\tPlotting Cook's distances\n")
tmp <- ggplot(data = cooks, aes(x = ind, y = values)) +
          # Draws the dashed whiskers
          geom_boxplot(linetype = "dashed", coef = 20) +
          # Draws the grey boxes
          geom_boxplot(fill = "grey", coef = 0, outlier.shape = NA) +
          # Draws the hinges
          stat_boxplot(geom = "errorbar", aes(ymin = after_stat(ymax)),
                       coef = 20, width = 0.3) +
          stat_boxplot(geom = "errorbar", aes(ymax = after_stat(ymin)),
                       coef = 20, width = 0.3) +
          scale_y_continuous(breaks = pretty_breaks()) +
          coord_flip() +
          labs(title = bquote("Cook's distances in S. cerevisiae samples"),
               y = "log10(Cook's distance)") +
          theme(plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 18,
                                          margin = margin(t = 0, r = 0,
                                                          b = 10, l = 0)),
                plot.margin = margin(t = 10, r = 25, b = 5, l = 10),
                axis.title.x = element_text(margin = margin(t = 5, r = 0,
                                                            b = 0, l = 0)),
                axis.title.y = element_blank(),
                axis.text = element_text(size = 11),
                axis.title = element_text(size = 13),
                axis.text.x = element_text(margin = margin(t = 5, r = 0,
                                                           b = 0, l = 0)),
                axis.text.y = element_text(margin = margin(t = 0, r = 5,
                                                           b = 0, l = 5)),
                panel.background = element_rect(color = "black",
                                                fill = "grey97"),
                panel.grid.major.x = element_line(color = "grey85",
                                                  linewidth = 0.5,
                                                  linetype = 1),
                panel.grid.major.y = element_line(color = "grey75",
                                                  linewidth = 0.5,
                                                  linetype = 3),
                panel.grid.minor = element_blank())
plot_list <- c(plot_list, list(tmp))
rm(tmp)


# Creates a contrast variable in a DESeq2 style;
# displaying it in log is mainly for debugging purpose
contrast <- paste(condition, target, "vs", control, sep = "_")
cat("Current DESeq2 contrast is <", contrast, ">\n", sep = "")

cat("Extracting DE results for this contrast\n")
## Tests for a log fold change of log2(threshold); if threshold is
## negative/null/missing, log2(threshold) is set to 0 (= no statistical
## test on the LFC)
res <- results(object = dds, name = contrast, lfcThreshold = log2(threshold),
               altHypothesis = "greaterAbs", alpha = alpha, filterFun = ihw,
               independentFiltering = filter, pAdjustMethod = "BH",
               cooksCutoff = T, parallel = TRUE,
               BPPARAM = MulticoreParam(threads))
cat("\nResults summary:")
summary(res)

cat("Shrinking log2 fold changes\n")
res_LFC <- lfcShrink(dds = dds, res = res, coef = contrast, type = "apeglm",
                     svalue = F, parallel = TRUE,
                     BPPARAM = MulticoreParam(threads))

cat("\tAdding shrunken lfc and lfcSE to result table\n")
new_cols <- cbind(res_LFC$log2FoldChange, res_LFC$lfcSE)
res[, c("shrunken_log2FoldChange", "shrunken_lfcSE")] <- new_cols

cat("\tReordering table\n")
res <- res[, c("baseMean", "stat", "weight",
               "log2FoldChange", "lfcSE",
               "shrunken_log2FoldChange", "shrunken_lfcSE",
               "pvalue", "padj")]
names(res)[2] <- "Wald_stat"


# MA plot of log2 fold-changes and shrunken log2 fold-changes
cat("Extracting log ratio and mean averages\n")
if (threshold == 1) {
  isDE <- ifelse(is.na(res[["padj"]]), FALSE,
                 res[["padj"]] < alpha &
                   abs(res[["shrunken_log2FoldChange"]]) > log2(aposteriori))
  res_MA <- data.frame(mean = res[["baseMean"]],
                       lfc = res[["log2FoldChange"]],
                       isDE = isDE)
  res_MA_shrunken <- data.frame(mean = res[["baseMean"]],
                                lfc = res[["shrunken_log2FoldChange"]],
                                isDE = isDE)
} else {
  res_MA <- DESeq2::plotMA(res, alpha = alpha, returnData = T)
  res_MA_shrunken <- DESeq2::plotMA(res_LFC, alpha = alpha, returnData = T)
}

cat("\tPlotting MA values\n")
tmp <- ggplot(data = res_MA, aes(x = mean, y = lfc, color = isDE)) +
          geom_point(size = 1) +
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          scale_y_continuous(limits = c(-3, 3), oob = squish) +
          geom_hline(yintercept = 0, color = "grey40",
                     linewidth = 0.75, linetype = 2) +
          labs(title = bquote("MA plot -"~.(target)~"vs"~.(control)~
                              "in S. cerevisiae samples"),
               x = "Mean of normalised counts",
               y = bquote("log"[2]~"fold-change")) +
          theme(plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 18,
                                          margin = margin(t = 0, r = 0,
                                                          b = 10, l = 0)),
                plot.margin = margin(t = 10, r = 25, b = 5, l = 10),
                axis.title.x = element_text(margin = margin(t = 10, r = 0,
                                                            b = 0, l = 0)),
                axis.title.y = element_text(margin = margin(t = 0, r = 10,
                                                            b = 0, l = 0)),
                axis.text = element_text(size = 13),
                axis.title = element_text(size = 14),
                panel.background = element_rect(color = "black",
                                                fill = "grey97"),
                panel.grid.major = element_line(color = "black",
                                                linewidth = 0.5,
                                                linetype = 3),
                panel.grid.minor = element_blank()) +
          scale_color_manual(name = "DE status", labels = c("No", "Yes"),
                             values = c("grey70", "blue"))
plot_list <- c(plot_list, list(tmp))
rm(tmp)

cat("\tPlotting shrunken MA values\n")
tmp <- ggplot(data = res_MA_shrunken, aes(x = mean, y = lfc, color = isDE)) +
          geom_point(size = 1) +
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
          scale_y_continuous(limits = c(-3, 3), oob = squish) +
          geom_hline(yintercept = 0, color = "grey40",
                     linewidth = 0.75, linetype = 2) +
          labs(title = bquote("MA plot -"~.(target)~"vs"~.(control)~"in"~
                              italic(.("S. cerevisiae"))~"samples"),
               x = "Mean of normalised counts",
               y = bquote("shrunken log"[2]~"fold-change")) +
          theme(plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 18,
                                          margin = margin(t = 0, r = 0,
                                                          b = 10, l = 0)),
                plot.margin = margin(t = 10, r = 25, b = 5, l = 10),
                axis.title.x = element_text(margin = margin(t = 10, r = 0,
                                                            b = 0, l = 0)),
                axis.title.y = element_text(margin = margin(t = 0, r = 10,
                                                            b = 0, l = 0)),
                axis.text = element_text(size = 13),
                axis.title = element_text(size = 14),
                panel.background = element_rect(color = "black",
                                                fill = "grey97"),
                panel.grid.major = element_line(color = "black",
                                                linewidth = 0.5,
                                                linetype = 3),
                panel.grid.minor = element_blank()) +
          scale_color_manual(name = "DE status", labels = c("No", "Yes"),
                             values = c("grey70", "blue"))
plot_list <- c(plot_list, list(tmp))
rm(tmp)


# Volcano plot of shrunken LFC and adjusted pvalues
cat("Preparing data for volcano plot\n")
# Creates a copy of the data
volcano <- suppressWarnings(as.data.table(res[order(res$padj), ],
                                          keep.rownames = "Geneid"))
# Adds a column about differential expression and sorts genes
volcano[, isDE := "Not DE"]
volcano[shrunken_log2FoldChange >= log2(threshold) & padj < alpha,
        isDE := "Up-regulated"]
volcano[shrunken_log2FoldChange <= -log2(threshold) & padj < alpha,
        isDE := "Down-regulated"]
# Adds a column for labels for DEG
volcano[volcano[, isDE != "Not DE"], label := volcano[isDE != "Not DE",
                                                      Geneid]]
cat("\tCreating volcano plot\n")
tmp <- ggplot(data = volcano, aes(x = shrunken_log2FoldChange,
                                  y = -log10(padj), shape = isDE)) +
          # Sets color for point outline and point filling, depending of padj
          geom_point(aes(fill = padj, color = padj)) +
          scale_colour_gradient(low = "red", high = "blue",
                                guide = "colorbar",
                                name = "Adjusted p-value") +
          scale_fill_gradient(low = "red", high = "blue",
                              guide = "colorbar",
                              name = "Adjusted p-value") +
          scale_shape_manual(name = "DE status", values = c(24, 25, 21),
                             drop = F,
                             # Sets point shape, depending of DE status
                             breaks = c("Up-regulated",
                                        "Down-regulated",
                                        "Not DE")) +
          # Adds labels
          geom_text_repel(aes(label = label), size = 3, na.rm = T) +
          # Adds an horizontal line for the p-value threshold and a legend
          geom_hline(yintercept = -log10(alpha), col = "red", linetype = 2) +
          geom_text(x = -Inf, y = -log10(alpha), label = "p = 0.05",
                    color = "red", hjust = 1.1) +
          # Adds vertical lines for the log2FoldChange thresholds if
          # it's different from 0
          geom_vline(xintercept = c(-log2(threshold), log2(threshold)),
                     col = "red", linetype = 2) +
          geom_text(x = -log2(threshold), y = -Inf,
                    label = -round(log2(threshold), 2),
                    color = "red", vjust = 1.5) +
          geom_text(x = log2(threshold), y = -Inf,
                    label = round(log2(threshold), 2),
                    color = "red", vjust = 1.5) +
          # Rest is graphical improvement
          scale_x_continuous(breaks = pretty_breaks()) +
          scale_y_continuous(breaks = pretty_breaks()) +
          labs(title = bquote("Volcano plot of"~.(target)~"vs"~.(control)~
                              "in"~italic(.("S. cerevisiae"))~"samples"),
               x = bquote("Log"[2]~"fold change (shrunken)"),
               y = bquote("-log"[10]~"(padj)")) +
          theme(plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 18,
                                          margin = margin(t = 0, r = 0,
                                                          b = 10, l = 0)),
                plot.margin = margin(t = 10, r = 25, b = 5, l = 10),
                axis.title.x = element_text(margin = margin(t = 10, r = 0,
                                                            b = 0, l = 0)),
                axis.title.y = element_text(margin = margin(t = 0, r = 10,
                                                            b = 0, l = 0)),
                axis.text = element_text(size = 13),
                axis.title = element_text(size = 14),
                legend.title = element_text(size = 13),
                panel.background = element_rect(color = "black",
                                                fill = "grey97"),
                panel.grid.major = element_line(color = "grey75",
                                                linewidth = 0.5,
                                                linetype = 1),
                panel.grid.minor = element_blank()) +
          guides(shape = guide_legend(order = 1,
                                      override.aes = list(fill = "black",
                                                          size = 2.5))) +
          coord_cartesian(clip = "off")
plot_list <- c(plot_list, list(tmp))
rm(tmp)


cat("Sorting results table by smallest adjusted pvalues\n")
final <- suppressWarnings(as.data.table(res[order(res$padj), ],
                                        keep.rownames = "Geneid"))


cat("Finding missing genes if there are any\n")
missing <- !(cts$Geneid %in% rownames(counts_filtered))

# Adds missing genes to the final table
# (removed before because they're lowly expressed)
if (any(missing)) {
  cat("\tAdding missing genes to result table\n")
  unexpressed_genes <- as.data.table(matrix(data = NA, ncol = ncol(final),
                                            nrow = nrow(cts[missing, "Geneid"])))
  colnames(unexpressed_genes) <- colnames(final)
  unexpressed_genes[, Geneid := cts[missing, "Geneid"]]
  unexpressed_genes[, baseMean := rep(x = 0, nrow(unexpressed_genes))]
  final <- rbind(final, unexpressed_genes)
} else {

  cat("\tNo missing genes found\n")

}


cat("Saving result tables:\n")
deg_final <- final[abs(shrunken_log2FoldChange) >= ifelse(threshold > 1,
                                                          log2(threshold),
                                                          log2(aposteriori))
                   & padj < alpha, ]
fwrite(x = deg_final, file = DEG_table, quote = F, sep = "\t",
       col.names = T, na = NA)
cat("\tDEG table saved in <", DEG_table, ">\n", sep = "")


# Prints all previous plots in a pdf
cat("Saving plots\n")
pdf(file = plots_pdf, width = 9, height = 7, onefile = T, title = plots_pdf)
for (i in 1:length(plot_list)){
  # pheatmap objects sometimes stack on other objects in pdf, this avoids it
  if (data.class(plot_list[[i]]) == "pheatmap") {
    grid::grid.newpage()
    print(plot_list[[i]])
  }else{
    print(plot_list[[i]])
  }
}
invisible(dev.off())
cat("All plots saved in <", plots_pdf, ">", sep = "")
