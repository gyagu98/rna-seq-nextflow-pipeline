#!/usr/bin/env Rscript
# Two-group, unpaired gene-level differential expression from Salmon.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag) {
  i <- which(args == flag)
  if (length(i) != 1L || i + 1L > length(args)) stop("Required argument: ", flag)
  args[i + 1L]
}
samples <- read.csv(arg("--samplesheet"), stringsAsFactors = FALSE)
required <- c("sample", "condition")
if (!all(required %in% names(samples)) || nrow(samples) == 0L)
  stop("Samplesheet requires sample and condition columns and at least one row")
if (anyNA(samples[, required]) || any(!nzchar(samples$sample)) ||
    any(!nzchar(samples$condition)) || anyDuplicated(samples$sample))
  stop("Sample IDs must be unique; sample and condition must be nonempty")
if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", samples$sample)))
  stop("Sample IDs may contain letters, digits, underscores, dots and hyphens")
reference <- arg("--reference")
comparison <- arg("--comparison")
if (identical(reference, comparison) ||
    !setequal(unique(samples$condition), c(reference, comparison)))
  stop("Exactly two conditions matching reference and comparison are required")
if (any(table(samples$condition) < 2L))
  stop("At least two independent biological replicates per condition are required")
files <- file.path(arg("--quants"), samples$sample, "quant.sf")
if (!all(file.exists(files))) stop("Missing quant.sf: ", paste(files[!file.exists(files)], collapse=", "))
names(files) <- samples$sample
mapping <- read.delim(arg("--tx2gene"), stringsAsFactors=FALSE)
if (!identical(names(mapping), c("transcript", "gene")) || nrow(mapping) == 0L ||
    anyNA(mapping) || any(!nzchar(mapping$transcript)) || any(!nzchar(mapping$gene)) ||
    anyDuplicated(mapping$transcript))
  stop("tx2gene must have unique transcript IDs and nonempty transcript/gene columns")
for (f in files) {
  q <- read.delim(f, stringsAsFactors=FALSE)
  if (!all(c("Name", "Length", "EffectiveLength", "TPM", "NumReads") %in% names(q)) ||
      nrow(q) == 0L || anyDuplicated(q$Name))
    stop("Invalid Salmon table: ", f)
  if (!all(q$Name %in% mapping$transcript))
    stop("Every quantified transcript must have an exact tx2gene mapping: ", f)
}
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})
samples$condition <- factor(samples$condition, levels=c(reference, comparison))
rownames(samples) <- samples$sample
txi <- tximport(files, type="salmon", tx2gene=mapping, dropInfReps=TRUE)
dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)
keep <- rowSums(counts(dds)) >= 10
if (sum(keep) < 2L) stop("Fewer than two genes pass the total count >=10 filter")
dds <- DESeq(dds[keep, ])
res <- results(dds, contrast=c("condition", comparison, reference), alpha=0.05)
out <- arg("--out")
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
table <- data.frame(gene=rownames(res), as.data.frame(res), row.names=NULL)
table <- table[order(table$padj, na.last=TRUE), ]
write.csv(table, out, row.names=FALSE)
write.csv(data.frame(gene=rownames(dds), counts(dds, normalized=TRUE)),
          file.path(dirname(out), "normalized_counts.csv"), row.names=FALSE)
pdf(file.path(dirname(out), "MA_plot.pdf"))
plotMA(res, alpha=0.05)
dev.off()
transformed <- varianceStabilizingTransformation(dds, blind=TRUE)
pdf(file.path(dirname(out), "PCA.pdf"))
print(plotPCA(transformed, intgroup="condition"))
dev.off()
writeLines(c("Design: ~ condition (unpaired; no batch adjustment)",
             paste("Contrast:", comparison, "versus", reference),
             "Filter: total estimated count >=10; BH adjusted p-values; alpha=0.05",
             "Log2 fold changes are unshrunken. NA adjusted p-values must not be treated as significant.",
             "Replicate independence and confounding must be assessed by the investigator.",
             capture.output(sessionInfo())),
           file.path(dirname(out), "analysis_session.txt"))
