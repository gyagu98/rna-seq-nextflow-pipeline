# Standalone integration test: real DESeq2/tximport on synthetic Salmon tables.
script <- normalizePath("bin/deseq2_analysis.R")
root <- tempfile("deseq-test-")
dir.create(root)
set.seed(21)
samples <- data.frame(sample=paste0("s", 1:8),
                      condition=rep(c("control", "treated"), each=4))
write.csv(samples, file.path(root, "samples.csv"), row.names=FALSE)
mapping <- data.frame(transcript=paste0("tx", 1:500), gene=paste0("gene", 1:500))
write.table(mapping, file.path(root, "map.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
for (i in 1:8) {
  dir.create(file.path(root, samples$sample[i]))
  mu <- rep(100, 500)
  if (i > 4) mu[1:30] <- 1600
  n <- rnbinom(500, mu=mu, size=20) + 1
  q <- data.frame(Name=mapping$transcript, Length=1000, EffectiveLength=800,
                  TPM=n/sum(n)*1e6, NumReads=n)
  write.table(q, file.path(root, samples$sample[i], "quant.sf"),
              sep="\t", quote=FALSE, row.names=FALSE)
}
run <- function(sheet, out) {
  system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(script, "--samplesheet", sheet, "--quants", root,
       "--tx2gene", file.path(root, "map.tsv"), "--reference", "control",
       "--comparison", "treated", "--out", out)))
}
out <- file.path(root, "positive", "results.csv")
stopifnot(run(file.path(root, "samples.csv"), out) == 0)
res <- read.csv(out)
stopifnot(all(c("gene", "log2FoldChange", "pvalue", "padj") %in% names(res)),
          median(res$log2FoldChange[res$gene %in% mapping$gene[1:30]]) > 2,
          any(res$padj < 0.05, na.rm=TRUE),
          all(file.exists(file.path(dirname(out),
            c("normalized_counts.csv", "PCA.pdf", "MA_plot.pdf", "analysis_session.txt")))))
# No replication must fail and must not create fabricated results.
write.csv(samples[c(1,5), ], file.path(root, "unreplicated.csv"), row.names=FALSE)
bad <- file.path(root, "unreplicated", "results.csv")
stopifnot(run(file.path(root, "unreplicated.csv"), bad) != 0, !file.exists(bad))
unlink(file.path(root, "s1", "quant.sf"))
missing <- file.path(root, "missing", "results.csv")
stopifnot(run(file.path(root, "samples.csv"), missing) != 0, !file.exists(missing))
unlink(root, recursive=TRUE)
cat("Synthetic integration and failure-path tests passed.\n")
