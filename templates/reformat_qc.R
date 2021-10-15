#!/usr/bin/env Rscript

library(monocle3)

dir.create("temp_fold")
cds <- readRDS("$cds_object")
cell_qc <- read.csv("$cell_qc")

if(nrow(pData(cds)) > 0) {
    if("$params.skip_doublet_detect" == 'false') {
        scrublet_out <- read.csv("$scrub_csv", header=F)
        pData(cds)\$scrublet_score <- scrublet_out\$V1
        pData(cds)\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
        cell_qc\$scrublet_score <- scrublet_out\$V1
        cell_qc\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
    }
}

write.csv(cell_qc, quote=FALSE, file="temp_fold/$cell_qc")

dup_stats <- read.table(paste0("$key", ".duplication_rate_stats.txt"))

df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3,
                    duplication_rate = dup_stats\$V4,
                    doublet_count = sum(cell_qc\$scrublet_call == "Doublet", na.rm=TRUE),
                    doublet_perc = paste0(round(sum(cell_qc\$scrublet_call == "Doublet",
                                                    na.rm=TRUE)/nrow(cell_qc) * 100, 1), "%"),
                    doublet_NAs=sum(is.na(cell_qc\$scrublet_call)))

write.csv(df, file=paste0("$key", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)
saveRDS(cds, file="temp_fold/$cds_object")

if ("$key" == "Barnyard") {
    fData(cds)\$mouse <- grepl("ENSMUSG", fData(cds)\$id)
    fData(cds)\$human <- grepl("ENSG", fData(cds)\$id)

    pData(cds)\$mouse_reads <- Matrix::colSums(exprs(cds)[fData(cds)\$mouse,])
    pData(cds)\$human_reads <- Matrix::colSums(exprs(cds)[fData(cds)\$human,])
    pData(cds)\$total_reads <- pData(cds)\$mouse_reads + pData(cds)\$human_reads
    pData(cds)\$human_perc <- pData(cds)\$human_reads/pData(cds)\$total_reads
    pData(cds)\$mouse_perc <- pData(cds)\$mouse_reads/pData(cds)\$total_reads
    pData(cds)\$collision <- ifelse(pData(cds)\$human_perc >= .9 | pData(cds)\$mouse_perc >= .9, FALSE, TRUE)

    collision_rate <- round(sum(pData(cds)\$collision/nrow(pData(cds))) * 200, 1)
    fileConn<-file("Barn_collision.txt")
    writeLines(paste0("$key", "\t", collision_rate, "%"), fileConn)
    close(fileConn)
} else {
    fileConn<-file("${key}_no_collision.txt")
    writeLines(paste0("$key", "\t", "NA"), fileConn)
    close(fileConn)
}