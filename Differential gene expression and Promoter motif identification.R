library(tidyverse)
library(devtools)
library(goseq)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

gff = import.gff("Brapa_reference/Brapa_gene_v1.5.gff")
gff$gene_id = ifelse(is.na(gff$ID),gff$Parent,gff$ID)
export(gff,"Brapa_reference/Brapa_gene_v1.5.gtf",format="gtf")
library(Rsubread)
readCounts = featureCounts(
  files=c("tophat_out-IMB211_All_A01_INTERNODE.fq/accepted_hits_A01.bam",
          "tophat_out-R500_All_A01_INTERNODE.fq/accepted_hits_A01.bam"),
  annot.ext="Brapa_reference/Brapa_gene_v1.5.gtf", 
  isGTFAnnotationFile=TRUE,
  GTF.featureType="CDS",
  GTF.attrType="gene_id"
)

counts.data = read_tsv("gh_internode_counts2.tsv")
counts.data = counts.data %>% filter(gene_id!="*")
counts.data[is.na(counts.data)] = 0
colnames(counts.data) = colnames(counts.data) %>% str_remove(., fixed(".1_matched.merged.fq.bam"))
counts.data = counts.data[rowSums(counts.data[,-1] > 10) >= 3,]
sample.description = tibble("sample" = colnames(counts.data)[-1])
sample.description$gt = str_replace(sample.description$sample, "(.*)_.*_.*_.*", "\\1")
sample.description$trt = str_replace(sample.description$sample, ".*_(.*)_.*_.*", "\\1")
sample.description$group = str_c(sample.description$gt, "_", sample.description$trt)
sample.description = sample.description %>%
  mutate(gt=factor(gt), 
         trt=factor(trt,levels = c("NDP","DP")))

# MDS with BCV
library(edgeR)
counts.matrix = counts.data %>% select(-gene_id) %>% as.matrix()
rownames(counts.matrix) = counts.data$gene_id
dge.data = DGEList(counts=counts.matrix, 
                    group=sample.description$group)
dge.data = calcNormFactors(dge.data, method = "TMM")
plotMDS(dge.data, method = "bcv") 

design = model.matrix(~gt+trt,data = sample.description)
rownames(design) = sample.description$sample
dge.data = estimateGLMCommonDisp(dge.data,design,verbose = TRUE)
dge.data = estimateGLMTrendedDisp(dge.data,design)
dge.data = estimateGLMTagwiseDisp(dge.data,design)
fit = glmFit(dge.data, design)
gt.lrt = glmLRT(fit,coef = "gtR500")
DEgene.gt = topTags(gt.lrt,n = Inf,p.value = 0.01)$table
write.csv(DEgene.gt,"../output/DEgenes.gt.csv")
DEgene.gt.all = topTags(gt.lrt,n = Inf, p.value = 1)$table
write.csv(DEgene.gt.all,"../output/DEgenes.gt.all.csv")

# differentially expressed genes in response to DP treatment (at a FDR < 0.01)
trtDP.lrt = glmLRT(fit, coef = "trtDP")
summary(decideTestsDGE(trtDP.lrt, p.value = 0.01))

# top 9 plot
plotDE = function(genes, dge, sample.description) {
  require(ggplot2)
  tmp.data = t(log2(cpm(dge[genes,])+1))
  tmp.data = tmp.data %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(sample.description,by="sample")
  tmp.data = tmp.data %>%
    pivot_longer(cols=starts_with("Bra"), values_to = "log2_cpm", names_to = "gene")
  pl = ggplot(tmp.data,aes(x=gt,y=log2_cpm,fill=trt))
  pl = pl + facet_wrap( ~ gene)
  pl = pl + ylab("log2(cpm)") + xlab("genotype")
  pl = pl + geom_boxplot()
  pl + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
}
DEgene.gt = topTags(trtDP.lrt,n = Inf,p.value = 0.01)$table
plotDE(rownames(DEgene.gt)[1:9],dge.data,sample.description)

# Gene by treatment interaction    
design.interaction = model.matrix(~gt*trt,data = sample.description)
rownames(design.interaction) = sample.description$sample
dge.data.int = estimateGLMCommonDisp(dge.data,design.interaction,verbose = TRUE)
dge.data.int = estimateGLMTrendedDisp(dge.data,design.interaction)
dge.data.int = estimateGLMTagwiseDisp(dge.data,design.interaction)
fit = glmFit(dge.data.int, design.interaction)
int.lrt = glmLRT(fit, coef = "gtR500:trtDP")
summary(decideTestsDGE(int.lrt, p.value=0.01))
int.gene.gt = topTags(int.lrt,n = Inf,p.value = 0.01)$table
write.csv(int.gene.gt,"interaction_genes.gt.csv")

#top 9 genes that have a significantly different response to treatment in IMB211 as compared to R500.
plotDE(rownames(int.gene.gt)[1:9],dge.data,sample.description)

DEgene.trt = read_csv("DEgenes.trt.csv")
colnames(DEgene.trt)[1] = "GeneID"
DEgene.trt.int = read_csv("DEgenes.interaction.csv")
colnames(DEgene.trt.int)[1] = "GeneID"
B.rapa = read_tsv("FileS9.txt")
colnames(B.rapa)[1] = "GeneID"
DEgene.trt.ann = left_join(DEgene.trt, B.rapa, by = "GeneID")
DEgene.trt.int.ann = left_join(DEgene.trt.int, B.rapa, by = "GeneID")

go.terms = read_tsv("FileS11.txt",col_names=FALSE)
colnames(go.terms) = c("GeneID","GO")
expressed.genes = read_tsv("internode_expressed_genes.txt")
names(expressed.genes) = "GeneID"
gene.lengths = read_tsv("Brapa_CDS_lengths.txt")
gene.lengths.vector = gene.lengths$Length[gene.lengths$GeneID %in% expressed.genes$GeneID]
names(gene.lengths.vector) = gene.lengths$GeneID[gene.lengths$GeneID %in% expressed.genes$GeneID]
expressed.genes.match = expressed.genes[expressed.genes$GeneID %in% names(gene.lengths.vector),]
go.list = strsplit(go.terms$GO,split=",")
names(go.list) = go.terms$GeneID
DE.trt = expressed.genes.match$GeneID %in% DEgene.trt$GeneID
names(DE.trt) = expressed.genes.match$GeneID
DE.trt = as.numeric(DE.trt)
nullp.result = nullp(DEgenes = DE.trt,bias.data = gene.lengths.vector)
rownames(nullp.result) = names(gene.lengths.vector)
GO.out = goseq(pwf = nullp.result,gene2cat = go.list,test.cats=("GO:BP"))
write.table(GO.out[GO.out$over_represented_pvalue < 0.05,1:2],row.names=FALSE,file="GO_terms.txt", quote = FALSE,col.names = FALSE)

gff = import.gff("Brapa_reference/Brapa_gene_v1.5.gff")
mRNAranges = gff[gff$type=="mRNA",c("type", "ID")]
mRNAranges = mRNAranges[str_detect(seqnames(mRNAranges), "Scaffold", negate = TRUE), ]
promoterRanges = flank(mRNAranges, 1500)

Brapaseq = readDNAStringSet("../../assigment_08-imwhaling/input/Brapa_reference/BrapaV1.5_chrom_only.fa")
names(Brapaseq) = str_remove(names(Brapaseq), " \\[.*")
promoters = Brapaseq[promoterRanges]
names(promoters) = promoterRanges$ID
promoters = DNAStringSet(gsub("N","-",promoters))
motifs = read.delim("../input/element_name_and_motif_IUPACsupp.txt",header=FALSE,as.is=TRUE)
motifsV = as.character(motifs[,2])
names(motifsV) = motifs[,1]
motifsSS = DNAStringSet(motifsV)
DEgene.trt.match = DEgene.trt$GeneID[DEgene.trt$GeneID %in% names(promoters)]
expressed.genes.match = expressed.genes$GeneID[expressed.genes$GeneID %in% names(promoters)]
universe.promoters = promoters[expressed.genes.match]
target.promoters = promoters[DEgene.trt.match]
non.target.promoters = universe.promoters[
  ! names(universe.promoters) %in% names(target.promoters)] # all expressed genes not in target set
myMotif = DNAStringSet("YACGTGGC")
target.counts = vcountPDict(myMotif, 
                             subject = target.promoters, 
                             fixed = FALSE) #fixed = FALSE allows IUPAC ambiguity codes
target.counts = 
  vcountPDict(motifsSS[2], 
              subject = target.promoters, 
              fixed = FALSE) +
  vcountPDict(motifsSS[2], 
              subject = reverseComplement(target.promoters), 
              fixed = FALSE)
target.total = rowSums(target.counts>0)
non.target.counts = vcountPDict(motifsSS[2], 
                                 subject = non.target.promoters, 
                                 fixed = FALSE) +
  vcountPDict(motifsSS[2], 
              subject = reverseComplement(non.target.promoters), 
              fixed = FALSE)
non.target.total = rowSums(non.target.counts > 0)
m = matrix(c(
  target.total,
  length(target.promoters) - target.total,
  non.target.total,
  length(non.target.promoters) - non.target.total ),
  ncol=2, 
  dimnames = list(motif=c("Y", "N"), target=c("Y", "N")))
fisher.test(m)
motifEnrichment = function(target.promoters,non.target.promoters,all.counts=F,motifs=motifsSS) {
  
  target.counts = vcountPDict(motifs,target.promoters,fixed=F) + 
    vcountPDict(motifsSS,reverseComplement(target.promoters),fixed=F)
  non.target.counts = vcountPDict(motifs,non.target.promoters,fixed=F) + 
    vcountPDict(motifsSS,reverseComplement(non.target.promoters),fixed=F)
  
  if (all.counts) { 
    target.counts.sum = rowSums(target.counts)
    non.target.counts.sum = rowSums(non.target.counts)
  } else {
    target.counts.sum = rowSums(target.counts > 0)
    non.target.counts.sum = rowSums(non.target.counts > 0)
  }
  n.motifs = length(target.counts.sum)
  results = vector(mode="numeric",length=n.motifs)
  for (i in 1:n.motifs) {
    if (all.counts) {
      m = matrix(c(
        target.counts.sum[i],
        dim(target.counts)[2],
        non.target.counts.sum[i],
        dim(non.target.counts)[2]
      ),ncol=2)
    } else {
      m = matrix(c(
        target.counts.sum[i],
        dim(target.counts)[2]-target.counts.sum[i],
        non.target.counts.sum[i],
        dim(non.target.counts)[2]-non.target.counts.sum[i]
      ),ncol=2)
    }
    results[i] = fisher.test(m,alternative="greater")$p.value
  }
  results.table = data.frame(
    motif=names(motifs),
    non.target.percent = round(non.target.counts.sum/dim(non.target.counts)[2],3)*100,
    target.percent = round(target.counts.sum/dim(target.counts)[2],3)*100,
    p.value =  results)
  results.table = results.table[order(results.table$p.value),]
  results.table
}
motif.results = motifEnrichment(target.promoters, non.target.promoters)
head(motif.results)
