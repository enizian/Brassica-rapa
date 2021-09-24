trt.genes = read_csv("DEgenes.trt.csv") 
brass_voom_E = read_csv("voom_transform_brassica.csv") %>%
  select(GeneID, matches("INTERNODE|PETIOLE|LEAF"))
head(brass_voom_E[,1:6])
brass_voom_E_trt = brass_voom_E %>%
  semi_join(trt.genes)
E_matrix_trt = brass_voom_E_trt %>%
  as.data.frame() %>%
  column_to_rownames("GeneID") %>%
  as.matrix()   
annotation = read_tsv("FileS9.txt",col_names = c("GeneID","description")) 
E_matrix_5 = E_matrix_trt[11:15,]
E_matrix_5_cor = cor(t(E_matrix_5))
E_matrix_5_cor %>% round(3)
diag(E_matrix_5_cor) = 0
E_matrix_5_cor %>% round(2)
E_matrix_5_rank = apply(E_matrix_5_cor,2,function(x) rank(-abs(x)))
E_matrix_5_MR = sqrt(E_matrix_5_rank * t(E_matrix_5_rank))
E_matrix_5_MR %>% round(3)
genes_adj_MR2 = E_matrix_5_MR <= 2
genes_adj_MR2 = 1 * genes_adj_MR2
E_matrix_trt.cor = cor(t(E_matrix_trt))
diag(E_matrix_trt.cor) = 0
E_matrix_trt.cor = apply(E_matrix_trt.cor,2,function(x) rank(-abs(x)))
E_matrix_trt.cor = sqrt(E_matrix_trt.cor * t(E_matrix_trt.cor))
genes_adj_MR4 = E_matrix_trt.cor <= 4
genes_adj_MR4 = genes_adj_MR4 * 1
genes_adj_MR10 = E_matrix_trt.cor <= 10
genes_adj_MR10 = genes_adj_MR10 * 1

gene_graphMR4 = graph.adjacency(genes_adj_MR4, mode = "undirected")
compsMR4 = clusters(gene_graphMR4)$membership
colbar = rainbow(max(compsMR4)+1)
V(gene_graphMR4)$color = colbar[compsMR4+1]
plot(gene_graphMR4, layout = layout_with_fr, vertex.size = 4, vertex.label = NA, main="MR 4")

gene_graphMR10 = graph.adjacency(genes_adj_MR10, mode = "undirected")
compsMR10 = clusters(gene_graphMR10)$membership
colbar = rainbow(max(compsMR10)+1)
V(gene_graphMR10)$color = colbar[compsMR10+1]
plot(gene_graphMR10, layout = layout_with_fr, vertex.size = 4, vertex.label = NA, main="MR 10")

main.cluster = which.max(clusters(gene_graphMR4)$csize)
non.main.vertices = clusters(gene_graphMR4)$membership %>%
  magrittr::extract(. != main.cluster) %>%
  names()
gene_graphMR4 = delete.vertices(gene_graphMR4, non.main.vertices)
distMatrix = shortest.paths(gene_graphMR4, v = V(gene_graphMR4), to = V(gene_graphMR4))
head(distMatrix)[,1:7]
gene1 = match("Bra007662", rownames(distMatrix))
gene2 = match("Bra036631", rownames(distMatrix))
pl = get.shortest.paths(gene_graphMR4, gene1, gene2)$vpath[[1]]
V(gene_graphMR4)[pl]$color = "green"
E(gene_graphMR4)$color = "grey"
E(gene_graphMR4, path = pl)$color = "blue"
E(gene_graphMR4)$width = 1
E(gene_graphMR4, path = pl)$width = 7
V(gene_graphMR4)["Bra007662"]$color = "orange"
V(gene_graphMR4)["Bra036631"]$color = "purple"
plot(gene_graphMR4, layout = layout_with_fr, vertex.size = 5, vertex.label = NA)
