brass_voom_E = read_csv("voom_transform_brassica.csv")
calc.cv = function(x, na.rm=TRUE) {
  if(na.rm==TRUE) x = na.omit(x)
  result = sd(x) / mean(x)
  result = abs(result)
  return(result)
}
brass_voom_E = brass_voom_E %>%
  rowwise() %>%
  mutate(cv = calc.cv(c_across(-GeneID))) %>%
  ungroup() %>%
  select(GeneID, cv, everything())
brass_voom_E_1000 = brass_voom_E %>% arrange(desc(cv)) %>% slice(1:1000) %>% arrange(GeneID)
E_matrix = brass_voom_E_1000 %>% 
  select(-GeneID, -cv) %>% 
  as.matrix() %>% 
  scale() 
gene_hclust_col = E_matrix %>% t()  %>% dist() %>% hclust()
ggdendrogram(gene_hclust_col)

plot(gene_hclust_col, cex=.6)
rect.hclust(gene_hclust_col, k = 3, border = "red")
plot(gene_hclust_col, cex=.6)
rect.hclust(gene_hclust_col, k = 5, border = "blue")
plot(gene_hclust_col, cex=.6)
rect.hclust(gene_hclust_col, k = 7, border = "green")

library(pvclust)
set.seed(12456)
fit = pvclust(E_matrix, method.hclust = "ward.D", method.dist = "euclidean", nboot = 50)
plot(fit, print.num=FALSE)
fit = pvclust(E_matrix, method.hclust = "ward.D", method.dist = "euclidean", nboot = 1000)
plot(fit, print.num=FALSE)

library(gplots)
heatmap.2(as.matrix(cities), Rowv=as.dendrogram(cities_hclust), scale="row", density.info="none", trace="none")
gene_hclust_row = E_matrix %>% dist() %>% hclust()
heatmap.2(E_matrix, Rowv = as.dendrogram(gene_hclust_row),  density.info="none", trace="none", margins = c(10,5))

library(ggplot2)
prcomp_counts = prcomp(t(E_matrix))
scores = as.data.frame(prcomp_counts$rotation)[,c(1,2)]
set.seed(25)
fit = kmeans(E_matrix, 9)
clus = as.data.frame(fit$cluster)
names(clus) = paste("cluster")
plotting = merge(clus, scores, by = "row.names")
plotting$cluster = as.factor(plotting$cluster)
# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 2, stat = "identity") 
fit = kmeans(E_matrix, 3)
clus = as.data.frame(fit$cluster)
names(clus) = paste("cluster")
plotting = merge(clus, scores, by = "row.names")
plotting$cluster = as.factor(plotting$cluster)
# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 2, stat = "identity")
fit = kmeans(E_matrix, 15)
clus = as.data.frame(fit$cluster)
names(clus) = paste("cluster")
plotting = merge(clus, scores, by = "row.names")
plotting$cluster = as.factor(plotting$cluster)
# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 2, stat = "identity") 

# Gap statistic clustering
library(cluster)
set.seed(125)
gap = clusGap(E_matrix, FUN = kmeans, iter.max = 30, K.max = 20, B = 100)
plot(gap, main = "Gap Statistic")
#diminishing returns at around k = 7

with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
fit = kmeans(E_matrix, 7)
clus = as.data.frame(fit$cluster)
names(clus) = paste("cluster")
plotting = merge(clus, scores, by = "row.names")
plotting$cluster = as.factor(plotting$cluster)
# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 2, stat = "identity") 
with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
fit = kmeans(E_matrix, 8)
clus = as.data.frame(fit$cluster)
names(clus) = paste("cluster")
plotting = merge(clus, scores, by = "row.names")
plotting$cluster = as.factor(plotting$cluster)
# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 2, stat = "identity") 

