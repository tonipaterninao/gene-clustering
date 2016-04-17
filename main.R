
library("kernlab")

# Load the gene expression dataset (YORF column is rownames)
gene_data = read.delim("DNAcomplete_dataset.cdt", header = TRUE, row.names = 1)

# Separate the gene description from the numerical data - store in a separate file
# erase to save memory space
gene_info = gene_data[,1]
gene_data = gene_data[,2:ncol(gene_data)]
#write(levels(gene_info), sep="\t", file = "gene_names")

## Incorporate gene information (after editing the file with a text editor using regular expressions)
gene_info = read.delim("./gene_names", fill = TRUE, na.strings = "?");

gene_data_info = cbind(rownames(gene_data_na),gene_data_na)
colnames(gene_data_info)[1] = "Gene"
gene_data_info = merge(gene_info,gene_data_info,by="Gene")

## select biological processes with a representation higher than 10
## (at least 5 genes per category)

func.count = table(gene_data_info$Bio..process)
length(func.count)        ## Visualize the total biological processes in the subset

func.count = func.count[func.count > 10]
length(func.count)       ## Visuallize the nnumber of processes after filtering

## Visualize the fractions of biological processes
barplot(sort(func.count, decreasing = TRUE),
        col = colors()[c(61,57,60,58)],
        ylim = c(0,1000),
        xlab = "Biological function",
        ylab = "Gene count",
        main = "Gene distribution according to biological process")

###########################################################################################
################################### Spectral clustering ###################################
###########################################################################################

########### Find the optimal number of centers (k) with the eigengap criterion ############

## Compute the adjacency matrix - use the default kernel used by method specc

## Following calculations need a dataset with no NA values

gene_data_na = na.exclude(gene_data)
gene_data_na = gene_data_na[,2:ncol(gene_data_na)]            # remove gweight column
gene_data_na = gene_data_na[-c(8,16,20,29,34,43,47,51,55,59,62,63)]  # remove cols. of zeros

## Compute the distances

gaussian_kernel = rbfdot(sigma = 1);  # default value of sigma = 1

A = kernelMatrix(gaussian_kernel, as.matrix(gene_data_na))
diag(A) = rep(0, nrow(x = gene_data_na))

## Compute de degree matrix

D = matrix(rep(0, nrow(gene_data_na) ** 2), ncol = nrow(gene_data_na))
d = c()    ## intermediate vector of degrees

for (i in 1:nrow(gene_data_na)) {
  degree = sum(A[i])
  d = c(d,degree)
}

diag(D) = d

## Compute the Laplacian matrix

L = D - A

## erase intermediate matrices to economize resources

rm(A,D,d)

## Compute the eigenvalues and eigenvectors

p = eigen(L)
p$values =sort(p$values)   ## sort in increasing order
rm(L)

## Compute the optimal number of clusters (k) with the eigengap value

max_eigengap = 0
for (j in 3:length(p$values)){
  eigengap = abs(p$values[j] - p$values[j-1])
  if (eigengap > max_eigengap){
    max_eigengap = eigengap
    k_opt = j
  }
}

###################### Apply spectral clustering with a single kernel #####################

single_kernel_specc = specc(as.matrix(gene_data_na), centers = k_opt)
plot(single_kernel_specc)

## Dimension reduction of the dataset by PCA
s_pca = prcomp(gene_data_na)

## Visuallization of the two components that allow better distinction of datapoints

plot(s_pca$x[,1],s_pca$x[,2],col = single_kernel_specc@.Data,
     pch = 16,
     main = "PC visualization of genes according to spectral clustering",
     ylab = "PC2", xlab= "PC1")
legend(-11.5,4,c(1:k_opt),fill = c(1:k_opt),title="Cluster")

############## Apply spectral clustering with a linear combination of kernels #############

## Linear combination of kernels using fixed rule
##    Using: Gaussian, polynomial, and Laplacian kernels

K = (kernelMatrix(kernel = rbfdot(), x = as.matrix(gene_data_na)) 
+ kernelMatrix(kernel = polydot(), x = as.matrix(gene_data_na)) 
+ kernelMatrix(kernel = laplacedot(), x = as.matrix(gene_data_na))) /3

## Compute de degree matrix

D = matrix(rep(0, nrow(gene_data_na) ** 2), ncol = nrow(gene_data_na))
d = c()    ## intermediate vector of degrees

for (i in 1:nrow(gene_data_na)) {
  degree = sum(K[i])
  d = c(d,degree)
}

diag(D) = d

## Compute the Laplacian matrix

L = D - K

## erase intermediate matrices to economize resources

rm(D,d)

## Compute the eigenvalues and eigenvectors

p = eigen(L)
p$values =sort(p$values)   ## sort in increasing order
rm(L)

## Compute the optimal number of clusters (k) with the eigengap value

max_eigengap = 0
for (j in 3:length(p$values)){
  eigengap = abs(p$values[j] - p$values[j-1])
  if (eigengap > max_eigengap){
    max_eigengap = eigengap
    k_opt_lc = j
  }
}

## Apply spectral clustering using the linear combination

linear_comb_specc = specc(as.matrix(gene_data_na), centers = k_opt_lc, kernel = K)

## Visualize the clustering

plot(s_pca$x[,1],s_pca$x[,2],col = linear_comb_specc@.Data,
     pch = 16,
     main = "PC visualization of genes according to spectral clustering\n(with linear combination of clusters)",
     ylab = "PC2", xlab= "PC1")
legend(-11.5,4,c(1:k_opt_lc),fill = c(1:k_opt_lc),title="Cluster")

## Create file to perform MMC with python

export = c()

for (i in 1:nrow(gene_data_na)){
  for (j in 1:ncol(gene_data_na)){
    val = gene_data_na[i,j]
    val = paste(j-1,val,sep=":")
    export = c(export,val)
  }
}

export = matrix(data = export,ncol = ncol(gene_data_na))

write(export, sep = "\t", file = "gene_data_na.txt", ncolumns = 52)

rm(export)

## Parse the python results from MMC

mmc_labels = read.table("./python_clustering.out")

## Convert to cluster labels

mmc = as.matrix(mmc_labels)            ## intermediate variable
mmc_labels = c()
clusters = c(1,2,3,4)
for (i in 1:nrow(mmc)){
  c = clusters[as.vector(mmc[i,]) == 1];
  mmc_labels = c(mmc_labels,c)
}

## Visualize the clustering with MMC
plot(s_pca$x[,1],s_pca$x[,2],col = mmc_labels,
     pch = 16,
     main = "PC visualization of genes according to MMC\n(Gaussian kernel)",
     ylab = "PC2", xlab= "PC1")
legend(-11.5,4,c(1,2,3,4),fill = c(1,2,3,4),title="Cluster")

## Plot the sizes of the clusters found by each method
mmc_size = c(232,1839,51,100)

size_table = matrix(c(single_kernel_specc@size,linear_comb_specc@size,0,mmc_size),
                    ncol = 4,
                    byrow = TRUE)
barplot(t(size_table), names.arg = c("spectral","mk-spectral","MMC"),
        ylim = c(0,2500), main = "Sizes of the clusters found by each method",
        xlab = "Clustering method")

###########################################################################################
################################### Biological processes ##################################
###########################################################################################

## Single kernel spectral clsutering

sk_cluster_1 = gene_data_info[single_kernel_specc@.Data == 1,]
sk_cluster_1 = table(sk_cluster_1$Bio..process)
sk_cluster_1 = sort(sk_cluster_1[sk_cluster_1 > 5], decreasing = TRUE)
barplot(sk_cluster_1[2:length(sk_cluster_1)], col = colors()[62],
        names.arg = c("transcr.","prot. synth.","secret.","transport",
                      "glycol.","c.cycle","c.wall.","nuc.prot"),
        cex.names = 0.7,
        main="Biological processes in cluster 1")

sk_cluster_2 = gene_data_info[single_kernel_specc@.Data == 2,]
sk_cluster_2 = table(sk_cluster_2$Bio..process)
sk_cluster_2 = sort(sk_cluster_2[sk_cluster_2 > 5], decreasing = TRUE)
barplot(sk_cluster_2[2:length(sk_cluster_2)], col = colors()[62],
        names.arg = c("p. synth","trnscr.","trnsp.","splic.","secret.",
                      "cytosk.","c. wall","p.degr."),
        cex.names = 0.7,
        ylim = c(0,40),
        main="Biological processes in cluster 2")

sk_cluster_3 = gene_data_info[single_kernel_specc@.Data == 3,]
sk_cluster_3 = table(sk_cluster_3$Bio..process)
sk_cluster_3 = sort(sk_cluster_3[sk_cluster_3 > 5], decreasing = TRUE)
barplot(sk_cluster_3[2:length(sk_cluster_3)], col = colors()[62],
        names.arg = c("prot. synth.","transcript.","secretion","cell cycle","glycolysis"),
        ylim = c(0,30),
        main="Biological processes in cluster 3")

sk_cluster_4 = gene_data_info[single_kernel_specc@.Data == 4,]
sk_cluster_4 = table(sk_cluster_4$Bio..process)
sk_cluster_4 = sort(sk_cluster_4[sk_cluster_4 > 10], decreasing = TRUE)
barplot(sk_cluster_4[2:length(sk_cluster_4)], col = colors()[62],
        names.arg = c("prot. synth.","secret.","transcript.","cytosk.","transport","pr. degr."),
        ylim = c(0,60),
        main="Biological processes in cluster 4")

## MK spectral clustering

mk_cluster_1 = gene_data_info[linear_comb_specc@.Data == 1,]
mk_cluster_1 = table(mk_cluster_1$Bio..process)
mk_cluster_1 = sort(mk_cluster_1[mk_cluster_1 > 10], decreasing = TRUE)
barplot(mk_cluster_1[2:length(mk_cluster_1)], col = colors()[62],
        names.arg = c("prot. synt.","transcr.","secretion",
                      "transp.","cytosk.","chromatin"),
        main="Biological processes in cluster 1")

mk_cluster_2 = gene_data_info[linear_comb_specc@.Data == 2,]
mk_cluster_2 = table(mk_cluster_2$Bio..process)
mk_cluster_2 = sort(mk_cluster_2[mk_cluster_2 > 5], decreasing = TRUE)
barplot(mk_cluster_2[2:length(mk_cluster_2)], col = colors()[62],
        names.arg = c("pr. synth.","trnsc.","secret.","transp.",
                      "c. cycle","cytosk.","glycol"),
        #ylim = c(0,40),
        main="Biological processes in cluster 2")

mk_cluster_3 = gene_data_info[linear_comb_specc@.Data == 3,]
mk_cluster_3 = table(mk_cluster_3$Bio..process)
mk_cluster_3 = sort(mk_cluster_3[mk_cluster_3 > 10], decreasing = TRUE)
barplot(mk_cluster_3[2:length(mk_cluster_3)], col = colors()[62],
        names.arg = c("p.synt.","trnsc.","secrt.","trnsp.",
                      "p. degr.","splc.","cytosk.","c. wall","glycosyl"),
        cex.names = 0.75,
        ylim = c(0,60),
        main="Biological processes in cluster 3")