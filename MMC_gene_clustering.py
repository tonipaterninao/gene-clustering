from rlscore.learner.mmc import MMC
from rlscore.reader import read_sparse

## Import the dataset
gene_data_na = read_sparse("./gene_data_na.txt")

## Build the model
kwargs = {}
kwargs["X"] = gene_data_na
kwargs["regparam"] = 1
kwargs["kernel"] = "GaussianKernel"
kwargs["number_of_clusters"] = 4    ## Set the number of clusters found with the eigengap method for this kernel
learner = MMC(**kwargs)
labels = learner.results

# Write the results in output file
# out = open("python_clustering.out","w")
# for label in labels["predicted_clusters_for_training_data"]:
#    out.write(str(label) + "\n")
# out.close()