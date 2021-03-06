# DCMD

## Data process
implementation/data_propocess.R is used to propocess the raw data and generate the training and test set. 

## Implementation
To implement DCMD to real data, run implementation/real_cancer.R for one fold. You will need to change the path to the path of your data and change out.save.dir according to the path you want to save your file.

To run 10-fold CV in parallel, submit /implementation/realdata_submit.sh in the cluster.

# DCMD: Distance-based Classification Using Mixture Distributions on Microbiome Data

Konstantin Shestopaloff, Mei Dong, Wei Xu

DCMD aims to improve classification performance when using microbiome community data, where the predictors are composed of sparse and heterogeneous count data. This approach models the inherent uncertainty in sparse counts by estimating a mixture distribution for the sample data, and representing each observation as a distribution, conditional on observed counts and the estimated mixture, which are then used as inputs for distance-based classification. The method is implemented into a k-means and k-nearest neighbours framework and we identify two distance metrics that produce optimal results. 
