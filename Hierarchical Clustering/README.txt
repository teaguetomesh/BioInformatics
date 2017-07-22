Teague Tomesh 12/9/2016

As part of the course work for CS 576: Introduction to Bioinformatics

The scripts in this folder can implement a variety of agglomerative hierachical clustering algorithms to cluster given yeast genes into similar groups. The purpose of this assignment is to become familiar with the structure and implementation of the single-link, average-link, and complete-link clustering algorithms and clustering techniques in general. This particular algorithm uses Euclidean distance to compute the distance between clusters.

To run the cluster script, specify the gene expression data, the type of hierarchical clustering ('S' single-link, 'C' complete-link, 'A' average-link), and the number of clusters to return:

>>> ./cluster.sh 'small-yeast-set' S 2

This will output the different genes in each of the clusters as well as the average expression value for that cluster.