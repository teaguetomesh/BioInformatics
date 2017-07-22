Teague Tomesh 11/22/2016

These scripts were written as part of the course work for CS 576: Introduction to Bioinformatics.

The purpose of this exercise is to implement a working Hidden Markov Model using the Viterbi Algorithm and Maximum Likelihood Estimations. This model is then used to predict the sequence of Exons (segment of DNA that is actually used to create a protein) and Introns (segment that is not used to create proteins) in a given sequence of DNA. This is a relevant problem in Bioinformatics where sequencing techniques can determine the actual sequence of A's, T's, G's, and C's, but it is more difficult to determine where the exons and introns reside. 

To use these scripts, begin by running predict_exons.py:

>>> ./predict_exons.sh 'path to TRAIN file' 'path to TEST file'

This will produce an OUT file which can be compared to the TRUTH file using exon_accuracy.py:

>>> ./exon_accuracy.sh 'path to TRUTH file' 'path to OUT file'

This will generate an output that represents the accuracy, recall, and precision of the generated exon and intron sequence.