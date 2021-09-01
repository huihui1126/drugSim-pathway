# drugSim-pathway
Drug similarity evaluation based on pathway fingerprints. We proposed an approach using pathway fingerprint similarity based on a “drug-target-pathway” heterogeneous network to explore the mechanisms of natural drugs XYPI using NSAIDs and GCs as reference agents. A heterogeneous network similarity algorithm used for social networks was applied to investigate the pathway fingerprint similarity between drugs in the heterogeneous "drug-target-pathway" network and to predict the drug-pathway associations of XYPI.

This programme is written by Feiei Guo(ffguo@icmm.ac.cn).

The main programe is extracRelation.pl, which call the pathsim.py and hclu.r for calculation for two drug similarity based on pathway, and hierachical clustering based on drug similarity matrix.

Usage: 
       
       Perl extracRelation.pl $suffixname $randtime $pathdb
       eg., Perl extracRelation.pl test 10 GO

Input: Combine_score.txt is a file including drug-target interaction of XYPI, Glucocorticoids (GCs) and Nonsteroidal Anti-inflammatory Drugs (NSAIDs). The first three columns (Drug, Target, Combine_score) are required.

Output1: simMatrix-$suffixname.txt is a similarity matrix of drugs,including XYPI, GCs and NSAIDs.

Output2: clu-$suffxiname.pdf is a hierachcial clustering of drugs(XYPI, GCs and NSAIDs) base on simMatrix-$suffixname.txt.

Requirement

1.Installation for Perl (v5.24.1 and above version)
https://www.activestate.com/products/perl/downloads/

Perl package

List::Util

2.Installation for python (v3.6.4 and above version)
https://www.python.org/downloads/

Python packages

pandas
numpy
sys

3.Installation for R
https://cran.r-project.org/bin/windows/base/

R package

ape
