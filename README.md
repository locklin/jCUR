jCUR
====

attempt at CUR decomposition

Following the rCUR package as well as the original Matlab (included).
This might eventually be useful for use with the Jd database.

Original PNAS paper:
http://www.pnas.org/content/106/3/697.full.pdf

If this works and I can make it do something useful, I'll add some other matrix decomposition methods
(NNMF, Nuclear Norm regularization, etc) for eventual inclusion into pacman.

Presently, it is completely untested, though the #2 column select technique
ties out with the rCUR package in R.


***********
The SVD function is not super efficient for large scale problems.
In another approximate mx decomp package, the PROPACK lib was 
a suggested approach in favor of the LAPACK methods. This uses
Lanczos bidiagonalization algorithm with partial reorthogonalization (BPRO).
Lanczos is simple anyway. Claimed factor of 4 speedup, and 1/2 memory load.
http://sun.stanford.edu/~rmunk/PROPACK/paper.ps.gz

R presently uses dgesdd and zgesdd, so as an intermediate step, building the
interface for that might help things.

***********
The test data set is taken from R, which is in turn taken from the Gene Expression 
Omnibus database (GSE3443) from the Stanford microarray DB.