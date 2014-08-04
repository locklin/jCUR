jCUR
====

attempt at CUR decomposition

Following the rCUR package as well as the original Matlab (included).
If I can get it to work, this might be useful for use with the Jd database.

Original PNAS paper:
http://www.pnas.org/content/106/3/697.full.pdf

if this works and is useful, I'll add some other matrix decomposition methods
(NNMF, Nuclear Norm regularization, etc).

Presently, it is completely untested, though the #2 column select technique
ties out with the rCUR package in R.


***********
The SVD function is not super efficient for large scale problems.
In another approximate mx decomp package, the PROPACK lib was 
a suggested approach in favor of the LAPACK methods. This uses
Lanczos bidiagonalization algorithm with partial reorthogonalization (BPRO).
Lanczos is simple anyway. Claimed factor of 4 speedup, and 1/2 memory load.
http://sun.stanford.edu/~rmunk/PROPACK/paper.ps.gz
