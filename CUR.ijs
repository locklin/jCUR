load'math/misc/svd'
NB. load '~addons/math/lapack/lapack.ijs'
NB. load'~addons/math/lapack/dgesvd.ijs'
NB. svdl =: 3 : 0
NB.  'u s v' =. dgesvd_jlapack_ y
NB.  u;(diag s);v
NB. )
NB. svd =: svdl NB. Use Lapack for now; more robust algorithm


mp=: +/ .*             NB. matrix product
diag =: (<0 1)&|:      NB. diagonal
fcols=: ] {"1~ [: i. [  NB. takes first x columns of y
qrd=: 128!:0            NB. built in QR decomp
DTOL =: (9!:18 '')^0.5



NB. coclass 'jCUR' NB. uncomment this namespace stuff when you're done
NB. CUR_z_ =: CUR_jCUR_


NB. OK: issue with the lapack svd; J is using dgesvd; dgesdd is better
NB. Oddly, the native svd has the right format, and
NB. *svds mx 
NB. -singular value decomposition 
NB. -Returns u;s;v, where s is diagonals
NB. svds =: 3 : 0
NB.  min=. <./ $ y
NB.  'u s v'=. svd y
NB.  (min fcols u);s ; min fcols v
NB. )
svds=: 4 : 0
 'u s v'=. svd y
 x&fcols each u;s;v
)

NB. *<pct> pinv mx
NB. -calculates the partial inverse matrix of mx. Cuts off top 5% of svds by 
NB. -default. This is a weird way of doing things, but this is how MASS::ginv
NB. -works.
NB. The matlab code originally used x=: 0.05  -no idea where that came from
NB. There should be a way of automatically finding sqrt(precision) in J, but 
NB. I have no idea how. This value ties out with the R version.
pinv =: 3 : 0
 DTOL pinv y
:
 'u s v'=. svd y
 posdx =. I. s > >./ 0, x*0{s
 if. 0=#posdx do.
   0 $~ $ y
 else.
   (posdx{"1 v) mp ((1 % posdx{s) * |: posdx {"1 u)
 end.
)

NB. *(k;c;r;t) CUR mx
NB. -mx: m x n matrix
NB. -k: rank parameter with k << min(m,n)
NB. -c: number of columns that we want to select from mx
NB. -r: number of rows that we want to select from A.   
NB. -t: type of column select: 0, for random, 1 random scores with exact number,
NB. -2 ordered scores, exact number
NB. -returns U;C;R
NB. -U: c' x r' matrix
NB. -C: m x c' matrix with c' columns from A. E(c') <= c
NB. -R: r' x n matrix with r' rows form A.    E(r') <= r
NB. -cdx: column index of C in mx
NB. -rdx: row index of R in mx
NB. -Original PNAS paper:
NB. - http://www.pnas.org/content/106/3/697.full.pdf
CUR =: 4 : 0
 'k c r t' =. x
 'u s v'=.  k svds y
 cdx =. (k;c;v;t) Csel y
 rdx =. (k;r;u;t) Csel |: y
 R =. rdx{"2 y
 C =. cdx{"1 y
 U=. (pinv C) mp y mp pinv R
 C;U;R;cdx;rdx
)

NB. *(k;c;v;t) Csel mx
NB. -mx: m x n matrix
NB. -k: rank parameter with k << min(m,n)
NB. -c: number of columns that we want to select from mx
NB. -v: right singular vectors of mx
NB. -t: type 0 1 2 as documented above
NB. -returns the subset columns
Csel=: 4 : 0
 'k c v t'=. x
 'm n'=. $y
 pi =. k %~ +/"1 *: (i.k){"1 v  NB. is the i.k necessary?
 pbj =. 1 <.~ c*pi
 select. t
  case. 0 do. (I. pbj > ?n#0)            NB. random, original
  case. 1 do. c&fcols \: pbj - ?n#0   NB. randomized scores exact num
  case. 2 do. c&fcols  \: pbj              NB. ordered scores exact num
 end.
)

NB. *r SKEL mx
NB. -r is number of columns to select from m x n mx
NB. -returns C;R;G skeleton approximation
NB. - S. A. Goreinov, E. E. Tyrtyshnikov, and N. L. Zamarashkin. 
NB. - A theory of pseudoskeleton approximations. 
NB. - Linear Algebra and its applications, Volume 261, Issues 1-3, August
NB. - 1997, Pages 1-21
NB. -C m x r matrix contains r cols from mx
NB. - R: r x n matrix containing r rows from mx
NB. - G: r x r matrix.
NB. NOT DONE; NEED QR DECOMP PERMUTE MX
NB. 
SKEL=: 4 : 0
 r =.x
 'U S V'=. r svds y
 'Q R'=. qrd |: V
 pv=. '' NB. permutation mx such that pv{mx = Q mp R; in this case, always diag
 C=. (r&focls pv){"1 y
 'Q R'=. qrd |: U
 pu=. '' NB. permutation mx such that pu{mx = Q mp R
 R=. (r&focls pu){"2 y
 G =. (pinv C) mp y mp pinv R
 C;R;G
)

