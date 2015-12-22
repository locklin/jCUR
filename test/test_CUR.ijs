NB. "GIST","GIST","GIST","GIST","GIST","GIST","GIST","GIST","GIST","GIST","LEIO",
NB. "LEIO","LEIO","LEIO","LEIO","LEIO","LEIO","LEIO","LEIO","LEIO","LEIO","LEIO",
NB. "SARC","SARC","SARC","SARC","SARC","SARC","SARC","SARC","SARC"
NB.
NB.##########################################
NB. compare to R
NB. CUR(STTm,c=31,r=21,k=3,method="top.scores") -> cur
loc=. 3 : '> (4!:4 <''y'') { 4!:3 $0'  
PATH=: getpath_j_ loc''
DATA=: PATH,'../data/'
DATA=:'/home/scott/src/jCUR/data/'

STTm =: ".;._2 fread DATA,'STTm.csv'
'C U R a b'=. (3;31;21;2)CUR STTm

bb=. 4531 4634 4610 4619 4693 4620 2124 5262 2818 2888 4633 2884 4499 4501 4628 4577 2122 4603 4632 2123 4576
aa=. 2 5 10 12 13 16 24 25 26 27 28 6 21 7 1 3 29 8 20 11 0 15 30 18 19 22 4 17 9 14 23



NB. c r k
NB. 31 21 3 
NB. k c r
