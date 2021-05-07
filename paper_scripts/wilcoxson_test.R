### t.test ###

### data ###
alg = c(
4.362,
2.953,
2515.934,
1118350,
1141.595,
1340.290,
153703251,
5556585
)

LR = c(
9.009,
11.410,
2588.491,
17700169,
6885.041,
2075.618,
168134153,
32326913
)

MAR = c(
35.657,
6.947,
25275.660,
1970072,
1764.455,
2088.561,
337368866,
31730271
)

MAR_logTrans = c(
55.187,
6.053,
60742.360,
1749733,
6885.041,
3264.232,
214827228,
53749339
)

cbind(alg,LR,MAR)
cbind(LR,MAR,LR-MAR)
cbind(alg,MAR,alg-MAR)


hist(MAR)
plot(alg[1:6],MAR[1:6])
plot(LR[1:6],MAR[1:6])
plot(LR[1:6],alg[1:6])

?wilcox.test

wilcox.test(LR,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: LR is equal or higher than MAR	
wilcox.test(alg,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: alg is iqual or higher than MAR
wilcox.test(alg,MAR_logTrans,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)	# H0: alg is iqual or higher than MAR_logTrans
wilcox.test(alg,LR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: alg is equal or higher than LR	
wilcox.test(MAR,MAR_logTrans,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)	# H0: MAR is iqual or higher than MAR_logTrans
wilcox.test(MAR_logTrans,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)	# H0: MAR_logTrans is iqual or higher than MAR


wilcox.test(LR,alg,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)




##### Pars absDif

parLR = c(
1.97538,
0.24117,
1.35692,
0.88350,
0.00021,
0.20832,
1.06975,
0.19183,
0.82403,
0.18794
)

parALG = c(
2.68019,
0.32537,
46.60082,
0.66838,
0.00038,
0.15534,
0.09820,
0.60795,
0.36212,
0.97794
)

wilcox.test(parLR,parALG,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: parLR is equal or higher than parALG	




##### noisy vs replicates

noisy = c(
9.009,
4.362,
35.657,
55.187,
6885.041,	
1141.595,	
1764.455,	
6885.041
)

reps = c(
11.410,
2.953,	
6.947,
6.053,
2075.618,	
1340.290,	
2088.561,	
3264.232
)

wilcox.test(reps,noisy,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: reps is equal or higher than noisy	



