library(BioGeoBEARS)

#
# Matrix from 4area_Qmat_v1.xlsx
#
tmpmat = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, -3, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, -3, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
0, 0, 0, -3, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0,
0, 0, 0, 0, -3, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 1, 0, 1, 0, 0,
0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 1, 1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 0, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

Qmat = matrix(tmpmat, nrow=16, ncol=16, byrow=TRUE)

# Make a series of t values
tvals = seq(0.68, 0.71, 0.0001)

# Exponentiate each with EXPOKIT's dgpadm (good for small dense matrices)
for (t in tvals)
	{
	Pmat = expokit_dgpadm_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
	cat("\n\nTime=", t, "\n", sep="")
	
	# Location of the AB->ABC transition
	startstate = 6
	endstate = 12
	
	# Starting in AB
	starting_vector = rep(0, times=16)
	starting_vector[startstate] = 1
	
	ending_vector = starting_vector %*% Pmat
	cat(ending_vector, sep=" ")
	cat("\n")
	#print(Pmat)
	}



