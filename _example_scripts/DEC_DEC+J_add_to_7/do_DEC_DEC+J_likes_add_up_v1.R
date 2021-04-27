#######################################################
# do_DEC_DEC+J_likes_add_up_v1.R
#
# The purpose of this script is to show that the 
# likelihoods of the DEC and DEC+J models "add up".
#
# That is, for discrete data, if the likelihood
# calculations are valid -- meaning that they are
# valid conditional probabilities, P(data|model) -- 
# then, when likelihoods are calculated for every 
# possible dataset, the likelihoods should add to 1,
# or to 1/(base frequency) if the base frequencies 
# (assumed root state frequencies) are not included
# in the model.
#######################################################


# Load the package (after installation, see above).
library(optimx)         # You need to have some version of optimx available
                        # as it is a BioGeoBEARS dependency; however, if you
                        # don't want to use optimx, and use optim() (from R core) 
                        # you can set:
                        # BioGeoBEARS_run_object$use_optimx = FALSE
                        # ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)

source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)


wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/DEC_DEC+J_add_to_7/write_areas/"
setwd(wd)

geogfn = "geog_3areas.data"

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges


# Make a list of all possible ranges at one tip
# (including the null range)
# 3 areas means 2^3 = 8 possible ranges
possible_ranges = c(
"000",
"100",
"010",
"001",
"110",
"011",
"101",
"111")

# Of course, for 2 tips, there will be 8x8=64 possible data patterns.  
# So, we need to make a geography file for each one.

new_fns = NULL
for (i in 1:length(possible_ranges))
	{
	for (j in 1:length(possible_ranges))
		{
		sp1_range = possible_ranges[i]
		sp2_range = possible_ranges[j]
		tipranges_new = tipranges
		tipranges_new@df[1,] = strsplit(sp1_range, split="")[[1]]
		tipranges_new@df[2,] = strsplit(sp2_range, split="")[[1]]
		new_fn = paste0("geog_", sp1_range, "_", sp2_range, ".data")
		new_fns = c(new_fns, new_fn)
		save_tipranges_to_LagrangePHYLIP(tipranges_object=tipranges_new, lgdata_fn=new_fn)
		} # END for (j in 1:length(possible_ranges))
	} # END for (i in 1:length(possible_ranges))


# Let's get the likelihoods under the DEC model,
# for some particular parameters
# (here, d and e are very low, allowing us to see 
# a pure cladogenetic-range-change model

# DEC likelihoods
trfn = "tree.newick"
lnLs = NULL
for (i in 1:length(new_fns))
	{
	max_range_size = 3
	trfn = trfn
	geogfn = new_fns[[i]]
	
	tr = read.tree(trfn)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
	
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$trfn = trfn
	BioGeoBEARS_run_object$geogfn = geogfn
	BioGeoBEARS_run_object$max_range_size = max_range_size
	BioGeoBEARS_run_object$min_branchlength = 0.000001
	BioGeoBEARS_run_object$include_null_range = TRUE
	BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
	BioGeoBEARS_run_object$return_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_ancprobs = TRUE
	BioGeoBEARS_run_object$num_cores_to_use = 1
	BioGeoBEARS_run_object$allow_null_tips = TRUE

	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = 0.00001
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = 0.00001
	
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = 0.00001
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = 0.00001
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "type"] = "fixed"
	 
	
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE)	
	lnL = res
	lnLs = c(lnLs, lnL)
	}

DEC_lnLs = lnLs
DEC_lnLs
exp(DEC_lnLs)

tab_DEC = cbind(new_fns, round(DEC_lnLs,4), round(exp(DEC_lnLs),3))
tab_DEC = as.data.frame(tab_DEC)
tab_DEC = dfnums_to_numeric(tab_DEC)
tab_DEC

table_DEC_dLow_eLow = tab_DEC
write.table(table_DEC_dLow_eLow, file="table_DEC_dLow_eLow.txt", sep="\t", quote=FALSE)

# What do the likelihoods of all the data add up to?
sum(exp(DEC_lnLs))
# 7
# (works -- 7 possible ancestral states, excluding null range)

# They don't add up to 8, because the probability of any possible tip data,
# when the ancestor has range null, is zero.  I guess technically I could
# have programmed the model such that an ancestral null range is "transmitted" 
# through the cladogenesis model and up the branches, producing null
# ranges at all tips with probability 1...but what's the point?



# Let's repeat for some other values of d and e
# d = 0.2
# e = 0.1
# DEC likelihoods
trfn = "tree.newick"
lnLs = NULL
for (i in 1:length(new_fns))
	{
	max_range_size = 3
	trfn = trfn
	geogfn = new_fns[[i]]
	
	tr = read.tree(trfn)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
	
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$trfn = trfn
	BioGeoBEARS_run_object$geogfn = geogfn
	BioGeoBEARS_run_object$max_range_size = max_range_size
	BioGeoBEARS_run_object$min_branchlength = 0.000001
	BioGeoBEARS_run_object$include_null_range = TRUE
	BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
	BioGeoBEARS_run_object$return_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_ancprobs = TRUE
	BioGeoBEARS_run_object$num_cores_to_use = 1
	BioGeoBEARS_run_object$allow_null_tips = TRUE

	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = 0.2
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = 0.2
	
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = 0.1
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = 0.1
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "type"] = "fixed"
	 
	
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE)	
	lnL = res
	lnLs = c(lnLs, lnL)
	}

DEC_lnLs = lnLs
DEC_lnLs
exp(DEC_lnLs)

tab_DEC = cbind(new_fns, round(DEC_lnLs,4), round(exp(DEC_lnLs),3))
tab_DEC = as.data.frame(tab_DEC)
tab_DEC = dfnums_to_numeric(tab_DEC)
tab_DEC

table_DEC_dMod_eMod = tab_DEC
write.table(table_DEC_dMod_eMod, file="table_DEC_dMod_eMod.txt", sep="\t", quote=FALSE)


# What do the likelihoods of all the data add up to?
sum(exp(DEC_lnLs))
# 7
# Yep, 7 again...






# DEC+J likelihoods, with an intermediate model
trfn = "tree.newick"
lnLs = NULL
for (i in 1:length(new_fns))
	{
	max_range_size = 3
	trfn = trfn
	geogfn = new_fns[[i]]
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$trfn = trfn
	BioGeoBEARS_run_object$geogfn = geogfn
	BioGeoBEARS_run_object$max_range_size = max_range_size
	BioGeoBEARS_run_object$min_branchlength = 0.000001
	BioGeoBEARS_run_object$include_null_range = TRUE
	BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
	BioGeoBEARS_run_object$return_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_ancprobs = TRUE
	BioGeoBEARS_run_object$num_cores_to_use = 1
	BioGeoBEARS_run_object$allow_null_tips = TRUE

	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = 0.2
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = 0.2
	
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = 0.1
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = 0.1
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "type"] = "fixed"
	 
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] = 0.15
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] = 0.15
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] = "fixed"
	
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE)	
	lnL = res
	lnLs = c(lnLs, lnL)
	}

DECj_lnLs = lnLs
exp(DECj_lnLs)

tab_DECj = cbind(new_fns, round(DECj_lnLs,4), round(exp(DECj_lnLs),3))
tab_DECj = as.data.frame(tab_DECj)
tab_DECj = dfnums_to_numeric(tab_DECj)
tab_DECj

table_DECj_dMod_eMod_jMod = tab_DECj
write.table(table_DECj_dMod_eMod_jMod, file="table_DECj_dMod_eMod_jMod.txt", sep="\t", quote=FALSE)


sum(exp(DECj_lnLs))
# 7
# 7, again...

sum(exp(DEC_lnLs))
sum(exp(DECj_lnLs))

write.table(tab_DEC, file="tab_DEC.txt", sep="\t")
write.table(tab_DECj, file="tab_DECj.txt", sep="\t")


# DEC+J likelihoods, with an almost-all-j model
trfn = "tree.newick"
lnLs = NULL
for (i in 1:length(new_fns))
	{
	max_range_size = 3
	trfn = trfn
	geogfn = new_fns[[i]]
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object$trfn = trfn
	BioGeoBEARS_run_object$geogfn = geogfn
	BioGeoBEARS_run_object$max_range_size = max_range_size
	BioGeoBEARS_run_object$min_branchlength = 0.000001
	BioGeoBEARS_run_object$include_null_range = TRUE
	BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
	BioGeoBEARS_run_object$return_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
	BioGeoBEARS_run_object$calc_ancprobs = TRUE
	BioGeoBEARS_run_object$num_cores_to_use = 1
	BioGeoBEARS_run_object$allow_null_tips = TRUE

	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = 0.00001
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = 0.00001
	
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = 0.00001
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = 0.00001
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "type"] = "fixed"
	 
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] = 2.9
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] = 2.9
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] = "fixed"
	
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE)	
	lnL = res
	lnLs = c(lnLs, lnL)
	}

DECj_lnLs = lnLs
exp(DECj_lnLs)

tab_DECj = cbind(new_fns, round(DECj_lnLs,4), round(exp(DECj_lnLs),3))
tab_DECj = as.data.frame(tab_DECj)
tab_DECj = dfnums_to_numeric(tab_DECj)
tab_DECj

table_DECj_dLow_eLow_jMod = tab_DECj
write.table(table_DECj_dLow_eLow_jMod, file="table_DECj_dLow_eLow_jMod.txt", sep="\t", quote=FALSE)


sum(exp(DECj_lnLs))
# 7



# CONCLUSION: DEC and DEC+J both produce the same total likelihood 
# when calculated over all possible datasets.  This adds up to 7, because 
# there are 7 possible ancestral states (excluding the eighth state, the 
# null range).
# If we added "equal base frequencies" to the root probability, and said 
# each of possible ancestral states had probability 1/7, then the likelihoods
# would add to 1.

# CONCLUSION #2: The BioGeoBEARS likelihoods, times 1/7, sum to 1 under both 
# models. Therefore the probabilities of any particular dataset are valid 
# conditional probabilities.  Therefore, the conditional probability that each 
# model confers on a single dataset can be compared using standard tools of 
# model choice (AICc, etc.).
sum(exp(DEC_lnLs))
sum(exp(DECj_lnLs))










