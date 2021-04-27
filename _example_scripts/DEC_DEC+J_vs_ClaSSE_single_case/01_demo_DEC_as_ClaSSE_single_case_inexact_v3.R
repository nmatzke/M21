
#######################################################
# 
# This script demonstrates how the biogeography models
# DEC and DEC+J can be shown to be special cases of 
# the ClaSSE model (Goldberg & Igic 2012). That is, an 
# appropriately parameterized ClaSSE will produce a
# likelihood calculation identical to that of DEC or DEC+J.
# 
# Specifically, if one assumes:
# 
# 1. A Yule process (pure-birth) produced the tree
# 
# 2. Speciation and extinction rate are state-independent
#
# Details include:
#
# 3. DEC and DEC+J do not calculate the likelihood of the tree
#    under a birth/death model. But, given assumptions 1 & 2,
#    it is trivial to calculate the tree likelihood under a 
#    Yule process and either add this to DEC/DEC+J, or 
#    subtract it from ClaSEE.
#
# 4. DEC and DEC+J do not use an extra likelihood calculation
#    for the states at the root. This either has to be added
#    to DEC/DEC+J, or subtracted from the ClaSSE likelihood.
# 
# 5. To make the models equivalent, we have to include the 
#    geographic range of "null" (living in no areas) as a 
#    state in the state space.  Peculiarly, this means that 
#    a lineage could keep existing even while living in no 
#    areas.  Practically, this doesn't affect the downpass likelihood 
#    calculations, since under the DEC anagenesis model, it 
#    is impossible to evolve out of the null state; therefore
#    it is impossible to have "null" as an ancestor state if 
#    all the tips have positive ranges; therefore the likelihood
#    of the data for state "null" is always 0, all the way 
#    down to the root during downpass.
# 
#######################################################



#######################################################
# OUTLINE
# 
# 1. Calculation of tip state likelihoods under DEC, in 
#    BioGeoBEARS
#
# 2. Calculation of tree likelihood under a Yule process
# 
# 3. Set up of equivalent ClaSSE model
#
# 4. ClaSSE likelihoods and comparison to DEC
# 
#######################################################

# Set your working directory to wherever you unzipped this directory
wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/DEC_DEC+J_vs_ClaSSE_single_case/"
setwd(wd)



#######################################################
#######################################################
# 1. Calculation of tip state likelihoods under DEC, in 
#    BioGeoBEARS
#######################################################
#######################################################
library(diversitree)

# After installation, load the package, dependencies, updates
library(optimx)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)
source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)


# ClaSSE helper functions
source("/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/DEC_DEC+J_vs_ClaSSE_single_case/ClaSSE_mods_v1.R")



# BioGeoBEARS contains various example files and scripts in its extdata directory.
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))


# Look at your phylogeny:
trfn = "tree.newick"
tr = read.tree(trfn)
tr
plot(tr)
title("Example phylogeny")
axisPhylo() # plots timescale

# Look your geography data at the tips:
geogfn = "geog.data"
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 2







#######################################################
# Single example:
#######################################################
d_val = 0.05
e_val = 0.03
j_val = 0.1

#######################################################
# Run DEC
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size
# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Set the initial values for the "d" and "e" parameters
dstart = d_val
estart = e_val
jstart = j_val

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 5

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 5


# Do the likelihood calculation under these set parameters
# (this is just 1 likelihood calculation, not an ML search)
BGB_DEC_init_LnL = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE)
BGB_DEC_init_LnL
# -2.075039

#######################################################
#######################################################
# 2. Calculation of tree likelihood under a Yule process
#######################################################
#######################################################
library(ape)
trfn = "tree.newick"
tr = read.tree(trfn)

# Get the ML estimates of birthRate (speciation rate) 
# and deathRate (extinction rate)
# ...using ape::birthdeath
BD =  birthdeath(tr)
BD
names(BD)

# Calculate the birthRate and deathRate from the outputs
x1 = unname(BD$para["d/b"])
x2 = unname(BD$para["b-d"])
deathRate = (x2*x1) / (1-x1)
birthRate = deathRate+x2
c(birthRate, deathRate)
#birthRate = 0.22222222222222222222222222222222222222222222
# [1] 0.2222222 0.0000000

# You should get:
# c(birthRate, deathRate)
# [1] 0.2222222 0.0000000

# Get the log-likelihood of the tree under the ML parameters
# Convert the deviance to likelihood
BD_LnL = -1 * BD$dev / 2
BD_LnL
# -3.216395


# Double-check with yule() function:
yule_inf = yule(phy=tr)
yule_inf
# $lambda
# [1] 0.2222222
# 
# $se
# [1] 0.1571348
# 
# $loglik
# [1] -3.216395
# 
# attr(,"class")
# [1] "yule"

#######################################################
#######################################################
# 3. Set up of equivalent ClaSSE model for DEC, starting values
#######################################################
#######################################################



##################################################
# Set up states for ClaSSE to match DEC states
##################################################

# States in the ClaSSE model, and how they correspond
# to the geographic ranges in DEC
# 
# ====================
# (statenum = range)
# ====================
# 1 = null range
# 2 = A = Africa
# 3 = B = Asia
# 4 = AB = both
# ====================

# States at the tips of the tree
# (ranges A, A, A, B)
states = c(2, 2, 2, 3)
names(states) = tr$tip.label
states


# Proportion of species in each state
# (Let's assume we have all species)
sampling.f = c(1,1,1,1)

# Number of states
k = 4

# Make the ClaSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
classe_2areas = make.classe(tree=tr, states=states, k=4, sampling.f=sampling.f, strict=FALSE)
# Look at all the parameters of this model!
# lambdas = speciation rates
# mus = extinction rates
# qs = anagenetic transition rates
classe_2areas



##################################################
# Set up parameters for ClaSSE to match DEC
##################################################
# The names of the parameters:
param_names = argnames(classe_2areas)

# Most parameters will be zero
classe_params = rep(0, times=length(param_names))
names(classe_params) = param_names

# All extinction rates are the same (state-independent)
# Here, deathRate is 0 for all states
classe_params[grepl(pattern="mu", x=param_names)] = deathRate

# Transition rates are mostly 0, except as specified
# by DEC anagenesis model
d = d_val
e = e_val
classe_params[grepl(pattern="q", x=param_names)] = 0
classe_params[param_names == "q21"] = e
classe_params[param_names == "q31"] = e
classe_params[param_names == "q24"] = d
classe_params[param_names == "q34"] = d
classe_params[param_names == "q42"] = e
classe_params[param_names == "q43"] = e
classe_params



# The birthRate (lambda) is state-independent.  However, 
# only certain cladogenesis scenarios are allowed under DEC.
#
# Disallowed cladogenesis scenarios have a rate of 0.
#
# If there is more than one cladogenesis scenario conditional 
# on a certain ancestor, DEC assigns each a weight of 1, and 
# then divides by the sum of the weights. I.e., if there are
# six possible cladogenetic range-inheritance events, they 
# each get a conditional probability of 1/6.
# 
# To translate to ClaSSE, if the speciation rate for a lineage 
# in a certain state is lambda, then the rate of each individual 
# allowed scenario would be lambda * 1/6
# 
y_val = (3-j_val)/3
total_of_weights = y_val + j_val + j_val
yprob = y_val / total_of_weights
jprob = j_val / total_of_weights


# Specifying the nonzero lambdas
# Null range cannot speciate (doesn't seem to matter anyway,
# as "null" cannot be an ancestor anyway)
classe_params[param_names=="lambda111"] = birthRate

# Narrow sympatry (ancestor A or B; rangesize of 1 area)
classe_params[param_names=="lambda222"] = yprob * birthRate
classe_params[param_names=="lambda333"] = yprob * birthRate

# Jump dispersal speciation
classe_params[param_names=="lambda223"] = jprob * birthRate
classe_params[param_names=="lambda323"] = jprob * birthRate

# Subset sympatry for state AB
classe_params[param_names=="lambda424"] = 1/6 * birthRate
classe_params[param_names=="lambda434"] = 1/6 * birthRate

# Vicariance for state AB
classe_params[param_names=="lambda423"] = 1/6 * birthRate

classe_params_DEC = classe_params


# Note: Under ClaSSE, these rates are implicitly assigned
# to the branch rotations:
#
# lambda424 is automatically also assigned to lambda442
# lambda434 is automatically also assigned to lambda443
# lambda423 is automatically also assigned to lambda432
# 
# Note that this symmetry assumption would exclude modeling
# of directional bias in founder-event speciation events

# Calculate the ClaSSE likelihood under these parameters
# Note that decision about the likelihood of the root states effects the LnL
res1 = classe_2areas(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2 = classe_2areas(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.1,0.2,0.3,0.4)
res3 = classe_2areas(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.0,1/3,1/3,1/3)
res4 = classe_2areas(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.25,0.25,0.25,0.25)
res5 = classe_2areas(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)

LnLs1 = get_classe_LnLs(res1)
LnLs2 = get_classe_LnLs(res2)
LnLs3 = get_classe_LnLs(res3)
LnLs4 = get_classe_LnLs(res4)
LnLs5 = get_classe_LnLs(res5)

LnLs = rbind(LnLs1, LnLs2, LnLs3, LnLs4, LnLs5)
print(LnLs)

# You should get:
# d=0.05, e=0.03, j=0.1
# LnLs1 -8.202212 -5.764771
# LnLs2 -8.785001 -5.764771
# LnLs3 -8.428936 -5.764771
# LnLs4 -8.427713 -5.764771
# LnLs5 -8.785001 -5.764771
# 
# The remaining bit of likelihood is:
# -8.785001--5.764771
# [1] -3.02023

# sum(attr(res5, "intermediates")$lq)
# -5.764771

# log(sum(attr(res5, "intermediates")$vals))
# -3.321764


# Calculate the likelihood under these parameters
res1 = classe_2areas(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = classe_2areas(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.1,0.2,0.3,0.4)
res3 = classe_2areas(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.0,1/3,1/3,1/3)
res4 = classe_2areas(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.25,0.25,0.25,0.25)
res5 = classe_2areas(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)

LnLs1 = get_classe_LnLs(res1)
LnLs2 = get_classe_LnLs(res2)
LnLs3 = get_classe_LnLs(res3)
LnLs4 = get_classe_LnLs(res4)
LnLs5 = get_classe_LnLs(res5)

LnLs = rbind(LnLs1, LnLs2, LnLs3, LnLs4, LnLs5)
print(LnLs)

# You should get:
# d=0.05, e=0.03, j=0.1
# LnLs1 -10.05694 -5.764771
# LnLs2 -10.47283 -5.764771
# LnLs3 -10.21122 -5.764771
# LnLs4 -10.18515 -5.764771
# LnLs5 -10.47283 -5.764771

classe_DEC_ML = classe_2areas(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
classe_DEC_ML

sum(attr(classe_DEC_ML, "intermediates")$lq)
# d=0.05, e=0.03, j=0.1
# -5.764771

classe_DEC_ML_minus_BD = sum(attr(classe_DEC_ML, "intermediates")$lq) - BD_LnL
# -2.548376

BGB_DEC_init_LnL
# -2.075039

# (Not an exact match, numerical calculation with multiple processes 
#  - look at correlation plot via a for-loop)


# How the root state probabilities contribute to the likelihood:
# d.root is the probability densities at the root(
# vals = the $init at the root, i.e. the 5th column
# init[1:4,] = the conditional probabilities of extinction given each state
# init[5:8,] = the conditional probabilities of the tree+data 
#              above given each state
# d.root = init[5:8,][,5] = D1t, D2t, D3t, D4t,
vals = attr(res5,"vals")
vals
k <- length(vals)/2
i <- seq_len(k)
d.root <- vals[-i]
d.root
# 0.000000000 0.009364228 0.008426473 0.018298426

root.p=attr(res5,"intermediates")$root.p
root.p
# 0.25 0.25 0.25 0.25

log(sum(root.p * d.root))
-4.708058

# The ClaSSE returned likelihood is the sum of these two:
-4.708058 + -5.764771
-10.47283



