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

########################################################
# TO GET THE OPTIMX/OPTIM FIX, AND THE UPPASS FIX, 
# SOURCE THE REVISED FUNCTIONS WITH THESE COMMANDS
#
# CRUCIAL CRUCIAL CRUCIAL: 
# YOU HAVE TO RUN THE SOURCE COMMANDS AFTER 
# *EVERY TIME* YOU DO library(BioGeoBEARS). THE CHANGES ARE NOT "PERMANENT", 
# THEY HAVE TO BE MADE EACH TIME.  IF YOU ARE GOING TO BE OFFLINE, 
# YOU CAN DOWNLOAD EACH .R FILE TO YOUR HARD DRIVE AND REFER THE source()
# COMMANDS TO THE FULL PATH AND FILENAME OF EACH FILE ON YOUR
# LOCAL SYSTEM INSTEAD.
########################################################
# library(BioGeoBEARS)
# source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
# calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
# calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
    # slight speedup hopefully

#######################################################
# Local source()-ing method -- uses BioGeoBEARS sourceall() function 
# on a directory of .R files, so you don't have to type them out.
# The directories here are on my machine, you would have to make a 
# directory, save the .R files there, and refer to them.
#
# NOTE: it's best to source the "cladoRcpp.R" update first, to avoid warnings like this:
##
## Note: possible error in 'rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs = tmpca_1, ': 
##         unused arguments (m = m, m_null_range = include_null_range, jts_matrix = jts_matrix) 
##
#
# TO USE: Delete or comment out the 'source("http://...")' commands above, and un-comment
#              the below...
########################################################################
# Un-comment (and fix directory paths) to use:
library(BioGeoBEARS)
source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
########################################################################

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

# Set your working directory for output files
# default here is your home directory ("~")
# Change this as you like
wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/9species_AvsB_ML_DECvDECj_BGB/"
setwd(wd)

# Double-check your working directory with getwd()
getwd()

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Example files are given below. To run your own data,
# make the below lines point to your own files, e.g.
# trfn = "/mydata/frogs/frogBGB/tree.newick"
# geogfn = "/mydata/frogs/frogBGB/geog.data"

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################
# This is the example Newick file for Hawaiian AvsB
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = "tree_9species.newick"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
tr
plot(tr)
title("Example 9-species phylogeny, with branchlengths of 1")
axisPhylo() # plots timescale

#######################################################
# Geography file
# Notes:
# 1. This is a PHYLIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
#    - after a tab, put the areas in parentheses, with spaces: (A B C D)
#
# 1.5. Example first line:
#    10    4    (A B C D)
# 
# 2. The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
# 
# 2.5a. This means a TAB between the species name and the area 0/1s
# 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
# 
# 3. See example files at:
#    http://phylo.wikidot.com/biogeobears#files
# 
# 4. Make you understand what a PLAIN-TEXT EDITOR is:
#    http://phylo.wikidot.com/biogeobears#texteditors
#
# 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
#
# 4. All names in the geography file must match names in the phylogeny file.
#
# 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#
# 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################

# This is the example geography file for Hawaiian AvsB
# (from Ree & Smith 2008)
geogfn = "geog_9species.data"

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 2

####################################################
####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
numstates_from_numareas(numareas=2, maxareas=2, include_null_range=TRUE)
numstates_from_numareas(numareas=2, maxareas=2, include_null_range=FALSE)
#######################################################


#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DEC" model is identical with 
# the Lagrange DEC model, and should return identical
# ML estimates of parameters, and the same 
# log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly 
# different, since BioGeoBEARS is reporting the 
# ancestral state probabilities under the global ML
# model, and Lagrange is reporting ancestral state
# probabilities after re-optimizing the likelihood
# after fixing the state at each node. These will 
# be similar, but not identical. See Matzke (2014),
# Systematic Biology, for discussion.
#
# Also see Matzke (2014) for presentation of the 
# DEC+J model.
#######################################################
#######################################################

#######################################################
#######################################################

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "AvsB_DEC_M0_unconstrained_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    

    save(res, file=resfn)
    resDEC = res
    } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
    }

resDEC$total_loglikelihood
exp(resDEC$total_loglikelihood)
resDEC$outputs@params_table["d","est"]
resDEC$outputs@params_table["e","est"]
resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node

# > resDEC$total_loglikelihood
# [1] -1.315831
# > exp(resDEC$total_loglikelihood)
# [1] 0.2682512
# > resDEC$outputs@params_table["d","est"]
# [1] 0.9040951
# > resDEC$outputs@params_table["e","est"]
# [1] 1e-12
# > resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node
#       [,1]      [,2]         [,3]      [,4]
#  [1,]    0 0.0000000 0.000000e+00 1.0000000
#  [2,]    0 0.0000000 0.000000e+00 1.0000000
#  [3,]    0 0.0000000 0.000000e+00 1.0000000
#  [4,]    0 0.0000000 0.000000e+00 1.0000000
#  [5,]    0 0.0000000 0.000000e+00 1.0000000
#  [6,]    0 0.0000000 0.000000e+00 1.0000000
#  [7,]    0 0.0000000 0.000000e+00 1.0000000
#  [8,]    0 0.0000000 0.000000e+00 1.0000000
#  [9,]    0 1.0000000 0.000000e+00 0.0000000
# [10,]    0 0.3338911 3.326952e-01 0.3334137
# [11,]    0 0.1806555 1.790608e-01 0.6402838
# [12,]    0 0.1607794 1.572923e-01 0.6819283
# [13,]    0 0.1601658 1.516948e-01 0.6881394
# [14,]    0 0.1652699 1.443367e-01 0.6903934
# [15,]    0 0.1791951 1.271244e-01 0.6936805
# [16,]    0 0.2171531 8.606677e-02 0.6967801
# [17,]    0 0.3410920 2.133868e-13 0.6589080


#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)






resfn = "AvsB_DEC+J_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res = bears_optim_run(BioGeoBEARS_run_object)
    res    

    save(res, file=resfn)

    resDECj = res
    } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
    }


resDECj$total_loglikelihood
exp(resDECj$total_loglikelihood)
resDECj$outputs@params_table["d","est"]
resDECj$outputs@params_table["e","est"]
resDECj$outputs@params_table["j","est"]
resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node

# > resDECj$total_loglikelihood
# [1] -1.315823
# > exp(resDECj$total_loglikelihood)
# [1] 0.2682535
# > resDECj$outputs@params_table["d","est"]
# [1] 0.9064843
# > resDECj$outputs@params_table["e","est"]
# [1] 1e-12
# > resDECj$outputs@params_table["j","est"]
# [1] 0.1808252
# > resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
#       [,1]      [,2]       [,3]      [,4]
#  [1,]    0 0.0000000 0.00000000 1.0000000
#  [2,]    0 0.0000000 0.00000000 1.0000000
#  [3,]    0 0.0000000 0.00000000 1.0000000
#  [4,]    0 0.0000000 0.00000000 1.0000000
#  [5,]    0 0.0000000 0.00000000 1.0000000
#  [6,]    0 0.0000000 0.00000000 1.0000000
#  [7,]    0 0.0000000 0.00000000 1.0000000
#  [8,]    0 0.0000000 0.00000000 1.0000000
#  [9,]    0 1.0000000 0.00000000 0.0000000
# [10,]    0 0.3333374 0.33325042 0.3334121
# [11,]    0 0.1795108 0.17935010 0.6411391
# [12,]    0 0.1588453 0.15835765 0.6827970
# [13,]    0 0.1563275 0.15468285 0.6889897
# [14,]    0 0.1572081 0.15156653 0.6912254
# [15,]    0 0.1624974 0.14301824 0.6944844
# [16,]    0 0.1852622 0.11719713 0.6975407
# [17,]    0 0.2931028 0.04729853 0.6595987



#######################################################
# PDF plots
#######################################################
pdffn = "AvsB_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on AvsB M0_unconstrained"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on AvsB M0_unconstrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it






#######################################################
# C++ Lagrange:
#######################################################
# 
# Optimizing (simplex) -ln likelihood.
# dis: 0.906435 ext: 3.52304e-08
# final -ln likelihood: 1.31582
# Ancestral splits for:	1
# 	A_B|A	0.413263	(2.19949)
# 	A|A	0.340416	(2.39341)
# 	B|A	0.246321	(2.71694)
# 
# Ancestral states for:	1
# 	A_B	0.659584	(1.73197)
# 	A	0.340416	(2.39341)
# 
# Ancestral splits for:	2
# 	A|A	0.216511	(2.84594)
# 	A_B|A	0.189495	(2.97922)
# 	B|A	0.158572	(3.15737)
# 	B|A_B	0.105631	(3.56363)
# 	B|B	0.0859642	(3.76965)
# 	A_B|B	0.0752375	(3.90293)
# 	A|B	0.0629599	(4.08108)
# 
# Ancestral states for:	2
# 	A_B	0.697525	(1.67604)
# 	A	0.216511	(2.84594)
# 
# Ancestral splits for:	3
# 	A|A	0.178637	(3.03823)
# 	A_B|A	0.140676	(3.27712)
# 	B|A	0.131403	(3.34531)
# 	B|B	0.126896	(3.38021)
# 	B|A_B	0.114558	(3.4825)
# 	A_B|B	0.0999304	(3.61911)
# 	A|B	0.0933429	(3.6873)
# 
# Ancestral states for:	3
# 	A_B	0.694468	(1.68043)
# 	A	0.178637	(3.03823)
# 	B	0.126896	(3.38021)
# 
# Ancestral splits for:	4
# 	A|A	0.164773	(3.11901)
# 	B|B	0.144019	(3.25363)
# 	A_B|A	0.124675	(3.39787)
# 	B|A	0.121355	(3.42486)
# 	B|A_B	0.115067	(3.47806)
# 	A_B|B	0.108973	(3.53248)
# 	A|B	0.106071	(3.55947)
# 
# Ancestral states for:	4
# 	A_B	0.691208	(1.68514)
# 	A	0.164773	(3.11901)
# 	B	0.144019	(3.25363)
# 
# Ancestral splits for:	5
# 	A|A	0.159703	(3.15026)
# 	B|B	0.151325	(3.20415)
# 	A_B|A	0.118567	(3.4481)
# 	B|A	0.117292	(3.45891)
# 	B|A_B	0.114814	(3.48026)
# 	A_B|B	0.112347	(3.50199)
# 	A|B	0.111138	(3.51281)
# 
# Ancestral states for:	5
# 	A_B	0.688972	(1.68838)
# 	A	0.159703	(3.15026)
# 	B	0.151325	(3.20415)
# 
# Ancestral splits for:	6
# 	A|A	0.160331	(3.14634)
# 	B|B	0.15689	(3.16803)
# 	A_B|A	0.115282	(3.4762)
# 	B|A	0.114781	(3.48055)
# 	B|A_B	0.113795	(3.48918)
# 	A_B|B	0.112808	(3.49789)
# 	A|B	0.112318	(3.50225)
# 
# Ancestral states for:	6
# 	A_B	0.682779	(1.69741)
# 	A	0.160331	(3.14634)
# 	B	0.15689	(3.16803)
# 
# Ancestral splits for:	7
# 	A|A	0.180225	(3.02937)
# 	B|B	0.178654	(3.03813)
# 	A_B|A	0.107416	(3.54687)
# 	B|A	0.107227	(3.54863)
# 	B|A_B	0.106853	(3.55212)
# 	A_B|B	0.106479	(3.55563)
# 	A|B	0.106293	(3.55738)
# 
# Ancestral states for:	7
# 	A_B	0.641121	(1.76036)
# 	A	0.180225	(3.02937)
# 	B	0.178654	(3.03813)
# 
# Ancestral splits for:	8
# 	A|A	0.333883	(2.41279)
# 	B|B	0.332705	(2.41632)
# 	A_B|A	0.0556866	(4.20384)
# 	B|A	0.0556471	(4.20455)
# 	B|A_B	0.0555687	(4.20596)
# 	A_B|B	0.0554902	(4.20737)
# 	A|B	0.0554508	(4.20808)
# 
# Ancestral states for:	8
# 	A	0.333883	(2.41279)
# 	A_B	0.333412	(2.4142)
# 	B	0.332705	(2.41632)
# 





#######################################################
# Python Lagrange:
#######################################################
# Newick tree with interior nodes labeled:
# (sp9:8.0,(sp8:7.0,(sp7:6.0,(sp6:5.0,(sp5:4.0,(sp4:3.0,(sp3:2.0,(sp2:1.0,sp1:1.0)N9:1.0)N10:1.0)N11:1.0)N12:1.0)N13:1.0)N14:1.0)N15:1.0)N16:0.0;
# 
# 
# Cladogram (branch lengths not to scale):
#    ----------------------------------------------------------------+ [AB] sp9 
# N16+                                                                          
#    :       --------------------------------------------------------+ [AB] sp8 
#    -----N15+                                                                  
#            :       ------------------------------------------------+ [AB] sp7 
#            -----N14+                                                          
#                    :       ----------------------------------------+ [AB] sp6 
#                    -----N13+                                                  
#                            :       --------------------------------+ [AB] sp5 
#                            -----N12+                                          
#                                    :       ------------------------+ [AB] sp4 
#                                    -----N11+                                  
#                                            :       ----------------+ [AB] sp3 
#                                            -----N10+                          
#                                                    :       --------+ [AB] sp2 
#                                                    ------N9+                  
#                                                            --------+ [A] sp1  
# 
# 
# 
# Optimization terminated successfully.
#          Current function value: 1.315823
#          Iterations: 3
#          Function evaluations: 140
# Global ML at root node:
#   -lnL = 1.316
#   dispersal = 0.9065
#   extinction = 4.285e-09
# 
# Ancestral range subdivision/inheritance scenarios ('splits') at
# internal nodes.
# 
# * Split format: [left|right], where 'left' and 'right' are the ranges
#   inherited by each descendant branch (on the printed tree, 'left' is
#   the upper branch, and 'right' the lower branch).
# 
# * Only splits within 2 log-likelihood units of the maximum for each
#   node are shown.  'Rel.Prob' is the relative probability (fraction of
#   the global likelihood) of a split.
# 
# At node N16:
#    split   lnL     Rel.Prob
#    [A|A]   -2.413  0.3339  
#    [B|B]   -2.416  0.3327  
#    [AB|A]  -4.204  0.05569 
#    [B|A]   -4.205  0.05565 
#    [B|AB]  -4.206  0.05557 
#    [A|AB]  -4.206  0.05557 
#    [AB|B]  -4.207  0.05549 
#    [A|B]   -4.208  0.05545 
# 
# At node N15:
#    split   lnL     Rel.Prob
#    [A|A]   -3.029  0.1802  
#    [B|B]   -3.038  0.1786  
#    [AB|A]  -3.547  0.1074  
#    [B|A]   -3.549  0.1072  
#    [B|AB]  -3.552  0.1069  
#    [A|AB]  -3.552  0.1069  
#    [AB|B]  -3.556  0.1065  
#    [A|B]   -3.557  0.1063  
# 
# At node N14:
#    split   lnL     Rel.Prob
#    [A|A]   -3.146  0.1603  
#    [B|B]   -3.168  0.1569  
#    [AB|A]  -3.476  0.1153  
#    [B|A]   -3.481  0.1148  
#    [B|AB]  -3.489  0.1138  
#    [A|AB]  -3.489  0.1138  
#    [AB|B]  -3.498  0.1128  
#    [A|B]   -3.502  0.1123  
# 
# At node N13:
#    split   lnL     Rel.Prob
#    [A|A]   -3.15   0.1597  
#    [B|B]   -3.204  0.1513  
#    [AB|A]  -3.448  0.1186  
#    [B|A]   -3.459  0.1173  
#    [A|AB]  -3.48   0.1148  
#    [B|AB]  -3.48   0.1148  
#    [AB|B]  -3.502  0.1123  
#    [A|B]   -3.513  0.1111  
# 
# At node N12:
#    split   lnL     Rel.Prob
#    [A|A]   -3.119  0.1648  
#    [B|B]   -3.254  0.144   
#    [AB|A]  -3.398  0.1247  
#    [B|A]   -3.425  0.1214  
#    [B|AB]  -3.478  0.1151  
#    [A|AB]  -3.478  0.1151  
#    [AB|B]  -3.532  0.109   
#    [A|B]   -3.559  0.1061  
# 
# At node N11:
#    split   lnL     Rel.Prob
#    [A|A]   -3.038  0.1786  
#    [AB|A]  -3.277  0.1407  
#    [B|A]   -3.345  0.1314  
#    [B|B]   -3.38   0.1269  
#    [B|AB]  -3.482  0.1146  
#    [A|AB]  -3.482  0.1146  
#    [AB|B]  -3.619  0.09993 
#    [A|B]   -3.687  0.09335 
# 
# At node N10:
#    split   lnL     Rel.Prob
#    [A|A]   -2.846  0.2165  
#    [AB|A]  -2.979  0.1895  
#    [B|A]   -3.157  0.1586  
#    [B|AB]  -3.564  0.1056  
#    [A|AB]  -3.564  0.1056  
#    [B|B]   -3.77   0.08596 
#    [AB|B]  -3.903  0.07524 
#    [A|B]   -4.081  0.06296 
# 
# At node N9:
#    split   lnL     Rel.Prob
#    [AB|A]  -2.199  0.4133  
#    [A|A]   -2.393  0.3404  
#    [B|A]   -2.717  0.2463  

