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
wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/4species_ABalternate_ML_DECvDECj_BGB/"
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
# This is the example Newick file for Hawaiian ABalternate
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = "tree_4species.newick"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
tr
# plot(tr)
# title("Example 4-species phylogeny, with branchlengths of 1")
# axisPhylo() # plots timescale

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

# This is the example geography file for Hawaiian ABalternate
# (from Ree & Smith 2008)
geogfn = "geog_4species.data"

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
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
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
resfn = "ABalternate_DEC_M0_unconstrained_v1.Rdata"
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
# [1] -5.196338
# > exp(resDEC$total_loglikelihood)
# [1] 0.005536806
# > resDEC$outputs@params_table["d","est"]
# [1] 0.06541686
# > resDEC$outputs@params_table["e","est"]
# [1] 1e-12
# > resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node
#      [,1]         [,2]         [,3]      [,4]
# [1,]    0 0.000000e+00 1.000000e+00 0.0000000
# [2,]    0 1.000000e+00 0.000000e+00 0.0000000
# [3,]    0 0.000000e+00 1.000000e+00 0.0000000
# [4,]    0 1.000000e+00 0.000000e+00 0.0000000
# [5,]    0 4.225776e-13 2.136949e-01 0.7863051
# [6,]    0 1.882430e-01 1.544321e-13 0.8117570
# [7,]    0 5.585723e-13 1.434443e-13 1.0000000


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
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
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

resfn = "ABalternate_DEC+J_M0_unconstrained_v1.Rdata"
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
# [1] -1.170586
# > exp(resDECj$total_loglikelihood)
# [1] 0.3101852
# > resDECj$outputs@params_table["d","est"]
# [1] 1e-12
# > resDECj$outputs@params_table["e","est"]
# [1] 1e-12
# > resDECj$outputs@params_table["j","est"]
# [1] 2.999965
# > resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
#      [,1]         [,2]         [,3]       [,4]
# [1,]    0 0.000000e+00 1.000000e+00 0.00000000
# [2,]    0 1.000000e+00 0.000000e+00 0.00000000
# [3,]    0 0.000000e+00 1.000000e+00 0.00000000
# [4,]    0 1.000000e+00 0.000000e+00 0.00000000
# [5,]    0 4.029843e-01 4.029859e-01 0.19402978
# [6,]    0 9.402970e-01 1.554478e-06 0.05970141
# [7,]    0 3.627115e-06 9.850710e-01 0.01492537





#######################################################
# PDF plots
#######################################################
pdffn = "ABalternate_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on ABalternate M0_unconstrained"

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
analysis_titletxt ="BioGeoBEARS DEC+J on ABalternate M0_unconstrained"

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







names(resDECj)
resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
resDECj$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
resDECj$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS

resDECj$computed_likelihoods_at_each_node
resDECj$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
resDECj$condlikes_of_each_state
resDECj$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
resDECj$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
resDECj$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
resDECj$ML_marginal_prob_each_state_at_branch_bottom_below_node
resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
resDECj$relative_probs_of_each_state_at_bottom_of_root_branch
resDECj$total_loglikelihood
resDECj$optim_result



#######################################################
# C++ Lagrange:
#######################################################
# 
# Optimizing (simplex) -ln likelihood.
# dis: 0.849427 ext: 7.1275e-11
# final -ln likelihood: 1.28778
# Ancestral splits for:	1
# 	A_B|A	0.404894	(2.19191)
# 	A|A	0.363369	(2.30012)
# 	B|A	0.231737	(2.74993)
# 
# Ancestral states for:	1
# 	A_B	0.636631	(1.73935)
# 	A	0.363369	(2.30012)
# 
# Ancestral splits for:	2
# 	A|A	0.265599	(2.61355)
# 	A_B|A	0.176341	(3.02311)
# 	B|A	0.14409	(3.2251)
# 	B|B	0.100914	(3.58127)
# 	B|A_B	0.0956541	(3.6348)
# 	A_B|B	0.0670006	(3.99083)
# 	A|B	0.0547467	(4.19282)
# 
# Ancestral states for:	2
# 	A_B	0.633487	(1.7443)
# 	A	0.265599	(2.61355)
# 	B	0.100914	(3.58127)
# 
# Ancestral splits for:	3
# 	A|A	0.390401	(2.22836)
# 	B|B	0.268183	(2.60386)
# 	A_B|A	0.0705879	(3.93868)
# 	B|A	0.0650668	(4.02012)
# 	B|A_B	0.0562871	(4.16507)
# 
# Ancestral states for:	3
# 	A	0.390401	(2.22836)
# 	A_B	0.341416	(2.36243)
# 	B	0.268183	(2.60386)


#######################################################
# Python Lagrange:
#######################################################
# Newick tree with interior nodes labeled:
# (sp4:3.0,(sp3:2.0,(sp2:1.0,sp1:1.0)N4:1.0)N5:1.0)N6:0.0;
# 
# 
# Cladogram (branch lengths not to scale):
#   ------------------------------------------------------------------+ [AB] sp4 
# N6+                                                                            
#   :                     --------------------------------------------+ [AB] sp3 
#   --------------------N5+                                                      
#                         :                     ----------------------+ [AB] sp2 
#                         --------------------N4+                                
#                                               ----------------------+ [A] sp1  
# 
# 
# 
# Optimization terminated successfully.
#          Current function value: 1.287780
#          Iterations: 3
#          Function evaluations: 140
# Global ML at root node:
#   -lnL = 1.288
#   dispersal = 0.8493
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
# At node N6:
#    split   lnL     Rel.Prob
#    [A|A]   -2.228  0.3904  
#    [B|B]   -2.604  0.2682  
#    [AB|A]  -3.939  0.07059 
#    [B|A]   -4.02   0.06507 
#    [B|AB]  -4.165  0.05629 
#    [A|AB]  -4.165  0.05629 
#    [AB|B]  -4.314  0.04849 
# 
# At node N5:
#    split   lnL     Rel.Prob
#    [A|A]   -2.613  0.2656  
#    [AB|A]  -3.023  0.1763  
#    [B|A]   -3.225  0.1441  
#    [B|B]   -3.581  0.1009  
#    [B|AB]  -3.635  0.09565 
#    [A|AB]  -3.635  0.09565 
#    [AB|B]  -3.991  0.067   
#    [A|B]   -4.193  0.05474 
# 
# At node N4:
#    split   lnL     Rel.Prob
#    [AB|A]  -2.192  0.4049  
#    [A|A]   -2.3    0.3634  
