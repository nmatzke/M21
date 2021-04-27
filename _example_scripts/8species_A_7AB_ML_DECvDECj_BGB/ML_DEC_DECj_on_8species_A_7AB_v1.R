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
wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/8species_A_7AB_ML_DECvDECj_BGB/"
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
# This is the example Newick file for Hawaiian A_7AB
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = "tree_8species.newick"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
# tr
# plot(tr)
# title("Example 8-species phylogeny, with branchlengths of 1")
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

# This is the example geography file for Hawaiian A_7AB
# (from Ree & Smith 2008)
geogfn = "geog_8species.data"

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
resfn = "A_7AB_DEC_M0_unconstrained_v1.Rdata"
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
# [1] -1.315462
# > exp(resDEC$total_loglikelihood)
# [1] 0.2683504
# > resDEC$outputs@params_table["d","est"]
# [1] 0.9045209
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
#  [8,]    0 1.0000000 0.000000e+00 0.0000000
#  [9,]    0 0.3347082 3.317606e-01 0.3335312
# [10,]    0 0.1815994 1.776683e-01 0.6407324
# [11,]    0 0.1629373 1.543323e-01 0.6827304
# [12,]    0 0.1655765 1.446259e-01 0.6897976
# [13,]    0 0.1791509 1.271237e-01 0.6937254
# [14,]    0 0.2170452 8.605174e-02 0.6969030
# [15,]    0 0.3409705 2.134464e-13 0.6590295

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






resfn = "A_7AB_DEC+J_M0_unconstrained_v1.Rdata"
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
# [1] -1.31546
# > exp(resDECj$total_loglikelihood)
# [1] 0.2683507
# > resDECj$outputs@params_table["d","est"]
# [1] 0.9054853
# > resDECj$outputs@params_table["e","est"]
# [1] 1e-12
# > resDECj$outputs@params_table["j","est"]
# [1] 0.3595307
# > resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
#       [,1]      [,2]       [,3]      [,4]
#  [1,]    0 0.0000000 0.00000000 1.0000000
#  [2,]    0 0.0000000 0.00000000 1.0000000
#  [3,]    0 0.0000000 0.00000000 1.0000000
#  [4,]    0 0.0000000 0.00000000 1.0000000
#  [5,]    0 0.0000000 0.00000000 1.0000000
#  [6,]    0 0.0000000 0.00000000 1.0000000
#  [7,]    0 0.0000000 0.00000000 1.0000000
#  [8,]    0 1.0000000 0.00000000 0.0000000
#  [9,]    0 0.3332575 0.33321263 0.3335299
# [10,]    0 0.1795163 0.17940765 0.6410760
# [11,]    0 0.1586772 0.15824467 0.6830781
# [12,]    0 0.1558903 0.15397490 0.6901348
# [13,]    0 0.1573002 0.14864956 0.6940502
# [14,]    0 0.1712004 0.13158943 0.6972101
# [15,]    0 0.2640985 0.07659314 0.6593083




#######################################################
# PDF plots
#######################################################
pdffn = "A_7AB_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on A_7AB M0_unconstrained"

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
analysis_titletxt ="BioGeoBEARS DEC+J on A_7AB M0_unconstrained"

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
# Optimizing (simplex) -ln likelihood.
# dis: 0.905516 ext: 4.44398e-10
# final -ln likelihood: 1.31546
# Ancestral splits for:	1
# 	A_B|A	0.413192	(2.1993)
# 	A|A	0.340683	(2.39226)
# 	B|A	0.246125	(2.71738)
# 
# Ancestral states for:	1
# 	A_B	0.659317	(1.73201)
# 	A	0.340683	(2.39226)
# 
# Ancestral splits for:	2
# 	A|A	0.216772	(2.84437)
# 	A_B|A	0.189484	(2.97891)
# 	B|A	0.158506	(3.15742)
# 	B|A_B	0.105579	(3.56375)
# 	B|B	0.0860081	(3.76877)
# 	A_B|B	0.0751811	(3.90332)
# 	A|B	0.0628901	(4.08183)
# 
# Ancestral states for:	2
# 	A_B	0.69722	(1.67611)
# 	A	0.216772	(2.84437)
# 
# Ancestral splits for:	3
# 	A|A	0.178913	(3.03632)
# 	A_B|A	0.140638	(3.27703)
# 	B|A	0.131341	(3.34542)
# 	B|B	0.127026	(3.37882)
# 	B|A_B	0.114489	(3.48273)
# 	A_B|B	0.0998514	(3.61953)
# 	A|B	0.093251	(3.68792)
# 
# Ancestral states for:	3
# 	A_B	0.694061	(1.68066)
# 	A	0.178913	(3.03632)
# 	B	0.127026	(3.37882)
# 
# Ancestral splits for:	4
# 	A|A	0.165364	(3.11507)
# 	B|B	0.14449	(3.25)
# 	A_B|A	0.124509	(3.39884)
# 	B|A	0.121181	(3.42593)
# 	B|A_B	0.11489	(3.47924)
# 	A_B|B	0.108792	(3.53378)
# 	A|B	0.105884	(3.56087)
# 
# Ancestral states for:	4
# 	A_B	0.690146	(1.68631)
# 	A	0.165364	(3.11507)
# 	B	0.14449	(3.25)
# 
# Ancestral splits for:	5
# 	A|A	0.162738	(3.13108)
# 	B|B	0.154173	(3.18514)
# 	A_B|A	0.117568	(3.4562)
# 	B|A	0.116297	(3.46707)
# 	B|A_B	0.113834	(3.48848)
# 	A_B|B	0.11138	(3.51027)
# 	A|B	0.110176	(3.52113)
# 
# Ancestral states for:	5
# 	A_B	0.683089	(1.69659)
# 	A	0.162738	(3.13108)
# 	B	0.154173	(3.18514)
# 
# Ancestral splits for:	6
# 	A|A	0.181411	(3.02245)
# 	B|B	0.177502	(3.04424)
# 	A_B|A	0.108249	(3.53878)
# 	B|A	0.107776	(3.54316)
# 	B|A_B	0.106847	(3.55182)
# 	A_B|B	0.105916	(3.56057)
# 	A|B	0.105453	(3.56495)
# 
# Ancestral states for:	6
# 	A_B	0.641087	(1.76005)
# 	A	0.181411	(3.02245)
# 	B	0.177502	(3.04424)
# 
# Ancestral splits for:	7
# 	A|A	0.334701	(2.40998)
# 	B|B	0.331769	(2.41878)
# 	A_B|A	0.0558822	(4.19997)
# 	B|A	0.0557835	(4.20174)
# 	B|A_B	0.0555883	(4.20524)
# 	A_B|B	0.0553928	(4.20877)
# 	A|B	0.0552949	(4.21054)
# 
# Ancestral states for:	7
# 	A	0.334701	(2.40998)
# 	A_B	0.33353	(2.41348)
# 	B	0.331769	(2.41878)




#######################################################
# Python Lagrange:
#######################################################
# Newick tree with interior nodes labeled:
# (sp8:7.0,(sp7:6.0,(sp6:5.0,(sp5:4.0,(sp4:3.0,(sp3:2.0,(sp2:1.0,sp1:1.0)N8:1.0)N9:1.0)N10:1.0)N11:1.0)N12:1.0)N13:1.0)N14:0.0;
# 
# 
# Cladogram (branch lengths not to scale):
#    ---------------------------------------------------------------+ [AB] sp8 
# N14+                                                                         
#    :        ------------------------------------------------------+ [AB] sp7 
#    ------N13+                                                                
#             :        ---------------------------------------------+ [AB] sp6 
#             ------N12+                                                       
#                      :        ------------------------------------+ [AB] sp5 
#                      ------N11+                                              
#                               :        ---------------------------+ [AB] sp4 
#                               ------N10+                                     
#                                        :        ------------------+ [AB] sp3 
#                                        -------N9+                            
#                                                 :        ---------+ [AB] sp2 
#                                                 -------N8+                   
#                                                          ---------+ [A] sp1  
# 
# 
# 
# Optimization terminated successfully.
#          Current function value: 1.315461
#          Iterations: 3
#          Function evaluations: 140
# Global ML at root node:
#   -lnL = 1.315
#   dispersal = 0.9053
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
# At node N14:
#    split   lnL     Rel.Prob
#    [A|A]   -2.41   0.3347  
#    [B|B]   -2.419  0.3318  
#    [AB|A]  -4.2    0.05588 
#    [B|A]   -4.202  0.05578 
#    [B|AB]  -4.205  0.05559 
#    [A|AB]  -4.205  0.05559 
#    [AB|B]  -4.209  0.05539 
#    [A|B]   -4.211  0.05529 
# 
# At node N13:
#    split   lnL     Rel.Prob
#    [A|A]   -3.022  0.1814  
#    [B|B]   -3.044  0.1775  
#    [AB|A]  -3.539  0.1082  
#    [B|A]   -3.543  0.1078  
#    [B|AB]  -3.552  0.1068  
#    [A|AB]  -3.552  0.1068  
#    [AB|B]  -3.561  0.1059  
#    [A|B]   -3.565  0.1054  
# 
# At node N12:
#    split   lnL     Rel.Prob
#    [A|A]   -3.131  0.1628  
#    [B|B]   -3.185  0.1542  
#    [AB|A]  -3.456  0.1176  
#    [B|A]   -3.467  0.1163  
#    [B|AB]  -3.489  0.1138  
#    [A|AB]  -3.489  0.1138  
#    [AB|B]  -3.51   0.1114  
#    [A|B]   -3.521  0.1102  
# 
# At node N11:
#    split   lnL     Rel.Prob
#    [A|A]   -3.115  0.1654  
#    [B|B]   -3.25   0.1445  
#    [AB|A]  -3.399  0.1245  
#    [B|A]   -3.426  0.1212  
#    [B|AB]  -3.479  0.1149  
#    [A|AB]  -3.479  0.1149  
#    [AB|B]  -3.534  0.1088  
#    [A|B]   -3.561  0.1059  
# 
# At node N10:
#    split   lnL     Rel.Prob
#    [A|A]   -3.036  0.179   
#    [AB|A]  -3.277  0.1406  
#    [B|A]   -3.345  0.1313  
#    [B|B]   -3.379  0.127   
#    [B|AB]  -3.483  0.1145  
#    [A|AB]  -3.483  0.1145  
#    [AB|B]  -3.62   0.09984 
#    [A|B]   -3.688  0.09323 
# 
# At node N9:
#    split   lnL     Rel.Prob
#    [A|A]   -2.844  0.2168  
#    [AB|A]  -2.979  0.1895  
#    [B|A]   -3.158  0.1585  
#    [B|AB]  -3.564  0.1056  
#    [A|AB]  -3.564  0.1056  
#    [B|B]   -3.769  0.08602 
#    [AB|B]  -3.903  0.07517 
#    [A|B]   -4.082  0.06288 
# 
# At node N8:
#    split   lnL     Rel.Prob
#    [AB|A]  -2.199  0.4132  
#    [A|A]   -2.392  0.3407  
#    [B|A]   -2.718  0.2461  

 