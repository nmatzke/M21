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
wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/2species_A_1AB_ML_DECvDECj_BGB/"
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
# This is the example Newick file for Hawaiian A_1AB
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = "tree_2species.newick"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
tr
plot(tr)
title("Example 2-species phylogeny, with branchlengths of 1")
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

# This is the example geography file for Hawaiian A_1AB
# (from Ree & Smith 2008)
geogfn = "geog_2species.data"

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
resfn = "A_1AB_DEC_M0_unconstrained_v1.Rdata"
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
# -1.791759

exp(resDEC$total_loglikelihood)
# 0.1666667
1/6
# 0.1666667


resDEC$outputs@params_table["d","est"]
resDEC$outputs@params_table["e","est"]
# 7.683213e-09
# 6.402844e-09


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

resfn = "A_1AB_DEC+J_M0_unconstrained_v1.Rdata"
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
# 0.1541502

exp(resDECj$total_loglikelihood)
# 1.166666
7/6
# 1.1666667


resDECj$outputs@params_table["d","est"]
resDECj$outputs@params_table["e","est"]
resDECj$outputs@params_table["j","est"]
# 1e-12
# 1e-12
# 2.99999


#######################################################
# PDF plots
#######################################################
pdffn = "AB_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on AB M0_unconstrained"

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
analysis_titletxt ="BioGeoBEARS DEC+J on AB M0_unconstrained"

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
# Dig into code
#######################################################
skip_optim=FALSE


if (is.null(BioGeoBEARS_run_object$allow_null_tips))
	{
	BioGeoBEARS_run_object$allow_null_tips = FALSE
	}

# 2017-11-29
# Error check for tr if not loaded elsewhere
if (exists("tr") == FALSE)
	{
	tr = read.tree(BioGeoBEARS_run_object$trfn)
	}


#######################################################
# Check for traits and trait model
#   - Need BOTH, or NEITHER
#######################################################
traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE

# Initialize m, if needed
m = NULL

if (is.null(BioGeoBEARS_run_object$min_branchlength) == FALSE)
	{
	min_branchlength = BioGeoBEARS_run_object$min_branchlength
	} else {
	min_branchlength = 0.000001
	}

# ERROR CHECKS FOR TRAIT MODEL
if (traitTF)
	{
	if (is.na(BioGeoBEARS_run_object$timesfn) == FALSE)
		{
		txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you have a timesfn, indicate a time-stratified analysis. Traits-based dispersal has not yet been implemented for time-stratified analyses."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)			
		}
	
	# Check for trait_Pmat_txt
	if (is.null(BioGeoBEARS_run_object$trait_Pmat_txt) == TRUE)
		{
		txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you are missing a BioGeoBEARS_run_object$trait_Pmat_txt"
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	# Check for trait transition rates
	# Check for t12, t21, etc.
	trait_transition_rates_TF = grepl(pattern="trait transition rate", x=BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$desc)
	if (sum(trait_transition_rates_TF) < 1)
		{
		txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you need one or more 'trait transition rates' in  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table"
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	# Check for trait-based dispersal rate multipliers		
	# Check for m1, m2, etc.
	numtraitstates = ncol(BioGeoBEARS_run_object$trait@df)
	traitbased_dispersal_Ms_TF = grepl(pattern="trait-based dispersal rate multiplier", x=BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$desc)

	if (sum(traitbased_dispersal_Ms_TF) != numtraitstates)
		{
		txt = paste0("STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, and it has ", numtraitstates, " states, so you need to have ", numtraitstates, " multipliers ('m1', 'm2', etc.) with 'desc' field 'trait-based dispersal rate multipliers...' in  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table. Instead, you have only this many: ", sum(traitbased_dispersal_Ms_TF))
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	} # END if (traitTF) # ERROR CHECK


# Load the trait as a (another) tipranges-class object
if (traitTF == TRUE)
	{
	trait = BioGeoBEARS_run_object$trait
	trait_as_tip_condlikes = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges=trait, phy=tr, states_list=NULL, maxareas=1, include_null_range=FALSE, useAmbiguities=TRUE, trait_as_tip_condlikes=NULL)
	
	# Number of traits
	ntrait_states = ncol(trait_as_tip_condlikes)
	
	# Extract these submatrices, just for dimensions, names etc.
	# ALWAYS extract parameter values from the main model_object
	# Trait modeling effect on dispersal
	# (m1, m2, etc.)
	BGB_trait_model_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[traitbased_dispersal_Ms_TF,]
	
	# Parameters of transition matrix for the trait
	# (t12, t21, etc.)
	BGB_trait_Pmat_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[trait_transition_rates_TF,]

	# Text transition matrix for the trait
	trait_Pmat_txt = BioGeoBEARS_run_object$trait_Pmat_txt
	} else {
	# No trait; set to NULL
	trait_as_tip_condlikes = NULL
	ntrait_states = NULL
	BGB_trait_model_params_table = NULL
	BGB_trait_Pmat_params_table = NULL
	} # END if (traitTF == TRUE)



#######################################################
# Load the model object
#######################################################
#inputs = BioGeoBEARS_run_object
BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object

# Should the optim run be printed?
print_optim = BioGeoBEARS_run_object$print_optim


# Get geographic ranges at tips
if (BioGeoBEARS_run_object$use_detection_model == FALSE)
	{
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
	}
if (BioGeoBEARS_run_object$use_detection_model == TRUE)
	{
	if (BioGeoBEARS_run_object$use_detection_model == TRUE)
		{
		tipranges = tipranges_from_detects_fn(detects_fn=BioGeoBEARS_run_object$detects_fn)
		} # END if (inputs$use_detection_model == TRUE)
	} # END if (BioGeoBEARS_run_object$use_detection_model == TRUE)


# Should we do optimx speedup?
speedup = BioGeoBEARS_run_object$speedup


# Get the list of geographic areas
#	print("print(tipranges):")
#	print(tipranges)
areas = getareas_from_tipranges_object(tipranges)
areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

# Change the names to tipranges@df:
# This converts the tipranges names to 0-based index
# this doesn't make sense if areas_list is 0-based indexes
# XXX - check at some point
# REMOVED: 2017-03-14
#names(tipranges@df) = areas_list

#######################################################
# Set the maximum range size (can be thought of as
# a free parameter
#######################################################
if (is.na(BioGeoBEARS_run_object$max_range_size))
	{
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		# Maximum range size is all areas
		max_range_size = length(areas)
		} else {
		# If not NA
		# Get max rangesize from states list
		max_range_size = max(sapply(X=BioGeoBEARS_run_object$states_list, FUN=length), na.rm=TRUE)
		}
	} else {
	# Maximum range size hard-coded
	max_range_size = BioGeoBEARS_run_object$max_range_size
	}
max_numareas = max_range_size

#######################################################
# Check that no tips have larger ranges than you allowed
#######################################################
#print("Here")
#print(tipranges@df)
# The dfnums_to_numeric fails, if you re-labeled the area names to 0, 1, 2, etc...
tipranges_df_tmp = tipranges@df
names(tipranges_df_tmp) = paste0("col", names(tipranges_df_tmp))
tipranges_df_tmp[tipranges_df_tmp=="?"] = 0
TF = (rowSums(dfnums_to_numeric(tipranges_df_tmp))) > max_range_size
if (sum(TF, na.rm=TRUE) > 0)
	{
	cat("\n\nERROR: Tips with ranges too big:\n", sep="")
	print(dfnums_to_numeric(tipranges_df_tmp)[TF, ])
	cat("\n\nCheck your input geography file!\n", sep="")
	txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", max_range_size, sep="")
	stop(txt)
	}






# Old/slow way of getting the list of states and speciation matrix (slow)
# old_states_list = areas_list_to_states_list(areas, include_null_range=BioGeoBEARS_run_object$include_null_range)
# old_states_list
# spmat = make_relprob_matrix_bi(old_states_list)
# spmat

# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
# max_numareas


# Take the list of areas, and get list of possible states
# (the user can manually input states if they like)
if (is.null(BioGeoBEARS_run_object$states_list))
	{
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=BioGeoBEARS_run_object$include_null_range)
	states_list
	#BioGeoBEARS_run_object$states_list = states_list
	#inputs$states_list = states_list
	} else {
	states_list = BioGeoBEARS_run_object$states_list
	#BioGeoBEARS_run_object$states_list = states_list
	#inputs$states_list = states_list
	}

#######################################################
# NON-STRATIFIED: Modify the states_list if needed
#######################################################
if ( is.numeric(BioGeoBEARS_run_object$timeperiods) == FALSE )
	{
	#######################################################
	# If needed, modify the states_list by areas_allowed_mat
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == FALSE))
		{
		# Take the first areas_allowed matrix (non-stratified)
		areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[1]]
	
		# Cut down the states accordingly (hopefully not slow!)
		original_states_list = states_list
		states_list = prune_states_list(states_list_0based_index=states_list, areas_allowed_mat=areas_allowed_mat)
		BioGeoBEARS_run_object$states_list = states_list

		print("Limiting original_states_list using an areas_allowed matrix")
		print("original_states_list")
		print(original_states_list)
		cat("\nlength(original_states_list) = ", length(original_states_list), " states/ranges.\n")
		cat("\n")

		print("states_list")
		print(states_list)
		cat("\nlength(original_states_list) = ", length(original_states_list), " states/ranges.")
		cat("\nlength(states_list) = ", length(states_list), " states/ranges.\n")


		} else {
		# Make no change
		pass = 1
		# states_list = states_list
		}

	#######################################################
	# If needed, modify the states_list by areas_adjacency_mat
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == FALSE))
		{
		# Take the first areas_adjacency matrix (non-stratified)
		areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[1]]
	
		# Cut down the states accordingly (hopefully not slow!)
		original_states_list = states_list
		states_list = prune_states_list_by_adjacency(states_list_0based_index=states_list, areas_adjacency_mat=areas_adjacency_mat)
		BioGeoBEARS_run_object$states_list = states_list
		
		print("Limiting original_states_list using an area adjacency matrix")
		print("original_states_list")
		print(original_states_list)
		print(length(original_states_list))
		cat("\n")

		print("states_list")
		print(states_list)
		print("length(states_list)")
		print(length(states_list))
		
		} else {
		# Make no change
		pass = 1
		# states_list = states_list
		}
	} # END if ( is.numeric(BioGeoBEARS_run_object$timeperiods) == FALSE )

# Change the states_list by traits, if needed
# (non-stratified)
if (traitTF == TRUE)
	{
	states_list_ORIG = states_list
	#states_list_wTrait = 
	
	#trait_as_tip_condlikes
	}



#######################################################
# STRATIFIED: Modify the states_list if needed
# (this is ONLY if the state-space is changing in the
#  different time-slices)
#######################################################
# Will the state space be changing?
TF1 = (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == FALSE)
TF2 = (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == FALSE)
state_space_changing_TF = (TF1 + TF2) > 0
need_to_print_list_of_states_list = TRUE
master_states_list = states_list	# store the master list of all states;
									# check that this includes all, at some point
									# if not, warn user to change it manually

if ( (is.numeric(BioGeoBEARS_run_object$timeperiods) == TRUE) && (state_space_changing_TF == TRUE) && (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == TRUE) )
	{
	need_to_print_list_of_states_list = FALSE
	ntimes = length(BioGeoBEARS_run_object$timeperiods)
	lists_of_states_lists_0based = list()
	
	# Go through each time bin, and make the state space different in each time bin
	for (ti in 1:ntimes)
		{
		# Initialize
		states_list_for_this_stratum = states_list
		
		
		#######################################################
		# If needed, modify the states_list by areas_allowed_mat
		#######################################################
		# Areas allowed matrix
		if (TF1 == TRUE)
			{
			# Take the first areas_allowed matrix (non-stratified)
			areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[ti]]
	
			# Cut down the states accordingly (hopefully not slow!)
			states_list_for_this_stratum = prune_states_list(states_list_0based_index=states_list_for_this_stratum, areas_allowed_mat=areas_allowed_mat)
			} else {
			# Make no change
			pass = 1
			# states_list = states_list
			}

		# Message to user
		timeslice_num = ti
		if (timeslice_num == 1)
			{
			toptime = 0
			} else {
			toptime = BioGeoBEARS_run_object$timeperiods[ti-1]
			}
		if (timeslice_num == ntimes)
			{
			bottime = BioGeoBEARS_run_object$timeperiods[ti]
			catend = "\n\n"
			} else {
			bottime = BioGeoBEARS_run_object$timeperiods[ti]
			catend = ""
			}
		txt = paste0("bears_optim_run() note: overall states_list has ", length(master_states_list), " states/ranges. In stratum #", ti, " (", toptime, "-", bottime, " mya), states_list_for_this_stratum has ", length(states_list_for_this_stratum), " states/ranges, due to areas_allowed and/or areas_adjacency matrices. See BioGeoBEARS_run_object$lists_of_states_lists_0based.")
		cat("\n")
		cat(txt)
		cat(catend)


		#######################################################
		# If needed, modify the states_list by areas_adjacency_mat
		#######################################################
		# Areas adjacency matrix
		if (TF2 == TRUE)
			{
			# Take the first areas_adjacency matrix (non-stratified)
			areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[ti]]
	
			# Cut down the states accordingly (hopefully not slow!)
			states_list_for_this_stratum = prune_states_list_by_adjacency(states_list_0based_index=states_list_for_this_stratum, areas_adjacency_mat=areas_adjacency_mat)
			} else {
			# Make no change
			pass = 1
			# states_list = states_list
			}

		# Store in the list of states_lists
		lists_of_states_lists_0based[[ti]] = states_list_for_this_stratum
		} # END for (ti in 1:ntimes)
	
	# Store the time-stratified list of states_lists in the BioGeoBEARS_run_object
	BioGeoBEARS_run_object$lists_of_states_lists_0based = lists_of_states_lists_0based
	
	} # END if ( (is.numeric(BioGeoBEARS_run_object$timeperiods) == TRUE) && (state_space_changing_TF == TRUE) )


# Or, if the time-stratified stats list is pre-specified
if (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == FALSE)
	{
	ntimes = length(BioGeoBEARS_run_object$timeperiods)
	txt = paste0("bears_optim_run() note: BioGeoBEARS_run_object$lists_of_states_lists_0based has been specified. This means there is a different state space in each timebin / stratum / epoch.")
	cat("\n")
	cat(txt)
	cat("\n")
	
	# Check that number of lists of states matches the number of timebins
	number_of_lists_of_states = length(BioGeoBEARS_run_object$lists_of_states_lists_0based)
	if (ntimes == number_of_lists_of_states)
		{
		txt = paste0("bears_optim_run() note: BioGeoBEARS_run_object has ", ntimes, " timebins and ", number_of_lists_of_states, " lists of states ranges. Check passed.")
		cat("\n")
		cat(txt)
		cat("\n")
		} else {
		txt = paste0("bears_optim_run() STOP ERROR: BioGeoBEARS_run_object has ", ntimes, " timebins and ", number_of_lists_of_states, " lists of states ranges. Check FAILED.")
		cat("\n")
		cat(txt)
		cat("\n")
		stop(txt)
		} # END if (ntimes = number_of_lists_of_states)

	
	
	# Go through each time bin, and make the state 
	# space different in each time bin
	if (need_to_print_list_of_states_list == TRUE)
		{
		for (ti in 1:ntimes)
			{
			# Extract the states list in this time-stratum
			states_list_for_this_stratum = BioGeoBEARS_run_object$lists_of_states_lists_0based[[ti]]
			
			# Message to user
			timeslice_num = ti
			if (timeslice_num == 1)
				{
				toptime = 0
				} else {
				toptime = BioGeoBEARS_run_object$timeperiods[ti-1]
				}
			if (timeslice_num == ntimes)
				{
				bottime = BioGeoBEARS_run_object$timeperiods[ti]
				catend = "\n\n"
				} else {
				bottime = BioGeoBEARS_run_object$timeperiods[ti]
				catend = ""
				} # END if (timeslice_num == ntimes)				
				
				
			txt = paste0("bears_optim_run() note: overall states_list has ", length(master_states_list), " states/ranges. In stratum #", ti, " (", toptime, "-", bottime, " mya), states_list_for_this_stratum has ", length(states_list_for_this_stratum), " states/ranges, due to user-specified states_lists. See BioGeoBEARS_run_object$lists_of_states_lists_0based.")
			cat("\n")
			cat(txt)
			cat(catend)
			} # END for (ti in 1:ntimes)
		} # END if (need_to_print_list_of_states_list == TRUE)
	} # END if (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == TRUE)
	# END printing user-specified list of states_lists



#######################################################
# Sparse matrix exponentiation, if desired (dubious)
#######################################################
if (is.na(BioGeoBEARS_run_object$force_sparse))
	{
	if (length(states_list) > 128)
		{
		force_sparse = TRUE
		cat("\nNote: force_sparse being set to TRUE, as length(states_list) > 128\n", sep="")
		} else {
		force_sparse = FALSE
		}
	} else {
	force_sparse = BioGeoBEARS_run_object$force_sparse
	}

if (force_sparse == TRUE)
	{
	cat("\nNote: force_sparse is set to TRUE; length(states_list)=", length(states_list), "\n", sep="")
	}



#######################################################
# Load the phylogenetic tree
#######################################################
trfn = np(BioGeoBEARS_run_object$trfn)
#phy = read.tree(file=trfn)
phy = check_trfn(trfn=trfn)

# The likelihood of each state at the tips
# Change this, if you have observations instead of presence/absence at the tips

# Options:
# 1. Use tipranges_to_tip_condlikes_of_data_on_each_state ()
# 2. Use detection model to generate tip likelihoods if desired; or 
# 3. Take pre-specified tip likelihoods. 
if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == TRUE)
	{
	if (BioGeoBEARS_run_object$use_detection_model == FALSE)
		{
		#print("here2")
		#print(states_list)
		
		tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas, include_null_range=BioGeoBEARS_run_object$include_null_range, useAmbiguities=BioGeoBEARS_run_object$useAmbiguities, trait_as_tip_condlikes=trait_as_tip_condlikes, allow_null_tips=BioGeoBEARS_run_object$allow_null_tips)
		} else {
		# Calculate the initial tip likelihoods, using the detection model
		# Assumes correct order, double-check this
		numareas = length(areas)
		detects_df = BioGeoBEARS_run_object$detects_df
		controls_df = BioGeoBEARS_run_object$controls_df
		mean_frequency = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mf", "init"]
		dp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["dp", "init"]
		fdp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["fdp", "init"]
	
		# return_LnLs=TRUE ensures no under-flow
		tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index=states_list, phy=phy, numareas=numareas, detects_df=detects_df, controls_df=controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, null_range_gets_0_like=TRUE, return_LnLs=TRUE, relative_LnLs=TRUE, exp_LnLs=TRUE, error_check=TRUE)
		}
	} else {
	# Or, use pre-specified tip conditional likelihoods
	# Pre-specified (custom) tip-likelihoods
	tip_condlikes_of_data_on_each_state = BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state
	} # END if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == FALSE)

#print(tip_condlikes_of_data_on_each_state)

if (is.null(BioGeoBEARS_run_object$printlevel))
	{
	BioGeoBEARS_run_object$printlevel = 1
	}
printlevel = BioGeoBEARS_run_object$printlevel


#######################################################
# Read the stratification/distances input files, if any
#######################################################
#inputs = readfiles_BioGeoBEARS_run(inputs=BioGeoBEARS_run_object)

#######################################################
# Check for problems in the input files; will throw stop() if there are problems
#######################################################
#check_result = check_BioGeoBEARS_run(inputs=BioGeoBEARS_run_object)
#check_result


#######################################################
# Set up the function for optimization
#######################################################	
# params are a list of the values of the FREE parameters; but everything is contained in the 
# BioGeoBEARS_model object at all times
# (moved to separate function)


# defaults for optimization
# We are using "L-BFGS-B", which is:
#####################################################################################################
# Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each
# variable can be given a lower and/or upper bound. The initial value must satisfy the constraints.
# This uses a limited-memory modification of the BFGS quasi-Newton method. If non-trivial bounds
# are supplied, this method will be selected, with a warning.
# 
# [...]
#
# Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995) A limited memory algorithm for bound constrained
# optimization. SIAM J. Scientific Computing, 16, 1190–1208.
#####################################################################################################
#
# "BGFS" refers to: 4 articles, Broyden, Fletcher, Goldfarb and Shanno (1970).

# Run check, before rescaling
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#######################################################
# 2016-03-23_NJM: adding rescaling
#######################################################
if (BioGeoBEARS_run_object$rescale_params == TRUE)
	{
	BioGeoBEARS_model_object@params_table = scale_BGB_params(orig_params_table=BioGeoBEARS_model_object@params_table, add_smin=0, add_smax=1)
	
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	}



params = BioGeoBEARS_model_object_to_init_params(BioGeoBEARS_model_object)
minj = 1e-05
# start on Lagrange results
#params = c(3.11882,  2.51741)
lower = BioGeoBEARS_model_object_to_params_lower(BioGeoBEARS_model_object)
upper = BioGeoBEARS_model_object_to_params_upper(BioGeoBEARS_model_object)

# High performance computing
# HPC using parallel package in R 2.15 or higher, which allows
# mcmapply (multicore apply)
# Don't use multicore if using R.app ('AQUA')
num_cores_to_use = BioGeoBEARS_run_object$num_cores_to_use
cluster_already_open = BioGeoBEARS_run_object$cluster_already_open

cluster_was_open = FALSE
if (.Platform$GUI != "AQUA" && ((is.na(num_cores_to_use) == TRUE) || ( (is.na(num_cores_to_use)==FALSE) && (num_cores_to_use > 1))) )
	{
	# We are doing manual, optional processing on several cores;
	# this seems to have less overhead/hassle/incompatibility issues
	# than using mcmapply, mclapply, etc...
	#require("parallel") #<- do this higher up

	num_cores_computer_has = detectCores()
	txt = paste0("Your computer has ", num_cores_computer_has, " cores.")
	cat("\n")
	cat(txt)
	cat("\n")		

	if (is.null(num_cores_to_use))
		{
		num_cores_to_use = num_cores_computer_has
		}

	if (num_cores_to_use > num_cores_computer_has)
		{
		txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", num_cores_to_use, " cores, but your computer only has ", num_cores_computer_has, ". Resetting to ", num_cores_computer_has, ".")
		cat("\n")
		cat(txt)
		cat("\n")
		warning(txt)
		num_cores_to_use = num_cores_computer_has
		}
	
	# Don't do this, if the cluster is already open
	cat("\nYour computer has ", num_cores_computer_has, " cores. You have chosen to use:\nnum_cores_to_use = ", num_cores_to_use, " cores for the matrix exponentiations in the likelihood calculations.\n", sep="")
	
	if ( is.logical(cluster_already_open) == TRUE )
		{
		if (cluster_already_open == FALSE)
			{
			cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")
			cat("Started cluster with ", num_cores_to_use, " cores.\n\n", sep="")

			# Flag so that you remember to close cluster at the end
			cluster_open=TRUE
			cluster_was_open = FALSE
			}
		} else {
		cluster_was_open = TRUE
		cat("Cluster with ", num_cores_to_use, " cores already open.\n\n", sep="")
		}
	} else {
	# You are using R.app and clusters don't work...

	num_cores_computer_has = detectCores()
	txt = paste0("Your computer has ", num_cores_computer_has, " cores.")
	cat("\n")
	cat(txt)
	cat("\n")	
	
	if (num_cores_to_use > 1)
		{
		txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", num_cores_to_use, " cores, but in R.app, multicore functionality doesn't work. Resetting num_cores_to_use=1.")
		cat("\n")
		cat(txt)
		cat("\n")
		warning(txt)
		num_cores_to_use = num_cores_computer_has
		}
	
	
	cluster_already_open = NULL
	cluster_was_open = FALSE
	}





#######################################################
#######################################################
# NON-stratified analysis
#######################################################
#######################################################
# Run optimization on a SINGLE tree
use_optimx = BioGeoBEARS_run_object$use_optimx

# IN OPTIMX, bobyqa method:
# no control on tolerance

# Try parscale:
# parscale: A vector of scaling values for the parameters. 
# Optimization is performed on par/parscale and these should 
# be comparable in the sense that a unit change in any element 
# produces about a unit change in the scaled value.For optim.
# https://www.mail-archive.com/r-help@r-project.org/msg152890.html
# "(optimx includes parscale on all methods)."
parscale = (upper - lower) / min(upper - lower)
print("parscale:")
print(parscale)
control_list = list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE)
# old:
# factr=0.0001)#, reltol=0.001)#, maxit=100)

# Bogus note (NJM):
# This causes pathology: reltol=0.001
# Actually, this was fine, it was 
# force_sparse = TRUE that was the problem
# (leads to different results!!  probably rounding errors)

# Limit the number of iterations so it 
# doesn't go on forever
num_free_params = sum(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$"type" == "free")
num_free_params
itnmax = 50 * num_free_params

	# Un-comment only for error checking, then re-comment!!!!!!!!!!!!!!
	cat("\n\nNOTE: Before running optimx(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside optimx().\n\n", sep="")

	
	
	loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)

params2 = params
params2[1] = 1.000000e-12
params2[2] = 1.000000e-12
params2[3] = 2.99999

loglike = calc_loglike_for_optim(params=params2, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
loglike


params3 = params
params3[1] = 1.000000e-12
params3[2] = 1.000000e-12
params3[3] = 1.5

loglike = calc_loglike_for_optim(params=params3, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
loglike


params4 = params
params4[1] = 1.000000e-12
params4[2] = 1.000000e-12
params4[3] = 0.05

loglike = calc_loglike_for_optim(params=params4, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
loglike






	
	cat("\ncalc_loglike_for_optim() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optimx()...\n\n", sep="")

	cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim():\n\n")
	print(warnings())


	if (skip_optim == TRUE)
		{
		# Skip optimization
		cat("Skipping ML search as skip_optim==TRUE.\n\n", sep="")
		return(loglike)
		}
	
	# optimx 2012 versus 2013
	if (packageVersion("optimx") < 2013)
		{
	# optimx 2012
	optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, gr=NULL, hess=NULL, lower=lower, upper=upper, method=c("bobyqa"), itnmax=itnmax, hessian=NULL, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	# old:
	# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
		} else {
	# Run optimx scalecheck
	scalecheck_results = optimx:::scalecheck(par=params, lower=lower, upper=upper)
	
	cat("\n\nResults of optimx:::scalecheck() below. Note: sometimes rescaling parameters may be helpful for ML searches, when the parameters have much different absolute sizes. This can be attempted by setting BioGeoBEARS_run_object$rescale_params = TRUE.\n\n")
	print(scalecheck_results)
	

	# optimx 2013
	optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, gr=NULL, hess=NULL, lower=lower, upper=upper, method=c("bobyqa"), itnmax=itnmax, hessian=FALSE, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))				
		} # end packageVersion



	# Run with all methods, for testing:
	# optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

	#######################################################
	# Compare optimization routines
	#######################################################
	
	# BEARS_results_7areas_2param
	#                        par   fvalues   method fns grs itns conv  KKT1 KKT2  xtimes
	# 6 0.010165217, 0.009422923 -57.81254   bobyqa  30  NA NULL    0 FALSE TRUE  29.761
	# 1 0.010180679, 0.009492254 -57.81255 L-BFGS-B  33  33 NULL    0 FALSE TRUE 151.137
	# 4   0.01017895, 0.00970706 -57.81263   Rcgmin  32   7 NULL    0 FALSE TRUE  42.461
	# 2 0.010242504, 0.009822486 -57.81284   nlminb  40   6    3    1 FALSE TRUE  43.399
	# 3   0.01017850, 0.01001962 -57.81293      spg  68  NA   47    0 FALSE TRUE 150.932
	# 5               0.01, 0.01 -57.81366   Rvmmin  20   1 NULL    0 FALSE TRUE  21.456


	#return (optim_result2)
	} # END if ( (use_optimx == TRUE) || (use_optimx == "optimx") )


# Try GenSA for more complex optimization problems (4+ parameters, or 
# wildly different parameter scalings)
if (use_optimx == "GenSA")
	{
	require(GenSA)
	
	cat("\n\nNOTE: You are optimizing with GenSA() ('Generalized Simulated Annealing') instead of optimx() or optim(). GenSA may be better for more complex problems (4+ parameters, wildly different scalings), but has not been extensively tested for BioGeoBEARS yet. And it may be slower.")
	
	# Un-comment only for error checking, then re-comment!!!!!!!!!!!!!!
	cat("\n\nNOTE: Before running GenSA(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside GenSA().\n\n", sep="")

	loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	
	cat("\ncalc_loglike_for_optim() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with GenSA()...\n\n", sep="")

	cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim():\n\n")
	print(warnings())

	if (skip_optim == TRUE)
		{
		# Skip optimization
		cat("Skipping ML search as skip_optim==TRUE.\n\n", sep="")
		return(loglike)
		}
	
	control_list = list(nb.stop.improvement=50, simple.function=TRUE, trace.mat=TRUE)			

	# NJM: I am assuming that the functions are fairly smooth in BioGeoBEARS analyses
	if (is.null(BioGeoBEARS_run_object$temperature) == FALSE)
		{
		temperature = BioGeoBEARS_run_object$temperature
		control_list = c(control_list, list(temperature=temperature))
		}
	if (is.null(BioGeoBEARS_run_object$max.call) == FALSE)
		{
		max.call = BioGeoBEARS_run_object$max.call
		control_list = c(control_list, list(max.call=max.call))
		} else {
		max.call = length(params) * 250
		control_list = c(control_list, list(max.call=max.call))
		}
				
	optim_result2 = GenSA(par=params, fn=calc_loglike_for_optim_neg, lower=lower, upper=upper, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	} # END if ( (use_optimx == TRUE) || (use_optimx == "optimx") )



# optimx 2013
optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, gr=NULL, hess=NULL, lower=lower, upper=upper, method=c("bobyqa"), itnmax=itnmax, hessian=FALSE, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)

if (skip_optim == TRUE)
	{
	# Skip optimization
	cat("Just returning initial loglike as skip_optim==TRUE.\n\n", sep="")
	return(loglike)
	}

# Update the parameter values in the output BioGeoBEARS_model_object using
# the ML results
optimx_result = optim_result2
use_optimx = BioGeoBEARS_run_object$use_optimx

cat("\n\nThis is the output from optim, optimx, or GenSA. Check the help on those functions to\ninterpret this output and check for convergence issues:\n\n")
print(optimx_result)

cat("\n\nReading the optim/optimx/GenSA output into the BioGeoBEARS_model object:\n\nBioGeoBEARS_model_object =\n\n")

BioGeoBEARS_model_object = update_BioGeoBEARS_model_object_w_optimx_result(BioGeoBEARS_model_object, optimx_result, use_optimx)


print(BioGeoBEARS_model_object)

cat("\n\n...successful.\n\n")

# Update the output
BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object

params = BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_model_object)
#print(params)
# Calculate the log-likelihood of the data, given the model parameters during this iteration	
#model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=1, calc_ancprobs=TRUE, include_null_range=BioGeoBEARS_run_object$include_null_range)

# We need to put the params back into the inputs 
# to get the reconstructed ancestors etc.
# Note that fixlikes SHOULD be included here in the 
# final results, if specified by the user at the beginning
# (thanks to Julien for pointing out issue)
calc_ancprobs = BioGeoBEARS_run_object$calc_ancprobs
# (originally, to do local ancestral states, you would set
#  calc_ancprobs to FALSE and use the subsequent optim_result
#  $ fvalue to get the LnL optimal on that node state.
#  I am now changing it to always use the fixlikes.)

#print("Calculating final LnL...")
model_results = calc_loglike_for_optim(params=params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="all", calc_ancprobs=calc_ancprobs)





	if (is.null(BioGeoBEARS_run_object$printlevel))
		{
		BioGeoBEARS_run_object$printlevel = 1
		}
	printlevel = BioGeoBEARS_run_object$printlevel
	

	# Is this a traits-based analysis?
	traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE
	
	# Initialize m
	m = NULL

	# Initialize jts_matrix, matrix of t12, t23, etc., during a j event
	jts_matrix = NULL

	
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	#print(params)
	#print(BioGeoBEARS_model_object)
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)

	######################################################
	# 2016-03-23_NJM: adding rescaling
	# (unscale params, if they were used before)
	######################################################
	if (BioGeoBEARS_run_object$rescale_params == TRUE)
		{
		BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table=BioGeoBEARS_model_object@params_table)
		
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_model_object@params_table
		}

	
	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Update to the run object, just to be SURE
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]


	#######################################################
	#######################################################
	# Do branch-length exponentiation if desired
	#######################################################
	#######################################################
	b = BioGeoBEARS_model_object@params_table["b","est"]
	# Modify the edge.lengths
	phy$edge.length = phy$edge.length ^ b
	# Make sure this doesn't duplicate a previous "^b", e.g.
	# the summarization step in bears_optim_run

	#######################################################
	#######################################################
	# Do distance-dependence and dispersal multipliers matrix
	#######################################################
	#######################################################
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
	
	# If there is a distance matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_distances_mats) == FALSE))
		{
		distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on distance, apply to distances matrix
	x = BioGeoBEARS_model_object@params_table["x","est"]
	dispersal_multipliers_matrix = distances_mat ^ x

	# Environmental distances
	if ( (is.null(BioGeoBEARS_run_object$list_of_envdistances_mats) == FALSE))
		{
		envdistances_mat = BioGeoBEARS_run_object$list_of_envdistances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		envdistances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on environmental distance, apply to distances matrix
	n = BioGeoBEARS_model_object@params_table["n","est"]
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * (envdistances_mat ^ n)

	# Apply manual dispersal multipliers, if any
	# If there is a manual dispersal multipliers matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == FALSE))
		{
		manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		manual_dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))
		}
	
	# Get the exponent on manual dispersal multipliers
	w = BioGeoBEARS_model_object@params_table["w","est"]

	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * (manual_dispersal_multipliers_matrix ^ w)

	#######################################################
	# ALTERNATIVE DISTANCE FUNCTIONS
	#
	# The exponential function of distance is heavy-tailed.
	# Other functions might be empirically better
	# 
	# NOTE: in order to keep the base dispersal rate
	# (base dispersal rate = rate at distance=0) comparable
	# across models, the probability density function should 
	# be normalized (/divided by) the pdf evaluated at zero.
	#
	# In other words, the dispersal multiplier should always
	# come out as 1, when distance=0. I.e. pdf(x=0) = 1.
	# 
	# These will each be identified by optional parameters
	# in the model:
	# WALD
	# HNORM
	#
	#######################################################
	
	#######################################################
	#	WALD distribution
	#######################################################
	#
	# If one has a model where dispersal rate is modified as 
	# a function of the WALD distribution (inverse gaussian),
	# this will require parameters:
	# WALD_mu - mu, the mean in the WALD probability density function
	# WALD_lambda - lambda, the shape parameter
	# See: https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
	#######################################################
	
	# Alternative distance model abbreviations
	alt_distance_models_abbr = c("WALD", "HNORM")
	tmp_param_names = row.names(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)

	# Calculate additional dispersal modifiers on distance
	# WALD distribution
	TF = grepl(pattern="WALD", x=tmp_param_names)
	if (sum(TF) > 0)
		{
		require(statmod)	# for dinvgauss
		tmpname = "WALD_mu"
		TF = grepl(pattern=tmpname, x=tmp_param_names)
		WALD_mu = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, "est"]

		tmpname = "WALD_lambda"
		TF = grepl(pattern=tmpname, x=tmp_param_names)
		WALD_lambda = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, "est"]
		
		# Calculate the multipliers under this model
		tmp_multipliers = dinvgauss(x=distances_mat, mean=WALD_mu, shape=WALD_lambda)
		normalizer = dinvgauss(x=0, mean=WALD_mu, shape=WALD_lambda)
		tmp_multipliers = tmp_multipliers / normalizer
		
		# Multiply element-wise
		dispersal_multipliers_matrix = dispersal_multipliers_matrix * tmp_multipliers
		} # END if (sum(TF) > 0)


	# Half-normal distribution
	# https://en.wikipedia.org/wiki/Half-normal_distribution
	# Theta is the scaled precision (inverse of the variance)
	# Theta = sqrt(pi) / (sigma * sqrt(2))
	TF = grepl(pattern="HNORM", x=tmp_param_names)
	if (sum(TF) > 0)
		{
		require(fdrtool)	# for dhalfnorm
		tmpname = "HNORM_theta"
		TF = grepl(pattern=tmpname, x=tmp_param_names)
		HNORM_theta = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, "est"]


		# Calculate the multipliers under this model
		tmp_multipliers = dhalfnorm(x=distances_mat, theta=HNORM_theta)
		normalizer = halfnorm(x=0, theta=HNORM_theta)
		tmp_multipliers = tmp_multipliers / normalizer
		
		# Multiply element-wise
		dispersal_multipliers_matrix = dispersal_multipliers_matrix * tmp_multipliers
		} # END if (sum(TF) > 0)



	#######################################################
	# multiply parameter d by dispersal_multipliers_matrix
	#######################################################
	dmat_times_d = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
	amat = dispersal_multipliers_matrix * matrix(a, nrow=length(areas), ncol=length(areas))




	#######################################################
	#######################################################
	# Do area-dependence and extinction multipliers list
	#######################################################
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == FALSE))
		{
		area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[1]]
		} else {
		# Default is all areas effectively equidistant
		area_of_areas = rep(1, length(areas))
		}
		
	# Get the exponent on extinction, apply to extinction modifiers	
	u = BioGeoBEARS_model_object@params_table["u","est"]
	extinction_modifier_list = area_of_areas ^ (1 * u)
	
	# Apply to extinction rate
	elist = extinction_modifier_list * rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	# someday we'll have to put "a" (anagenic range-switching) in here
	# (this was been done in 2014 - NJM)
	
	# Substitute here (for starters) if you want to have a custom Qmat
	if (is.null(BioGeoBEARS_run_object$custom_Qmat_fn_text) == TRUE)
		{
		# Standard analysis, no traits
		if (traitTF == FALSE)
			{
			Qmat = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=BioGeoBEARS_run_object$include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)
			}
		
		# Analysis with a trait modifying dispersal rate
		if (traitTF == TRUE)
			{
			res = modify_Qmat_with_trait(Qmat=NULL, BioGeoBEARS_run_object, areas_list=areas_list, states_list=states_list, dispersal_multipliers_matrix=dispersal_multipliers_matrix, elist=elist, force_sparse=force_sparse)
			Qmat = res$Qmat
			m = res$m

			# If the trait can change during jump events
			if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
				{
				jts_txt_matrix = BioGeoBEARS_run_object$jts_txt_matrix
				jts_matrix = matrix(data=0, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
				TF_matrix = matrix(data=TRUE, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
				diag(TF_matrix) = FALSE
				jts_txt_params = c(jts_txt_matrix[TF_matrix])
				jts_txt_params
			
				# Populate the numeric jts_matrix
				for (jts_i in 1:nrow(jts_txt_matrix))
					{
					diag_val = 1
					for (jts_j in 1:ncol(jts_txt_matrix))
						{
						if (jts_i == jts_j)
							{
							next()
							}
						jts_txt = jts_txt_matrix[jts_i,jts_j]
						newval = as.numeric(BioGeoBEARS_model_object@params_table[jts_txt, "est"])
						jts_matrix[jts_i,jts_j] = newval
						diag_val = 1-newval
						}
					# Populate the diagonal
					jts_matrix[jts_i,jts_i] = diag_val
					} # END for (jts_i in 1:nrow(jts_txt_matrix))
				} # END if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
			} # END if (if (traitTF == TRUE))
		
		} else {
		cat("\n\nNOTE: BioGeoBEARS is using a custom Qmat-generating function.\n\n")
		# Evaluates to "Qmat"
		eval(parse(text=BioGeoBEARS_run_object$custom_Qmat_fn_text))
		} # END if (is.null(BioGeoBEARS_run_object$custom_Qmat_fn_text) == TRUE)
	# Print the Qmat, if desired
#	print(round(Qmat, 4))
	
	
	
	#######################################################
	# Cladogenic model
	#######################################################
	j = BioGeoBEARS_model_object@params_table["j","est"]
	ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
	ys = BioGeoBEARS_model_object@params_table["ys","est"]
	v = BioGeoBEARS_model_object@params_table["v","est"]
	y = BioGeoBEARS_model_object@params_table["y","est"]
	s = BioGeoBEARS_model_object@params_table["s","est"]
	sum_SPweights = y + s + j + v


	maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
	
	# Text version of speciation matrix	
	maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
	maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
	maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
	maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]


	# Cladogenesis model inputs
	spPmat_inputs = NULL

	# Note that this gets the dispersal multipliers matrix, which is applied to 
	# e.g. the j events, NOT the dmat_times_d above which is d*dispersal_multipliers_matrix
	dmat = dispersal_multipliers_matrix
	spPmat_inputs$dmat = dmat

	states_indices = states_list
	
	# shorten the states_indices by 1 (cutting the 
	# null range state from the speciation matrix)
	if (BioGeoBEARS_run_object$include_null_range == TRUE)
		{
		states_indices[1] = NULL
		} # END if (BioGeoBEARS_run_object$include_null_range == TRUE)
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = s
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = y
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param

	# Use detection model to generate tip likelihoods if desired; or take
	# pre-specified tip likelihoods. Otherwise, use those input from bears_optim_run
	if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == TRUE)
		{
		# Usual strategy for calculating tip-likelihoods
		# Get the detection model
		if (BioGeoBEARS_run_object$use_detection_model == TRUE)
			{
			# Calculate the initial tip likelihoods, using the detection model
			# Assumes correct order, double-check this
			numareas = length(areas)
			detects_df = BioGeoBEARS_run_object$detects_df
			controls_df = BioGeoBEARS_run_object$controls_df
			mean_frequency = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mf", "init"]
			dp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["dp", "init"]
			fdp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["fdp", "init"]
		
			# return_LnLs=TRUE ensures no under-flow
			tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index=states_list, phy=phy, numareas=numareas, detects_df=detects_df, controls_df=controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, null_range_gets_0_like=TRUE, return_LnLs=TRUE, relative_LnLs=TRUE, exp_LnLs=TRUE, error_check=TRUE)
			}
		#print(tip_condlikes_of_data_on_each_state)

		} else {
		# Or, use pre-specified tip conditional likelihoods
		# Pre-specified (custom) tip-likelihoods
		tip_condlikes_of_data_on_each_state = BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state
		} # END if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == FALSE)



	if (print_optim == TRUE)
		{
		#outvars = as.data.frame(t(BioGeoBEARS_model_object@params_table$est))
		#names(outvars) = rownames(BioGeoBEARS_model_object@params_table)
		#outvars = c(BioGeoBEARS_model_object@params_table$est)
		
		#cat("\n")
		#cat(outvars, sep="	")
		
		# Before calculating the log likelihood, print it, in case there is e.g. a bug
		#cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent01=", maxent_constraint_01, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
		}



		# Fixing ancestral nodes
		fixnode = BioGeoBEARS_run_object$fixnode
		fixlikes = BioGeoBEARS_run_object$fixlikes
		
		# Print m
		# (m1, m2, etc.)
#		print("m:")
#		print(m)
		
		# Calculate EVERYTHING!
		#print(jts_matrix)
		model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=1, cluster_already_open=cluster_already_open, calc_ancprobs=TRUE, fixnode=fixnode, fixlikes=fixlikes, include_null_range=BioGeoBEARS_run_object$include_null_range, m=m, jts_matrix=jts_matrix, BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, on_NaN_error=BioGeoBEARS_run_object$on_NaN_error)









#######################################################
# calc_loglike_sp
#######################################################

sparse = FALSE

	#print("on_NaN_error:")
	#print(on_NaN_error)
	#print("Qmat#1")
	#print(Qmat)

	#######################################################
	# Error check on fixnode / fixlikes
	#######################################################
	if (!is.null(fixnode))
		{
		if (( is.null(dim(fixlikes)) == TRUE) && (length(fixnode)==1))
			{
			pass_fixlikes = TRUE
			} else {
			if ( (dim(fixlikes)[1]) == length(fixnode) )
				{
				pass_fixlikes = TRUE
				
				# Another error check: Multiple nodes in 'fixnode' MUST be sorted in increasing order
				if ( (all(c(order(fixnode) == 1:length(fixnode)))) == TRUE)
					{
					pass_fixlikes = TRUE
					} else {
					pass_fixlikes = FALSE
					error_msg = "ERROR in calc_loglike_sp(): \n             Multiple nodes in 'fixnode' MUST be sorted in increasing order.\n"
					cat(error_msg)
					stop(error_msg)
					} # end if order check

				} else {
				pass_fixlikes = FALSE
				error_msg = "ERROR: Either:\n             (1) fixnode must be a single node number, and fixlikes must be a vector, or\n             (2) fixlikes like must be a matrix with the # of rows equal to length(fixnode).\n"
				cat(error_msg)
				stop(error_msg)
				} # end 2nd if()
			} # end 1st if()
		}
	
	
	#######################################################
	# Calculating ancestral probs via an uppass requires the
	# downpass splitlikes be saved.
	#######################################################
	if (calc_ancprobs == TRUE)
		{
		# Also, cppSpMethod must be 3 (I'm not going to program the other ones!!)
		if (cppSpMethod != 3)
			{
			stop("ERROR: You must have cppSpMethod=3 if calc_ancprobs=TRUE")
			}
		}
	
#	print("m:")
#	print(m)
	
	#######################################################
	# ERROR CHECK
	#######################################################
	if (use_cpp == TRUE && !is.null(spPmat_inputs))
		{
		numstates_in_tip_condlikes_of_data_on_each_state = ncol(tip_condlikes_of_data_on_each_state)
		
		# Number of states in the spPmat, *including* the null range if desired
		if (include_null_range == TRUE)
			{
			numstates_in_spPmat_inputs_states_indices = 1 + length(spPmat_inputs$l)
			} else {
			numstates_in_spPmat_inputs_states_indices = 0 + length(spPmat_inputs$l)
			} # END if (include_null_range == TRUE)
		
		if (is.null(m) == FALSE)
			{
			numstates_in_spPmat_inputs_states_indices = numstates_in_spPmat_inputs_states_indices * length(m)
			}
		
		if (input_is_COO == TRUE)
			{
			numstates_in_Qmat = max( max(Qmat[,1]), max(Qmat[,2]) )
			} else {
			numstates_in_Qmat = ncol(Qmat)
			} # END if (input_is_COO == TRUE)
		
		if ( all(numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_spPmat_inputs_states_indices, numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_Qmat ) == TRUE )
			{
			# Continue
			all_inputs_correct_size = TRUE
			} else {
			
			# Stop if not everything is equal
			stop_function <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat)
				{
				cat("\n\nERROR: Some inputs have incorrect size -- \n")
				cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
				cat("numstates_in_spPmat_inputs_states_indices+1:	", numstates_in_spPmat_inputs_states_indices, "\n")
				cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
				cat("\n")
				cat("This means either\n(a) you have a conflict between the length of states_list\nand the dimensions of the tip likelihoods, Qmat, etc.; or\n(b) it means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\ndecided they were effectively 0; default min_precision=1e-16.\n(c) include_null_range setting is wrong, you have include_null_range=", include_null_range, "\n...Probably...(a) or (c).", sep="")
				cat("\n")
				} # END stop_function
			stop(stop_function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat))
			}
		} else {
		numstates_in_tip_condlikes_of_data_on_each_state = ncol(tip_condlikes_of_data_on_each_state)
		
		# Number of states in the spPmat, *including* the null range if desired
		if (include_null_range == TRUE)
			{
			numstates_in_spPmat_ie_nrows = 1 + nrow(spPmat)
			} else {
			numstates_in_spPmat_ie_nrows = 0 + nrow(spPmat)
			} # END if (include_null_range == TRUE)

		if (is.null(m) == FALSE)
			{
			numstates_in_spPmat_ie_nrows = numstates_in_spPmat_ie_nrows * length(m)
			}

		
		if (input_is_COO == TRUE)
			{
			numstates_in_Qmat = max( max(Qmat[,1]), max(Qmat[,2]) )
			} else {
			numstates_in_Qmat = ncol(Qmat)
			} # END if (input_is_COO == TRUE)
		
		if ( all(numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_spPmat_ie_nrows, numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_Qmat ) == TRUE )
			{
			# Continue
			all_inputs_correct_size = TRUE
			} else {
			
			# Stop if not everything is equal
			stop_function2 <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat)
				{
				cat("\n\nERROR: Some inputs have incorrect size -- \n")
				cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
				cat("numstates_in_spPmat_ie_nrows+1:				", numstates_in_spPmat_ie_nrows, "\n")
				cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
				cat("\n")
				cat("This means either\n(a) you have a conflict between the length of states_list\nand the dimensions of the tip likelihoods, Qmat, etc.; or\n(b) it means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\ndecided they were effectively 0; default min_precision=1e-16.\n(c) include_null_range setting is wrong, you have include_null_range=", include_null_range, "\n...Probably...(a) or (c).", sep="")
				cat("\n")
				} # END stop_function2
			stop(stop_function2(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat))
			}
		} # END if (use_cpp == TRUE && !is.null(spPmat_inputs))
		# END error check
	
	
	
	
	
	
	# Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	if ( is.null(spPmat_inputs)==FALSE )
		{
		spPmat_inputs$l[spPmat_inputs$l == c("_")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-1")] = NULL
		#spPmat_inputs$l[spPmat_inputs$l == c(-1)] = NULL
		}
	
	
	
	# Calculate likelihoods down tree
	#numstates = nrow(Qmat)
	numstates = ncol(tip_condlikes_of_data_on_each_state)

	
	#######################################################
	# Check if the phylogeny is actually just a number (i.e., a branch length)
	#######################################################
	# Hmm, this seems harder...


	#######################################################
	# The rest assumes that you've got a phylo3 object, in default order
	#######################################################
	
	
	
	
	edgelengths = phy$edge.length
	num_branches_below_min = sum(edgelengths < min_branchlength)

	if ( (printlevel >= 1) && (num_branches_below_min > 0) )
		{
		cat("\n")
		cat("Running calc_loglike_sp():\n")
		cat("This run of calc_loglike_sp() has a min_branchlength of: ", min_branchlength, "\n", sep="")
		cat("Branches shorter than this will be assumed to be connected to the tree with\n")
		cat("sympatric events (i.e., members of fossil lineages on ~0 length branches.)\n")
		cat("This tree has ", num_branches_below_min, " branches < ", min_branchlength, ".\n", sep="")
		cat("\n")
		}	

	
	num_internal_nodes = phy$Nnode
	numtips = length(phy$tip.label)
	
	# likelihoods are computed at all nodes
	# make a list to store the 
	numnodes = numtips + num_internal_nodes
	computed_likelihoods_at_each_node = numeric(length=numnodes)
	
	
	
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	phy2 <- reorder(phy, "pruningwise")
	
	
	tipnums <- 1:numtips

	# Put in the sums of the probabilities of the states at each tip
	# (This only works if the tip data are 0000100000, etc...)
	#computed_likelihoods_at_each_node[tipnums] = rowSums(tip_condlikes_of_data_on_each_state)
	
	# 2014-05-07_NJM: change to divide by the sum
	computed_likelihoods_at_each_node[tipnums] = rowSums(tip_condlikes_of_data_on_each_state) / rowSums(tip_condlikes_of_data_on_each_state)
	
	#######################################################
	# Initialize matrices for downpass and uppass storage of state probabilities
	#######################################################
	
	# This is all you need for a standard likelihood calculation
	# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = rel probs AT A NODE
	# BE SURE TO DISTINGUISH STATE PROBABILITIES (VS LIKELIHOODS, WHICH ARE NOT NORMALIZED)
	condlikes_of_each_state <- matrix(data=0, nrow=numnodes, ncol=numstates)
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS <- matrix(data=0, nrow=numnodes, ncol=numstates)
	
	# In THEORY, this should be OK to divide, but BECAUSE I pass relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS to 
	# the matrix exponentiation calculations, these need to be likelihoods at the tips;
	# SO DON'T DIVIDE, if you do it normalizes the probabilities at the tips, and this screws up things that are being passed 
	# likelihoods that are not 1/0000, e.g. stratification analyses, and causes M0 and M0 strat to disagree on Psychotria (e.g.)
	#
	# At some point we should more rigorously separate relprobs and condlikes, but don't screw with it now!
	# 
	# For e.g. uppass calculations of tip probs, just do a hack later.
	# 
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[tipnums, ] = tip_condlikes_of_data_on_each_state #/ rowSums(tip_condlikes_of_data_on_each_state)
	
	# BUT, DO NOT DIVIDE THIS BY rowSums(tip_condlikes_of_data_on_each_state), it forces normalization which prevents e.g.
	# passing down likelihoods during stratification
	condlikes_of_each_state[tipnums, ] = tip_condlikes_of_data_on_each_state
	#relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[tipnums, ] <- 1
	
	
	# But, if you want to do ancestral states properly, and get marginal
	# reconstructions, you've gotta store:
	# 1. Probability of the states at the bottom of each branch on downpass
	# 2. Probability of the states at the bottom of each branch on UP-pass
	# 3. Probability of the states at the top of each branch on UP-pass
	# Combine 0-3 to get the ML marginal probability of states at
	# 4. The bottom of each branch
	# 5. The top of each branch
	# Combine 4 & 5 to get:
	# 6. MAYBE The ML marginal probability of each split scenario at each node - i.e. #5(top) * #4(left bot) * #4(right bot)
	# 
	if (calc_ancprobs == TRUE)
		{
		# Every node (except maybe the root) has a branch below it, and there is also a 
		# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS at the bottom of this branch
		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)

		ML_marginal_prob_each_state_at_branch_bottom_below_node <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		ML_marginal_prob_each_state_at_branch_top_AT_node <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		ML_marginal_prob_each_split_at_branch_top_AT_node = list()
		}



	
	# If you want to use standard dense matrix exponentiation, rexpokit's dgpadm is best,
	# and you can run it through mapply on the branch lengths of all branches at once
	if (sparse==FALSE)
		{
		# Get probmats for each branch, put into a big array
		# Creates empty array to store results
		independent_likelihoods_on_each_branch = calc_independent_likelihoods_on_each_branch(phy2, Qmat, cluster_already_open, Qmat_is_sparse=FALSE)
		} else {
		# Sparse matrices
		# Here, if you want to use sparse matrix exponentiation, you are NOT
		# going to store the entire Pmat for each branch; you are going to directly
		# calculate the output probabilities w (w), via w=exp(Qt)v via myDMEXPV.
		#
		# This is most efficiently done if you transpose and COO-ify the matrix here,
		# ahead of time.
		if (input_is_COO==FALSE)
			{
			original_Qmat = Qmat
			
			# number of states in the original matrix
			coo_n = numstates
			anorm = as.numeric(norm(original_Qmat, type="O"))
			matvec = original_Qmat
			
			# DO NOT TRANSPOSE; we want to go BACKWARDS in time, NOT FORWARDS!
			#tmatvec = base::t(matvec)
			tmatvec = matvec
			tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
			} else {
			# The input Qmat is already in COO format (untransposed, as we are doing 
			# likelihoods backwards in time)
			
			# CHECK THAT ITS IN COO FORMAT
			if ( (class(Qmat) != "data.frame") || (ncol(Qmat) != 3) )
				{
				stop("ERROR: calc_loglike_sp is attempting to use a sparse COO-formated Q matrix, but you provided a regular square dense matrix")
				}
			coo_n = numstates
			anorm = 1
			tmpQmat_in_REXPOKIT_coo_fmt = Qmat
			#tmpQmat_in_REXPOKIT_coo_fmt = cooQmat
			} # END if (input_is_COO==FALSE)
		} # END if (sparse==FALSE)


	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)





use_cpp = TRUE
cppSpMethod = 3
include_null_range = BioGeoBEARS_run_object$include_null_range
include_null_range

	# Calculate the rowsums of the speciation matrix ONCE
	if ( use_cpp == TRUE )
		{
		if ( is.null(spPmat_inputs)==FALSE )
			{
			# Calculate the rowsums (for input into rcpp_calc_anclikes_sp()
			# Actually, just do this ONCE

			# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
			# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
			l = spPmat_inputs$l		# states_indices
			s = spPmat_inputs$s
			v = spPmat_inputs$v
			j = spPmat_inputs$j
			y = spPmat_inputs$y
			
			# dmat means combined dispersal multipliers matrix, NOT dmat_times_d
			dmat = spPmat_inputs$dmat
			
			# Take the max of the indices of the possible areas, and add 1
			# numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1 # old, bogus
			numareas = max(sapply(X=spPmat_inputs$l, FUN=length), na.rm=TRUE) + 0
			
			maxent01s_param = spPmat_inputs$maxent01s_param
			maxent01v_param = spPmat_inputs$maxent01v_param
			maxent01j_param = spPmat_inputs$maxent01j_param
			maxent01y_param = spPmat_inputs$maxent01y_param
			
			maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
			maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
			maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
			maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

			# You really need a list of sizes here:
			
			# Matrix of probs for each ancsize
			maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
			maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
			
			# Now, go through, and make a list of the max minsize for each decsize
			max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)

			# -1 for null range
			if (include_null_range == TRUE)
				{
				state_space_size_Qmat_to_cladoMat = -1
				} else {
				state_space_size_Qmat_to_cladoMat = 0
				}
			
			tmpca_1 = rep(1, (ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS) + state_space_size_Qmat_to_cladoMat))
			tmpcb_1 = rep(1, (ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS) + state_space_size_Qmat_to_cladoMat))

			# Print the matrix to screen from C++
			printmat = FALSE
			if (printlevel >= 2)
				{
				printmat = TRUE
				}
			
			# Calculate the rowSums of the speciation matrix, for future reference (i.e., setting all input likelihoods to 1)
			# Only need be done once.
			#
			# But, actually, what would make this all REALLY efficient would be to just calculate the conditional
			# probability matrix ONCE, storing it in a COO-like format.  The Rsp_rowsums would be easily derived
			# from that, and we wouldn't have to calculate the speciation model Nnodes times independently.
			#Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)
			
			# Get the speciation matrix conditional probabilities in a COO-like format
			# [[1]] = inums = indexes of left descendant state in speciation matrix, by ancestral rowsnums 1-15
			# [[2]] = jnums = indexes of right descendant state in speciation matrix, by ancestral rowsnums 1-15
			# [[3]] = probs = probvals of this left|right combination in speciation matrix, by ancestral rowsnums 1-15
			if (cppSpMethod == 2)
				{
				# Error check
				if (is.null(m) == FALSE)
					{
					txt = "STOP ERROR in calc_loglike_sp_prebyte(): if m is not NULL, only option cppSpMethod==3 can be used. You have cppSpMethod==2."
					cat("\n\n")
					cat(txt)
					cat("\n\n")
					stop(txt)
					} # END if (is.null(m) == FALSE)
				
				COO_probs_list = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)

				# Sum the probabilities (list [[3]]) in each row ([[3]][[1]] list of probs
				# through [[3]][[15]] list of probs)
				Rsp_rowsums = sapply(X=COO_probs_list[[3]], FUN=sum)
				} # END if (cppSpMethod == 2)
	
			if (cppSpMethod == 3)
				{
				if (printlevel >= 2)
					{
					params_to_print = c("tmpca_1", "tmpcb_1", "l", "s", "v", "j", "y", "dmat", "maxent01s", "maxent01v", "maxent01j", "maxent01y", "max_minsize_as_function_of_ancsize")

					for (tmppval in params_to_print)
						{
						cmdstr = paste(	"cat('", tmppval, "', ':\n', sep='')", sep="")
						eval(parse(text=cmdstr))
						
						# Get the value
						cmdstr = paste("tmppval = ", tmppval, sep="")
						eval(parse(text=cmdstr))
						
						# Print it
						print(tmppval)
						}
					} # END if (printlevel >= 2)
				
				#print(m)
				
				COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat, m=m, m_null_range=include_null_range, jts_matrix=jts_matrix)
				
				# combine with C++ function
				# This causes an error with spPmat=NULL; spPmat_inputs=NULL; use_cpp=TRUE; sparse=FALSE
				# i.e. gives 16 states with a 0 on the end, rather than 15 states
				#Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates)
				
				# This gives 15 states
				Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)
				} # END if (cppSpMethod == 3)
			} # END if ( is.null(spPmat_inputs)==FALSE )
		} # END if ( use_cpp == TRUE )
	


round(COO_weights_columnar[[4]], 4)
cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], COO_weights_columnar[[4]])





	min_branchlength=0.000001
	return_what = "all"
	probs_of_states_at_root=NULL
	rootedge=FALSE
	cppSpMethod=3
	cluster_already_open=NULL
	printlevel=1
	input_is_COO=FALSE
	
	stratified=FALSE
	use_cpp=TRUE
	#sparse=TRUE
	sparse=FALSE
	
	#include_null_range=TRUE
	include_null_range=TRUE
	m=NULL
	
	states_allowed_TF=NULL
	probs_of_states_at_root = NULL




	
	#######################################################
	#######################################################
	# THIS IS A DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	#######################################################
	for (i in edges_to_visit)
		{
		# First edge visited is i
		#print(i)
		
		# Its sister is j 
		j <- i + 1
		#print(j)

		# Get the node numbers at the tips of these two edges		
		left_desc_nodenum <- phy2$edge[i, 2]
		right_desc_nodenum <- phy2$edge[j, 2]
		
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		anc <- phy2$edge[i, 1]
		
		txt = paste("anc:", anc, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# Is sparse is FALSE, input the pre-calculated likelihoods;
		# If sparse is TRUE, dynamically calculate using expokit_dmexpv_Qmat
		if (sparse==FALSE)
			{
			
			#print("Checking Q matrix")
			#cat("p=", p, ", rate=", rate, "\n", sep=" ")
			#print(Q)
			#print(phy$edge.length[i])
			#condlikes_Left <- matexpo(Qmat * phy2$edge.length[i]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
			#condlikes_Right <- matexpo(Qmat * phy2$edge.length[j]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]


			if (printlevel >= 2)
				{
				print("dense matrix exponentiation")
				}
			
			
			if (is.null(cluster_already_open))
				{
				# Conditional likelihoods of states at the bottom of left branch
				condlikes_Left = independent_likelihoods_on_each_branch[,,i] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
							
				# Conditional likelihoods of states at the bottom of right branch
				condlikes_Right = independent_likelihoods_on_each_branch[,,j] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]
				} else {
				
				#cat("dim(independent_likelihoods_on_each_branch[[i]]):\n", sep="")
				#cat(dim(independent_likelihoods_on_each_branch))
				#cat("\n\n")
				
				#cat("length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]):\n", sep="")
				#cat(length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]))
				#cat("\n\n")
				
				condlikes_Left = independent_likelihoods_on_each_branch[[i]] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
				
				# Conditional likelihoods of states at the bottom of right branch
				condlikes_Right = independent_likelihoods_on_each_branch[[j]] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]

				}


			# Error check
			if (any(is.nan(condlikes_Left)))
				{
				cat("\n\nindependent likelihood %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS returned NaNs:\n\n")

				cat("\n\nPrinting parameter values being attempted when this occurred:\n\n")
				if (is.null(BioGeoBEARS_model_object))
					{
					print("BioGeoBEARS_model_object was not passed to calc_loglike_sp() for printing.")
					} else {
					print(BioGeoBEARS_model_object)
					} # END if (is.null(BioGeoBEARS_model_object))
				
				cat("i=", i, "\n")
				cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
				cat("\n\nprint(condlikes_Left):\n\n")
				print(condlikes_Left)

				cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS):\n\n")
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)

				cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]):\n\n")
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,])
				#print(coo_n)
				#print(anorm)
				
				if (is.null(on_NaN_error))
					{
					stop("\n\nStopping on error: NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood,, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n")
					}
				
				if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					warning(paste0("\n\nWarning on error: NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
					
					return(on_NaN_error)
					}
				
				}

			if (any(is.nan(condlikes_Right)))
				{
				cat("\n\nindependent likelihood %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS returned NaNs:\n\n")
				
				cat("\n\nPrinting parameter values being attempted when this occurred:\n\n")
				if (is.null(BioGeoBEARS_model_object))
					{
					print("BioGeoBEARS_model_object was not passed to calc_loglike_sp() for printing.")
					} else {
					print(BioGeoBEARS_model_object)
					} # END if (is.null(BioGeoBEARS_model_object))
				
				cat("j=", j, "\n")
				cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
				cat("\n\nprint(condlikes_Right):\n\n")
				print(condlikes_Right)

				cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS):\n\n")
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)

				cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]):\n\n")
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,])
				#print(coo_n)
				#print(anorm)

				if (is.null(on_NaN_error))
					{
					stop("\n\nStopping on error: NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood,, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n")
					}
				
				print("print(on_NaN_error):")
				print(on_NaN_error)
				if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					warning(paste0("\n\nWarning on error: NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
					
					return(on_NaN_error)
					}
				} # END if (any(is.nan(condlikes_Right)))
			
			if (printlevel >= 2) {
			txt = paste("condlikes at bottom of L: ", paste(round(condlikes_Left, 4), collapse=" ", sep=""), sep="")
			print(txt)
			}
			
			if (printlevel >= 2) {
			txt = paste("condlikes at bottom of R: ", paste(round(condlikes_Right, 4), collapse=" ", sep=""), sep="")
			print(txt)
			}
			#condlikes_Left <- expm(Qmat * phy$edge.length[i], method="Ward77") %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum, 
			#  ]
			#condlikes_Right <- expm(Qmat * phy$edge.length[j], method="Ward77") %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum, 
			#  ]
			
			
			#condlikes_Left <- exp(Qmat * phy$edge.length[i]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum, 
			#  ]
			#condlikes_Right <- exp(Qmat * phy$edge.length[j]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum, 
			#  ]
			} else {
			#######################################################
			# Sparse matrix exponentiation
			#######################################################
			# sparse == TRUE

			# For rapid exponentiation of sparse matrices, we use myDMEXPV, and input
			# starting probabilities/likelihoods, and output ending probabilities/likelihoods
			#
			# When there ARE inputprobs_for_fast, myDMEXPV is invoked, which appears to do a 
			# forward-probability matrix exponentiation calculation, not a backward probability
			# matrix exponentiation calculation.  This produces a positive value for a state which is 
			# impossible in the ancestor (e.g. the null range is impossible in an ancestor)
			#
			# To fix this, multiple the output probabilities by 0 if check_for_0_rows[i[ == TRUE

			if (printlevel >= 2)
				{
				print("sparse matrix exponentiation")
				}
			
			# Conditional likelihoods of data given the states at the bottom of left branch
			condlikes_Left = try (
			expokit_dmexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[i], inputprobs_for_fast=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,], transpose_needed=FALSE, transform_to_coo_TF=FALSE, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
			)
			
			# Error check
			if (class(condlikes_Left) == "try-error")
				{
				cat("\n\ntry-error on expokit_dmexpv_Qmat():\n\n")
				cat("i=", i, "\n")
				cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}
			
			if (any(is.nan(condlikes_Left)))
				{
				cat("\n\nexpokit_dmexpv_Qmat() returned NaNs:\n\n")
				cat("i=", i, "\n")
				cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}
	
			# Conditional likelihoods of data given the states at the bottom of right branch
			condlikes_Right = try(
			expokit_dmexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[j], inputprobs_for_fast=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,], transpose_needed=FALSE, transform_to_coo_TF=FALSE, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
			)

			# Error check
			if (class(condlikes_Right) == "try-error")
				{
				cat("\n\ntry-error on expokit_dmexpv_Qmat():\n\n")
				cat("j=", j, "\n")
				cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}


			if (any(is.nan(condlikes_Right)))
				{
				cat("\n\nexpokit_dmexpv_Qmat() returned NaNs:\n\n")
				cat("j=", j, "\n")
				cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
				print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}

			#print(c(condlikes_Left))
			#print(c(condlikes_Right))
			} # END if (sparse==FALSE) (ENDING MATRIX EXPONENTIATION)	
		
		# Zero out impossible states
		if (!is.null(states_allowed_TF))
			{
			condlikes_Left[states_allowed_TF==FALSE] = 0
			condlikes_Right[states_allowed_TF==FALSE] = 0
			}
	
		
		
		# Save the conditional likelihoods of the data at the bottoms of each branch 
		# (needed for marginal ancestral state probabilities)
		if (calc_ancprobs == TRUE)
			{
			# Every node (except maybe the root) has a branch below it, and there is also a 
			# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS at the bottom of this branch
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,] = condlikes_Left / sum(condlikes_Left)
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,] = condlikes_Right / sum(condlikes_Right)
			}			
			
		
		
		
		# If there is no speciational model, you are assuming 100% sympatry (range duplication)
		# at each speciation event
		#
		# In this case, you can just multiply the two conditional likelihood matrices together
		#
		# Also, if a branch is extremely short (a "hook"), this is essentially a zero-length
		# branch, we are assuming that this represents the range of a lineage at that 
		# point.  There is no speciation event here -- both "lineages" inherit
		# the same range.  This allows fossils to closely influence ancestral states.
		#
		# This was developed with Kaitlin Maguire over several years of screwing around.
		
		# Check for a short "hook" branch; if found, use just allopatric speciational model

		# get the correct edge
		left_edge_TF = phy2$edge[,2] == left_desc_nodenum
		right_edge_TF = phy2$edge[,2] == right_desc_nodenum
		
		# Check the branchlength of each edge
		is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
		is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
		hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0
		
# 		print(left_desc_nodenum)
# 		print(right_desc_nodenum)
# 		
# 		print(phy2$edge.length[left_desc_nodenum])
# 		print(phy2$edge.length[right_desc_nodenum])
# 		
# 		print(is_leftbranch_hook_TF)
# 		print(is_rightbranch_hook_TF)
# 		
# 		print(hooknode_TF)
		if (use_cpp == FALSE)
			{
			if ( (is.null(spPmat) || hooknode_TF==TRUE) )
				{
				##################################
				# no speciational model of range change
				##################################
				node_likelihood <- condlikes_Left * condlikes_Right
				if (printlevel >= 2)
					{
					print("use_cpp=FALSE, direct multiplication")
					}
				} else {
				##################################
				# combine the likelihoods from each branch bottom with this speciational model
				# of range change
				##################################

				if (printlevel >= 2)
					{
					print("use_cpp=FALSE, speciation model")
					}
				
				# for each ancestral state, get prob of branch pairs 
				outmat = matrix(0, nrow=nrow(spPmat), ncol=ncol(spPmat))
				
				# Get the probs of each pair, list by row
				# (This might be a slightly slow step for large matrices)
				
				# Exclude "_" ranges from this calculation, as if there is speciation,
				# they are not going extinct
				
				if (include_null_range == TRUE)
					{
					probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left[-1]), c(condlikes_Right[-1]), "*")))
					} else {
					probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left), c(condlikes_Right), "*")))
					}
				#sum(probs_of_each_desc_pair)
	
				# Make a matrix with the probability of each pair, in each cell
				# (duplicated for each ancestral state)
				for (spi in 1:nrow(spPmat))
					{
					outmat[spi,] = probs_of_each_desc_pair
					}
				# Multiply the probabilities of each pair, by the probability of each
				# type of speciation event, to get the probabilities at the 
				# ancestral node
				outmat2 = outmat * spPmat
				
				node_likelihood_with_speciation = rowSums(outmat2)
				
				
				# THIS ZERO IS ALREADY OVER-WRITING THE NULL STATE LIKELIHOOD!!
				
				# Add the 0 back in, representing 0 probability of "_"
				# range just below speciation event
				if (include_null_range == TRUE)
					{
					node_likelihood = c(0, node_likelihood_with_speciation)
					} else {
					node_likelihood = node_likelihood_with_speciation
					}
				
				print("node_likelihood:")
				print(node_likelihood)
				
				#######################################################
				# If the states/likelihoods have been fixed at a particular node
				#######################################################
				if (!is.null(fixnode))
					{
					# For multiple fixnodes
					if (length(fixnode) > 1)
						{
						# Get the matching node
						TF = (anc == fixnode)
						temporary_fixnode = fixnode[TF]
						temporary_fixlikes = c(fixlikes[TF,])
						} else {
						temporary_fixnode = fixnode
						temporary_fixlikes = c(fixlikes)
						}

					
					if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
						{
						# If the node is fixed, ignore the calculation for this node, and
						# instead use the fixed likelihoods (i.e., the "known" state) for
						# this node.
						# fix the likelihoods of the (NON-NULL) states
						node_likelihood = node_likelihood * temporary_fixlikes
						}
					} # end if (!is.null(fixnode))
				}
			} else {
			##################################
			# use_cpp == TRUE
			##################################
			if ( hooknode_TF==TRUE )
				{
				##################################
				# no speciational model of range change
				##################################
				if (printlevel >= 2) 
					{
					print("use_cpp=TRUE, direct multiplication")
					}
# 				print("use_cpp=TRUE, direct multiplication")
# 				print("condlikes_Left:")
# 				print(condlikes_Left)
# 				print("condlikes_Right:")
# 				print(condlikes_Right)
				node_likelihood <- condlikes_Left * condlikes_Right			
				} else {
				if ( is.null(spPmat_inputs)==TRUE )
					{
					if ( is.null(spPmat) == TRUE )
						{
						##################################
						# no speciational model of range change
						##################################
						node_likelihood <- condlikes_Left * condlikes_Right
						print ("WARNING #2: (use_cpp==TRUE && is.null(spPmat_inputs)==TRUE)")
						print("use_cpp=TRUE, no spPmat_inputs, direct multiplication")

						} else {
						# otherwise ...
						if (printlevel >= 2)
							{
							print("use_cpp=TRUE, no spPmat_inputs, speciation model")
							}
						##################################
						# combine the likelihoods from each branch bottom with this speciational model
						# of range change
						##################################

						# for each ancestral state, get prob of branch pairs 
						outmat = matrix(0, nrow=nrow(spPmat), ncol=ncol(spPmat))
						
						# Get the probs of each pair, list by row
						# (This might be a slightly slow step for large matrices)
						
						# Exclude "_" ranges from this calculation, as if there is speciation,
						# they are not going extinct
						
						if (include_null_range == TRUE)
							{
							probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left[-1]), c(condlikes_Right[-1]), "*")))
							} else {
							probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left), c(condlikes_Right), "*")))
							}
						#sum(probs_of_each_desc_pair)
			
						# Make a matrix with the probability of each pair, in each cell
						# (duplicated for each ancestral state)
						for (spi in 1:nrow(spPmat))
							{
							outmat[spi,] = probs_of_each_desc_pair
							}
						# Multiply the probabilities of each pair, by the probability of each
						# type of speciation event, to get the probabilities at the 
						# ancestral node
						outmat2 = outmat * spPmat
						
						node_likelihood_with_speciation = rowSums(outmat2)
						
						# THIS ZERO IS ALREADY OVER-WRITING THE NULL STATE LIKELIHOOD!!
						
						# Add the 0 back in, representing 0 probability of "_"
						# range just below speciation event
						if (include_null_range == TRUE)
							{
							node_likelihood = c(0, node_likelihood_with_speciation)
							} else {
							node_likelihood = node_likelihood_with_speciation
							}

						print("node_likelihood:")
						print(node_likelihood)


						#######################################################
						# If the states/likelihood have been fixed at a particular node
						#######################################################
						if (!is.null(fixnode))
							{
							# For multiple fixnodes
							if (length(fixnode) > 1)
								{
								# Get the matching node
								TF = (anc == fixnode)
								temporary_fixnode = fixnode[TF]
								temporary_fixlikes = c(fixlikes[TF,])
								} else {
								temporary_fixnode = fixnode
								temporary_fixlikes = c(fixlikes)
								}
					
							if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
								{
								# If the node is fixed, ignore the calculation for this node, and
								# instead use the fixed likelihoods (i.e., the "known" state) for
								# this node.
								# fix the likelihoods of the (NON-NULL) states
								node_likelihood = node_likelihood * temporary_fixlikes
								}
							}
						}
					} else {
					## if ( is.null(spPmat_inputs)==FALSE )
					if (printlevel >= 2)
						{
						print("use_cpp=TRUE, yes spPmat_inputs, speciation model")
						}
					################################
					# Use the C++ function!
					################################
					##################################
					# combine the likelihoods from each branch bottom with this speciational model
					# of range change
					##################################
					
					# Calculate the likelihoods at each node,
					# give the speciation model, and input probs for 
					# each branch
					if (include_null_range == TRUE)
						{
						ca = condlikes_Left[-1]
						cb = condlikes_Right[-1]
						} else {
						ca = condlikes_Left
						cb = condlikes_Right						
						}
					
					# Rcpp does weird alterations to the input variables, so use each only once!
					#tmpca_1 = ca
					#tmpcb_1 = cb
					tmpca_2 = ca
					tmpcb_2 = cb
					
					# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
					# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
					l = spPmat_inputs$l		# states_indices
					s = spPmat_inputs$s
					v = spPmat_inputs$v
					j = spPmat_inputs$j
					y = spPmat_inputs$y
					
					# dmat means combined dispersal multipliers matrix, NOT d
					dmat = spPmat_inputs$dmat
	
					# Print the matrix to screen from C++
					printmat = FALSE
					if (printlevel >= 2) {
					printmat = TRUE
					}
					
					# Calculate the rowsums (for input into rcpp_calc_anclikes_sp()
					# Actually, just do this ONCE (above)
					#Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s, maxent01v, maxent01j, maxent01y, printmat=printmat)
					
					if (printlevel >= 2)
						{
						print("Rsp_rowsums:")
						print(Rsp_rowsums)
						}
					
					#node_likelihood_with_speciation = rep(0.0, length(tmpca_2))
					# This calculates the speciation probabilities again at each node; this is inefficient
					if (cppSpMethod == 1)
						{
						node_likelihood_with_speciation = rcpp_calc_anclikes_sp(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, Rsp_rowsums=Rsp_rowsums, printmat=printmat)						
						#print(node_likelihood_with_speciation[1:5])
						} # END if (cppSpMethod == 1)
					
					# Really, we should just iterate through the COO_probs_list using a C++ function
					# 2012-12-04 NJM says: this KICKS ASS in C++!!
					
					if (cppSpMethod == 2)
						{
						node_likelihood_with_speciation2 = rep(0.0, length(tmpca_2))
						node_likelihood_with_speciation2 = rcpp_calc_anclikes_sp_using_COOprobs(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, RCOO_left_i_list=COO_probs_list[[1]], RCOO_right_j_list=COO_probs_list[[2]], RCOO_probs_list=COO_probs_list[[3]], Rsp_rowsums=Rsp_rowsums, printmat=printmat)
						#print(node_likelihood_with_speciation2[1:5])
						node_likelihood_with_speciation = node_likelihood_with_speciation2
						} # END if (cppSpMethod == 2)

					if (cppSpMethod == 3)
						{
						node_likelihood_with_speciation3 = rep(0.0, length(tmpca_2))
						node_likelihood_with_speciation3 = rcpp_calc_splitlikes_using_COOweights_columnar(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, COO_weights_columnar=COO_weights_columnar, Rsp_rowsums=Rsp_rowsums, printmat=printmat)
						#print(node_likelihood_with_speciation2[1:5])
						node_likelihood_with_speciation = node_likelihood_with_speciation3
						
						# Null ranges can leave NaNs
						# when the downpass is calculated on 
						# ranges x trait m's
						if (length(m) > 0)
							{
							nanTF = is.nan(node_likelihood_with_speciation)
							
							# IF this is the issue, the number of NaNs will be length(m)-1
							# (don't want to zero out any weird NaNs)
							if (sum(nanTF) == (length(m)-1))
								{
								node_likelihood_with_speciation[nanTF] = 0.0
								}
							}
						} # END if (cppSpMethod == 3)


					

					
					if (printlevel >= 2)
						{
						print("node_likelihood_with_speciation:")
						print(node_likelihood_with_speciation)
						}
	
					if (include_null_range == TRUE)
						{
						node_likelihood = c(0, node_likelihood_with_speciation)
						} else {
						node_likelihood = node_likelihood_with_speciation
						}



					#######################################################
					# If the states/likelihood have been fixed at a particular node
					#######################################################
					if (!is.null(fixnode))
						{
						# For multiple fixnodes
						if (length(fixnode) > 1)
							{
							# Get the matching node
							TF = (anc == fixnode)
							temporary_fixnode = fixnode[TF]
							temporary_fixlikes = c(fixlikes[TF,])
							} else {
							temporary_fixnode = fixnode
							temporary_fixlikes = c(fixlikes)
							}
					
						if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
							{
							# If the node is fixed, ignore the calculation for this node, and
							# instead use the fixed likelihoods (i.e., the "known" state) for
							# this node.
							# fix the likelihoods of the (NON-NULL) states
							node_likelihood = node_likelihood * temporary_fixlikes
							}
						} # END if (!is.null(fixnode))

					#formatted_table = conditional_format_table(t(node_likelihood))
					#print(formatted_table)

					#print("!is.null(states_allowed_TF)")
					#print(!is.null(states_allowed_TF))
					# Zero out impossible states					
					if (!is.null(states_allowed_TF))
						{
						
						node_likelihood[states_allowed_TF==FALSE] = 0
						} # END if (!is.null(states_allowed_TF))
					} # END if ( is.null(spPmat_inputs)==TRUE )
				} # END if ( hooknode_TF==TRUE )
			} # END if (use_cpp == FALSE)


		########################################################################
		# If you are at the root node, you also multiply by the probabilities
		# of the starting states
		########################################################################
		if (i == max(edges_to_visit))
			{
			# If not specified, you are assuming even probabilities of each state
			if (is.null(probs_of_states_at_root) == TRUE)
				{
				node_likelihood = node_likelihood
				} else {
				# Otherwise, use user-specified probs_of_states_at_root
				# These could be:
				# 1. User-specified base frequencies
				# 2. Observed frequencies at tips
				#    (WARNING: likely have states with 0 tip observations!)
				# 3. Observed frequencies in some reference set of species
				# 4. Some distribution on range sizes, even within range sizes
				# 5. Equilibrium probabilities assuming Qmat rate matrix is
				#    run to infinite time
				node_likelihood = probs_of_states_at_root * node_likelihood 
				}
			}

		
		#txt = paste("left node: ", left_desc_nodenum, ", right node:", right_desc_nodenum, sep="")
		#print(txt)
		#print(node_likelihood)
		
		nanTF = is.nan(node_likelihood)
		if (sum(nanTF) > 0)
			{
			anc <- phy2$edge[i, 1]

			errortxt = paste("ERROR in calc_loglike_sp at 'if (sum(nanTF) > 0)': you have ", sum(nanTF), "NaNs in the 'node_likelihood' at edges #", i, " & ", j, ", node #", anc, "!", sep="")
			cat("\n\n")
			cat(errortxt)
			cat("\n\n")
			
			cat("\n\nnanTF:\n\n")
			print(nanTF)
			cat("\n\nprint(node_likelihood):\n\n")
			print(node_likelihood)
			print(i)
			print(j)
			left_desc_nodenum
			right_desc_nodenum
		
			# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
			#anc <- phy2$edge[i, 1]
			#plot(phy2)
			#nodelabels()
			
			#printall(prt(phy2))
			
			#print(tip_condlikes_of_data_on_each_state[left_desc_nodenum,])
			#print(tip_condlikes_of_data_on_each_state[right_desc_nodenum,])
			
			stop()
			}
		
		total_likelihood_for_node = sum(node_likelihood)
		
		#print(total_likelihood_for_node)
		
		computed_likelihoods_at_each_node[anc] = total_likelihood_for_node
		
		relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ] = node_likelihood / total_likelihood_for_node
		condlikes_of_each_state[anc, ] = node_likelihood
		
		#print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
		} # end downpass
	#######################################################
	# END DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	rootnode = anc




node_likelihood_with_speciation3 = rep(0.0, length(tmpca_2))
node_likelihood_with_speciation3 = rcpp_calc_splitlikes_using_COOweights_columnar(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, COO_weights_columnar=COO_weights_columnar, Rsp_rowsums=Rsp_rowsums, printmat=printmat)

round(Rsp_rowsums,6)

round(COO_weights_columnar[[4]], 4)
cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], COO_weights_columnar[[4]])



#print(node_likelihood_with_speciation2[1:5])
node_likelihood_with_speciation = node_likelihood_with_speciation3



# Behavior of:
# rcpp_calc_splitlikes_using_COOweights_columnar

tmpca_2 = c(1, 5e-25, 1e-12)
tmpcb_2 = c(5e-25, 1, 1e-12)

COO_weights_columnar = list(c(0L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 0L, 1L, 0L, 1L), c(0L, 
1L, 2L, 2L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L), c(0L, 1L, 0L, 1L, 
2L, 2L, 0L, 1L, 1L, 0L, 0L, 1L), c(3.33333332491748e-06, 3.33333332491748e-06, 
3.33333332491748e-06, 3.33333332491748e-06, 3.33333332491748e-06, 
3.33333332491748e-06, 3.33333332491748e-06, 3.33333332491748e-06, 
2.99998998641968, 2.99998998641968, 2.99998998641968, 2.99998998641968
))

printmat = FALSE

Rsp_rowsums = c(5.999983, 5.999983, 0.000020)
#Rsp_rowsums = c(6, 6, 6)

Rcpp_leftprobs = tmpca_2
Rcpp_rightprobs = tmpcb_2


	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
	RCOO_left_i_list = COO_weights_columnar[[2]]
	RCOO_right_j_list = COO_weights_columnar[[3]]
	RCOO_probs_list = COO_weights_columnar[[4]]
	
	# Call the fast C++ function
	splitlikes = .Call( "cpp_calc_splitlikes_using_COOweights_columnar", leftprobs=as.numeric(Rcpp_leftprobs), rightprobs=as.numeric(Rcpp_rightprobs), RCOO_weights_columnar_anc_i_list=as.integer(RCOO_weights_columnar_anc_i_list), RCOO_left_i_list=as.integer(RCOO_left_i_list), RCOO_right_j_list=as.integer(RCOO_right_j_list), RCOO_probs_list=as.numeric(RCOO_probs_list), Rsp_rowsums=as.numeric(Rsp_rowsums), PACKAGE = "cladoRcpp" )
splitlikes


# R version of the above calculations:
# 
num_splits_scenarios = length(COO_weights_columnar[[1]])
num_ancestral_ranges = length(Rsp_rowsums)
tmp_split_likes = rep(0, times=num_ancestral_ranges)
cat("\n")
for (i in 1:num_splits_scenarios)
	{
	tmp_probval = Rcpp_leftprobs[1+COO_weights_columnar[[2]][i]] * Rcpp_rightprobs[1+COO_weights_columnar[[3]][i]] * COO_weights_columnar[[4]][i]
	tmp_probval
	cat(tmp_probval,"\n")
	tmp_split_likes[1+COO_weights_columnar[[1]][i]] = tmp_split_likes[1+COO_weights_columnar[[1]][i]] + tmp_probval
	}

split_likes = tmp_split_likes / Rsp_rowsums
split_likes


cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], COO_weights_columnar[[4]])
Rsp_rowsums

