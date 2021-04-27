

library(optimx)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

# Setup
wd = "/drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/5species_ABBBB_ML_DECvDECj_LGpy/"
setwd(wd)


#######################################################
# Parse Python LAGRANGE output
#######################################################

states_fn = "ABBBB.results.txt"

# Remove this file, as it will append otherwise
splits_fn = "ABBBB.results_splits00001.txt"
if (file.exists(splits_fn))
	{
	file.remove(splits_fn)
	}

sumstats = parse_lagrange_python_output(outfn=states_fn, outputfiles=TRUE, results_dir=getwd(), new_splits_fn=TRUE, new_states_fn=FALSE, filecount=0, append=TRUE)
sumstats

splits_table_fn = sumstats$splits_table_fn
moref(splits_table_fn)


LGpy_splits = LGpy_splits_fn_to_table(splits_fn=splits_table_fn)
LGpy_splits

MLsplits = LGpy_MLsplit_per_node(splits=LGpy_splits)
MLsplits




# Plot LGpy node numbers
trfn = "tree_5species.newick"
tr = read.tree(trfn)
ntips = length(tr$tip.label)
Rnodenums = (ntips+1):(ntips+tr$Nnode)



# Nodenums for Python Lagrange
downpass_node_matrix = get_lagrange_nodenums(tr)
downpass_node_matrix

# Nodenums for Python Lagrange include
# tip numbers N0 - N(ntips-1)
# E.g. in this example they start at node N4
downpass_node_matrix[,2] = downpass_node_matrix[,2] + ntips - 1

# Sort in LAGRANGE node order
# Sort MLsplits the same way
if (nrow(downpass_node_matrix) > 1)
	{
	downpass_node_matrix = downpass_node_matrix[order(downpass_node_matrix[,2]),]
	# Sort MLsplits the same way
	MLsplits = MLsplits[order(MLsplits$nodenum_LGpy),]
	} else {
	downpass_node_matrix = downpass_node_matrix
	MLsplits = MLsplits
	}



# Get tip states
geogfn = "geog_5species.data"
tipranges = getranges_from_LagrangePHYLIP(geogfn)
tipranges

tip_areastrings = tipranges_to_area_strings(tipranges=tipranges, areaabbr=NULL)
tip_areastrings



pdffn = "LGpy_plots_v1.pdf"
pdf(pdffn, height=11, width=8.5)
par(mfrow=c(2,1))

plot(tr, show.tip.label=TRUE, label.offset=0.1)
axisPhylo()
nodelabels(node=Rnodenums, text=Rnodenums)
tiplabels(tip=1:ntips, text=1:ntips)
title("Tree with APE node numbers")



plot(tr, show.tip.label=TRUE)
axisPhylo()
nodelabels(node=downpass_node_matrix[,1], text=downpass_node_matrix[,2])
#tiplabels()
title("Tree with Python LAGRANGE node numbers")


plot(tr, show.tip.label=TRUE, label.offset=0.1)
axisPhylo()
nodelabels(node=downpass_node_matrix[,1], text=MLsplits$splits)
tiplabels(text=tip_areastrings, tip=1:ntips)
title("Python LAGRANGE most probable splits")


plot(tr, show.tip.label=TRUE, label.offset=0.1)
axisPhylo()
nodelabels(node=downpass_node_matrix[,1], pie=MLsplits$relprob)
tiplabels(text=tip_areastrings, tip=1:ntips)
title("Python LAGRANGE state probabilities\n(P(best state), then all others)")



dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)




