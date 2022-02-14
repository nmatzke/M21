#!/usr/bin/env python
import os
import lagrange

notes = '''
cd /drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/9species_A_8AB_ML_DECvDECj_LGpy/

python A_8AB.lagrange.py

Rscript plot_LGpy_result_v1.R

'''


data = """\
### begin data
{'area_adjacency': [[1, 1], [1, 1]],
 'area_dispersal': [[[1.0, 1.0], [1.0, 1.0]]],
 'area_labels': 'AB',
 'base_rates': '__estimate__',
 'dispersal_durations': [10.0],
 'dm_symmetric_entry': True,
 'excluded_ranges': [],
 'lagrange_version': '20130526',
 'max_range_size': 2,
 'model_name': 'A_8AB',
 'newick_trees': [{'included': '__all__',
                   'name': 'Tree0',
                   'newick': '(sp9:8,(sp8:7,(sp7:6,(sp6:5,(sp5:4,(sp4:3, (sp3:2, (sp2:1, sp1:1)N9:1)N10:1)N11:1)N12:1)N13:1)N14:1)N15:1)N16;',
                   'root_age': 8.0}],
 'ranges': [(),
            (0,),
            (1,),
            (0, 1)],
 'taxa': ['sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8', 'sp9'],
 'taxon_range_data': {'sp1': (0,), 'sp2': (0,1), 'sp3': (0,1), 'sp4': (0,1), 'sp5': (0,1), 'sp6': (0,1), 'sp7': (0,1), 'sp8': (0,1), 'sp9': (0,1)}}
### end data
"""

i = 0
while 1:
    if not i:
        outfname = "A_8AB.results.txt"
    else:
        outfname = "A_8AB.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = open(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)
