#!/usr/bin/env python
import os
import lagrange

notes = '''
cd /drives/GDrive/__GDrive_projects/2018-01-02_Ree_Sanmartin/_example_scripts/6species_ABBBBB_ML_DECvDECj_LGpy/

python ABBBBB.lagrange.py

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
 'model_name': 'ABBBBB',
 'newick_trees': [{'included': '__all__',
                   'name': 'Tree0',
                   'newick': '(sp6:5,(sp5:4,(sp4:3, (sp3:2, (sp2:1, sp1:1)N6:1)N7:1)N8:1)N9:1)N10;',
                   'root_age': 5.0}],
 'ranges': [(),
            (0,),
            (1,),
            (0, 1)],
 'taxa': ['sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6'],
 'taxon_range_data': {'sp1': (0,), 'sp2': (1,), 'sp3': (1,), 'sp4': (1,), 'sp5': (1,), 'sp6': (1,)}}
### end data
"""

i = 0
while 1:
    if not i:
        outfname = "ABBBBB.results.txt"
    else:
        outfname = "ABBBBB.results-"+str(i)+".txt"
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
