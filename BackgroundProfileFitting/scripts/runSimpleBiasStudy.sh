#!/bin/tcsh

# ----- launch toys -------
scripts/sub_standardBias_jobs.py -D dat/standardBiasLegacy_massFac.dat -t 1000 -n 10 --skipPlots -q cmscaf1nh
# ----- make list of trees per job type  -------
scripts/prepareStdBiasToysTreeLists.csh /eos/cms/store/group/phys_higgs/meridian/BiasStudy/legacy_freeze_v2_massFac_c/ legacy_freeze_v2_massFac_c_trees
# ----- prepare histograms of pulls from trees -------
scripts/make_hists_simplebias_from_tree.py -d dat/bsConfig_2012_legacy_freeze_v2_massFac_hists.dat
# ----- to make final plots ------  
#scripts/draw_simplebias_plots.py
