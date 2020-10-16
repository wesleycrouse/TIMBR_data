Supplemental Material for Crouse, Kelada, and Valdar, 2020.pdf - the supplemental text, figures, and tables for the TIMBR manuscript

TIMBR-master.zip - frozen version of the 'TIMBR' R package, also available at https://github.com/wesleycrouse/TIMBR/
install_TIMBR.R - script for installing TIMBR, either locally or from GitHub

simulation/TIMBR_sim.R - R script for simulations with command line inputs for parallelization
simulation/TIMBR_sim.sh - Bash script for submitting TIMBR_sim.R in parallel 
simulation/TIMBR_sim_collapse_results.R - R script for collapsing and processing parallel sets of results
simulation/TIMBR_sim_plot.R - R script for generating plots
simulation/results/* - RData files containing the simulation results used in the manuscript

dspr/data/* - data files from King et al. (2014)
dspr/dspr_locus1.R - R script to generate diplotype probabilties for "multiallelic" example
dspr/dspr_locus1.R - R script to generate diplotype probabilties for "biallelic" example
dspr/TIMBR_dspr_locus1.R - R script to run TIMBR and generate plots for "multiallelic" example
dspr/TIMBR_dspr_locus2.R - R script to run TIMBR and generate plots for "biallelic" example

precc_cache/full.zip - compressed diplotype probabilities for PreCC data in HAPPY cache format; must be unzipped

eqtl/data/* - data files from Kelada et al. (2014)
eqtl/eqtl_list.R - R script to subset eQTL list for analysis
eqtl/TIMBR_eqtl.R - R script for TIMBR eQTL analysis with command line inputs for parallelization
eqtl/TIMBR_eqtl.sh - Bash script for submitting TIMBR_eqtl.R in parallel
eqtl/TIMBR_eqtl_collapse_results.R - R script for collapsing and processing parallel sets of results
eqtl/TIMBR_eqtl_plot_results.R - R script for generating plots

mcv_tree/tree_pipeline.sh - Bash script for generating trees for constructing pseudogenome and generating trees at MCV locus
mcv_tree/founders.txt - a list of founder strains for the CC
mcv_tree/four_gamete_test.R - R script for performing the four gamete test to identify genomic regions without evidence of recombination
mcv_tree/construct_pseudogenome.R - R script for constructing pseudogenomes for a genomic region
mcv_tree/chr7_3776.fa - FASTA for pseudogenome at MCV locus
mcv_tree/beast1_custom_template.xml - template for generating BEAST input files using BEASTGen
mcv_tree/chr7_3776.xml - BEAST input file at MCV locus
mcv_tree/chr7_3776.trees - BEAST tree samples at MCV locus
mcv_tree/chr7_3776.log - BEAST log file for tree samples at MCV locus
mcv_tree/convert_to_coalescent.R - R script to converting BEAST tree samples to coalescent units
mcv_tree/mcv_trees.RData - RData file of trees at MCV locus
mcv_tree/mcv_trees.trees - NEXUS file of trees at MCV locus; visualized with DensiTree
mcv_tree/random_coalescent.R - R file to sample random coalescent trees and output as NEXUS file

mcv/data/mcv_pheno.txt - data file from Kelada et al. (2012)
mcv/TIMBR_mcv.R - R script for to run TIMBR and generate plots for the MCV QTL
