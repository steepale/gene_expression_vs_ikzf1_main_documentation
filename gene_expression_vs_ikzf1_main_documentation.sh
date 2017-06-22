#===============================================================================
#
#         FILE: /Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1/gene_expression_vs_ikzf1_main_documentation.sh
#
#        USAGE: for documentation purposes, scripts inside
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for determining if IKZF1 perterbation affects gene expression
#                
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2017.06.22
#     REVISION:  
#===============================================================================

# PROJECT DIRECTORY (TUM Cluster)
proj_dir="/Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1"
cd $proj_dir

# Make proper directories
mkdir -p ./{data,scripts,analysis,jobs}

# Sources to cite:
# http://f1000researchdata.s3.amazonaws.com/manuscripts/9996/16888710-6725-433b-ab2a-13c04bfe6cb5_8987_-_gordon_smyth_v2.pdf
# http://2014-msu-rnaseq.readthedocs.io/en/latest/model-orgs.html
# https://wikis.utexas.edu/display/bioiteam/Differential+gene+expression+analysis#Differentialgeneexpressionanalysis-Optional:edgeR

# Generate the targets file
./data/targets_IKZF1.txt



### Run the majority of the analysis in R (interactive in R studio)
Rscript --vanilla ./scripts/gene_expression_vs_ikzf1.R

# ./scripts/gene_expression_vs_ikzf1.R
#######################################
