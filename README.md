# Summary
This script takes as input the combined Orthogroups.tsv and Orthogroups_Unassigned.tsv files from orthofinder and calcuates the core/shell/cloud of the investiaged genomes. 

# Requirements
pandas, python3 

# Definitions
core genes: genes present in all genomes in any ratio (N:N)
shell genes: genes present all but one genome
cloud genes: genes present in only one genome, single copy or multiple copies

# Example command
run orthofinder
cat Orthogroups.tsv Orthogroups_Unassigned.tsv > Orthogroups.combined.tsv

python ./calc_csc_OG.py

