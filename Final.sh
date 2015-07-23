#$1: get_iap_meth_stats.py: Python file giving out JSON
#$2: iap_whole_genome.merged.filtered.bed: whole genome iap/erv regions file
#$3: var_meth.bed: ragged regions validated by Anastasiya
#$4: BL6-F-108_bedgraph,BL6-F-138_bedgraph,BL6-M-88_bedgraph,BL6-M-98_bedgraph,C57B-M-7905954,C57B-M-7905953,C57B-F-7905948,C57B-F-7905947: replicates methylation data
#$5: 0: number of random regions or 0 for all
#$6: 20: min number of CpGs
#$7: 200: min number of bps
#$8: stdin_json.R: R file giving out variance data
#$9: result_23_02.csv: csv file where variance data is to be saved
#$10: hyperbola.py: Python file to optimise hyperbola parameters
#$11: checked.bed: training data file
#$12: best1.txt: file where hyperbola parameters are/will be saved
#$13: 0: if training data has been changed or not- 0/1
#$14: Plotting.R: R file plotting regions and hyperbola
#$15: regions.csv: csv file name to save ragged regions from hyperbola approach
#$16: svm.R: R file to perform svm
#$17: svm_regions.csv: csv file name to save ragged regions from svm

#takes in replicates methylation data and gives out calculated methylation variance in a csv file 
python $1 $2 $3 $4 $5 $6 $7 | Rscript $8 $9

#takes in training data set and gives out best hyperbola parameters, with number of incorrect regions 
python $10 $11 $12 $13

#takes in hyperbola parameters and gives out ragged regions file and also a graph
Rscript $14 $12 $9 $11 $15

#takes in training data and predicts ragged regions using svm
Rscript $16 $11 $9 $17

#bash Final.sh get_iap_meth_stats.py iap_whole_genome.merged.filtered.bed var_meth.bed BL6-F-108_bedgraph,BL6-F-138_bedgraph,BL6-M-88_bedgraph,BL6-M-98_bedgraph,C57B-M-7905954,C57B-M-7905953,C57B-F-7905948,C57B-F-7905947 0 20 200 stdin_json.R result_23_02.csv hyperbola.py checked.bed best1.txt 0 Plotting.R regions.csv svm.R svm_regions.csv