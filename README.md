# Simulated-annealing-lifting-QC-LDPC
Simulated annealing lifting for high girth QC-LDPC include EMD optimization of protograph, which allow to decrese error-floor (by eliminate harmfull Trapping Sets), or improve waterfall properties (by lifting of protograph with better treshold).
To use application call binary file with command:
binary -file Your_protograph_file -circulant size_of_circulant -upGirth check_condition_up_to_cycles -emd EMD_values -seed initial_value_for_random_generator -numberOfMatrices number_of_requirement_matrix -girth girth_size

Your_protograph_file
Contain number of columns(Variable nodes), rows (Check nodes)
1,0 value for circulant permutation block matrix, obtain by some optimization of LDPC codes ensemble (Density Evolution, Covariance Evolution, PEXIT chart and etc).
For Example:
simulatedAnnealingEMD.exe -file proto.txt -circulant 500 -upGirth 8 -emd 20 -seed 123 -numberOfMatrices 1 -girth 8
Proto.txt contain base matrix:
3 2
1 0 1
1 1 0
