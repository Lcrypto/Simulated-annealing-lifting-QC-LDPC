# Simulated Annealing lifting for QC-LDPC Codes and MET QC-LDPC Codes
The source code for the article "Simulated Annealing Method for Construction of High-Girth QC-LDPC Codes" by Vasiliy Usatyuk and Ilya Vorobyev, presented at the 41st International Conference on Telecommunications and Signal Processing (TSP) in Athens, Greece in July 2018, is now available. This code can be used to construct regular and irregular QC-LDPC codes with multiple edge type circulants from a protograph with the required minimal EMD value.

A short review about the related problem of Code on the Graph construction is available in English at http://www.mathnet.ru/php/presentation.phtml?option_lang=rus&presentid=17899.

This method outperforms all currently published algorithms such as PEG, QC-PEG, Fossorier-Declercq-Vasic Improved PEG, and Yedidia Hill-Climbing for constructing QC-LDPC codes from capability of cycles broken. For more details, please refer to the paper available at https://ieeexplore.ieee.org/document/8441303/ or https://www.researchgate.net/publication/327194285_Simulated_Annealing_Method_for_Construction_of_High-Girth_QC-LDPC_Codes.

In addition, when combined with code distance-based sieving (for short and moderate code lengths), this method allows for the construction of QC codes with very high coding gain, which is arguably the best for the current state-of-the-art.


![alt text](https://github.com/Lcrypto/Simulated-annealing-lifting-QC-LDPC/blob/master/Figure.png)


Table I. Minimal value of circulant for regular mother matrix with row number m=3 and column number n with girth 10

|Column number|Simulated annealing|Hill Climbing|Improved PEG|Lower Bound|
|-------------|-------------------|-------------|------------|-----------|
|4|37|39|37|37|
|5|61|63|61|61|
|6|91|103|91|91|
|7|155|160|155|127|
|8|215|233|227|168|
|9|304|329|323|217|
|10|412|439|429|271|
|11|545|577|571|331|
|12|709|758|-|-|

Table II. Minimal value of circulant for regular mother matrix with row number m=3 and column number n with girth 12

|Column number|Simulated annealing|Improved PEG|Table V|
|-------------|-------------------|-------------|------------|
|4|73|73|97|
|5|160|163|239|
|6|320|369|479|
|7|614|679|881|
|8|1060|1291|1493|
|9|1745|1963|2087|
|10|2734|-|-|
|11|4083|-|-|
|12|5964|-|-|

Constructed regular codes represented at tables contained in file "high-girth regular LDPC results.zip". 

Constructed 8 cyclic group decomposition MET QC-LDPC Codes families based on 5G eMBB Base Graph 2 contained in folder "SA results".

Simulated annealing lifting for high girth QC-LDPC include EMD optimization of protograph, which allow to decrease error-floor (by eliminate harmful Trapping Sets), or improve waterfall properties (by lifting of protograph with better threshold).
To use application call binary file with command:

*binary -file Your_protograph_file -circulant size_of_circulant -upGirth check_condition_up_to_cycles -emd EMD_values -seed initial_value_for_random_generator -numberOfMatrices number_of_requirement_matrix -girth girth_size*

Your_protograph_file:

*Contain number of columns(Variable nodes), rows (Check nodes)*

*1, 0 value for circulant permutation block matrix, obtain by some optimization of LDPC codes ensemble (Density Evolution, Covariance Evolution, PEXIT chart and etc).*


For Example:
simulatedAnnealingEMD.exe -file proto.txt -circulant 500 -upGirth 8 -emd 20 -seed 123 -numberOfMatrices 1 -girth 8
Proto.txt contain base matrix:

*3 2*

*1 0 1*

*1 1 0*

For construction multiple edge use 2,3,..., edges instead 1. 
Example:

*simulatedAnnealingEMD.exe -file proto.txt -circulant 500 -upGirth 6 -emd 2 -seed 123 -numberOfMatrices 1 -girth 8*

*16 6*

*1 0 0 0 0 1 0 1 0 1 1 0 1 0 0 2*

*1 1 0 0 0 0 0 1 1 0 0 1 0 1 1 2*

*0 1 1 0 0 0 0 1 0 1 0 1 0 1 0 1*

*0 0 1 1 0 0 1 0 1 0 1 0 0 1 1 1*

*0 0 0 1 1 0 1 0 0 1 1 0 1 0 1 2*

*0 0 0 0 1 1 1 0 1 0 0 1 1 0 0 2*

output:

*16	6	500*

*440	-1	-1	-1	-1	0	-1	258	-1	203	0	-1	237	-1	-1	329&215*

*53	75	-1	-1	-1	-1	-1	95	0	-1	-1	443	-1	0	238	284&257*

*-1	104	363	-1	-1	-1	-1	67	-1	466	-1	5	-1	174	-1	425*

*-1	-1	265	249	-1	-1	363	-1	59	-1	392	-1	-1	324	119	488*

*-1	-1	-1	260	64	-1	429	-1	-1	383	402	-1	421	-1	348	97&222*

*-1	-1	-1	-1	429	480	86	-1	234	-1	-1	114	41	-1	-1	392&402*


How to compile:

For Linux: compile source code from 'linux' folder by 'run.sh' compile script.  


For Windows: compile by MS VS 2015 project at 'simulatedAnnealingEMD' folder. Windows binary file, and shell scripts to run examples of lifting at folder 'simulatedAnnealingEMD/x64/Release'. If you not install MS VS 2015 to run binary file simulatedAnnealingEMD.exe, don't forget to download and install redistributed kit from https://www.microsoft.com/en-us/download/details.aspx?id=48145.


For EMD Spectrum calculation of QC-LDPC Code built by the simulated annealing method (or any othen method to compare cycle properties)  use  https://github.com/Lcrypto/EMD-Spectrum-LDPC .

For the QC MET-LDPC codes from above example, EMD Spectrum when evaluating the spectrum up to cycles of length 10  -upperGirth 10

*16	6	500*

*440	-1	-1	-1	-1	0	-1	258	-1	203	0	-1	237	-1	-1	329&215*

*53	75	-1	-1	-1	-1	-1	95	0	-1	-1	443	-1	0	238	284&257*

*-1	104	363	-1	-1	-1	-1	67	-1	466	-1	5	-1	174	-1	425*

*-1	-1	265	249	-1	-1	363	-1	59	-1	392	-1	-1	324	119	488*

*-1	-1	-1	260	64	-1	429	-1	-1	383	402	-1	421	-1	348	97&222*

*-1	-1	-1	-1	429	480	86	-1	234	-1	-1	114	41	-1	-1	392&402*

girth = 8	31000 cycles, cycle = 10	419000 cycles

|EMD value|8 cycles number| 10 cycles number|
|---|-------|------|
|3|	1000|	500|
|4|	500|	4500|
|5|	500|	2000|
|6|	3000|	4000|
|7|	2000|	15000|
|8|	1000|	19000|
|9|	1000|	14000|
|10|	2500|	26000|
|11|	2000|	20500|
|12|	4000|	25000|
|13|	1500|	38000|
|14|	1000|	26500|
|15|	2000|	24500|
|16|	2000|	37500|
|17|	1000|	19500|
|18|	1000|	43500|
|19|	500|	20000|
|20|	2000|	17500|
|21|	500|	7000|
|22|	0| 24500|
|24|	1500|10500|
|25|0|500|
|26|0|1500|
|27|0|	5500|
|28|	500|1500|
|29| 0|1500|




P.S. Compare Simulated Annealing lifting for 8 girth with method from paper "A. Kharin, A. Dryakhlov, E. Mirokhin, K. Zavertkin, A. Ovinnikov and E. Likhobabin, "An Approach to the Generation of Regular QC-LDPC Codes with Girth 8," 2020 9th Mediterranean Conference on Embedded Computing (MECO), Budva, Montenegro, 2020, pp. 1-4".




Table III. Minimal value of circulant for regular mother matrix with row number m=3 and column number n with girth 8 with running time constrain less than 24 hours

|Column number|Simulated annealing|Ovinnikov et al|
|-------------|-------------------|-------------|
|4|9|9|
|5|13|13|
|6|18|18|
|7|21|21|
|8|25|25|
|9|30|30|
|10|35|35|
|11|40|40|
|12|42|45|



The SA method was used to construct a 3x12 regular code with girth 8 in less than an hour on a multitread (AMD Ryzen 3950X) and less than 21 hours on a single thread (Intel i7700K). The files "12_3_42girth8upGirth6emd0seed11protograph_from_proto.txt_matrix629558.txt" and "12_3_42girth8upGirth6emd0seed333protograph_from_proto.txt_matrix828012.txt" are attached to GitHub.


Overall, the SA lifting method remains one of the best QC-LDPC protograph lifting methods published for breaking "flat" models (pure girth and girth+EMD maximization) of Trapping sets.


If you need to construct regular codes with maximal girth, consider using the kroneker-based approach results from Alireza Tasdighi and Emmanuel Boutillon's paper titled "Integer Ring Sieve for Constructing Compact QC-LDPC Codes with Girths 8, 10, and 12.," which has been submitted to IEEE Transactions on Information Theory in February 2021 ( http://www-labsticc.univ-ubs.fr/~boutillon/ldpc/ldpc.htm ). However, it is important to remember that maximum girth does not necessarily mean good Trapping set and low weight codeword spectrum (good EMD spectrum) performance/complexity, especially when compared to irregular LDPC codes which allow for trade-offs between waterfall and error-floor.



It is strongly recommended to improve not only the EMD spectrum but also the code (Hamming) distance, see https://github.com/Lcrypto/Length-und-Rate-adaptive-code . 



To achieve this goal, we implement the Lattice-based method in practice using the Kannan embedding, SVP(Shortest Vector Problem), and SBP (Block Korkin-Zolotarev, BKZ for solution Shortest Basis Problem) techniques. However, it is also possible to use the Dumer or Brouwer-Zimmerman algorithms implementation from GAP/MAGMA.
According to our (Usatyuk Vasiliy) results in the code distance challenge at https://decodingchallenge.org/low-weight/, Lattice methods are superior:
![alt text](https://github.com/Lcrypto/Length-und-Rate-adaptive-code/blob/master/Code_distance_challenge.png)

Lattice method for code distance estimation:


1.V. S. Usatyuk and S. I. Egorov, "Heuristic Number Geometry Method for Finding Low-Weight Codeword in Linear Block Codes," 2024 26th International Conference on Digital Signal Processing and its Applications (DSPA), Moscow, Russian Federation, 2024, pp. 1-6  https://ieeexplore.ieee.org/document/10510086



2.  Usatyuk V.S., E. Sergey and G. Svistunov, "Construction of Length and Rate Adaptive MET QC-LDPC Codes by Cyclic Group Decomposition," 2019 IEEE East-West Design & Test Symposium (EWDTS), Batumi, Georgia, 2019, pp. 1-5 https://ieeexplore.ieee.org/document/8884427

   
3. Usatyuk V.S., Yegorov S.I. Construction of quasi-cyclic non-binary LDPC codes, based on joint optimization of distance properties and connection spectra of codes. Telecommunications and Radio Engineering №8, 2016, pp. 32-40

   
4. Usatjuk V.S. Computing the minimum distance of nonbinary ldpc codes using block korkin-zolotarev method. Proceedings of the Southwest State University Series: Control, Computer Engineering, Information Science. Medical Instruments Engineering. - 2015. - № 3 (16). - pp. 76-85.



With BR,
Vasiliy.

