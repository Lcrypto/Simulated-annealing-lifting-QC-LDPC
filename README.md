# Simulated Annealing lifting for QC-LDPC Codes and MET QC-LDPC Codes
Source code for article
Vasiliy Usatyuk und Ilya Vorobyev 
Simulated Annealing Method for Construction of High-Girth QC-LDPC Codes
 41st International Conference on Telecommunications and Signal Processing (TSP) 2018, 4-6 Jule, Athens, Greece.
 
It construct regular and irregular QC-LDPC codes with  multiple edge type circulants from protograph with required minimal EMD value. Short review (ENG) about Code on the Graph construction related problem http://www.mathnet.ru/php/presentation.phtml?option_lang=rus&presentid=17899.
 
Simulated annealing superior all currently published algoriths (PEG,QC-PEG, Fossorier-Declercq-Vasic Improved PEG, Yedidia Hill-Climbing) and  for construction QC-LDPC codes from capability of cycles broken, for detail read paper https://ieeexplore.ieee.org/document/8441303/  https://www.researchgate.net/publication/327194285_Simulated_Annealing_Method_for_Construction_of_High-Girth_QC-LDPC_Codes. With combination of code distance based sieving (for short and moderate code length) it allow to contruct QC codes with very good (probably best for current state of art) coding gain.


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


For EMD Spectrum calculation of QC-LDPC Code built by the simulated annealing method  use  https://github.com/Lcrypto/EMD-Spectrum-LDPC .

For the QC MET-LDPC codes from above 16	6	500  EMD Spectrum when evaluating the spectrum up to cycles of length twelve -upperGirth 12  

*girth = 8	31000 cycles
EMD	Number
3	1000
4	500
5	500
6	3000
7	2000
8	1000
9	1000
10	2500
11	2000
12	4000
13	1500
14	1000
15	2000
16	2000
17	1000
18	1000
19	500
20	2000
21	500
24	1500
28	500
girth = 10	419000 cycles
EMD	Number
3	500
4	4500
5	2000
6	4000
7	15000
8	19000
9	14000
10	26000
11	20500
12	25000
13	38000
14	26500
15	24500
16	37500
17	19500
18	43500
19	20000
20	17500
21	7000
22	24500
23	9000
24	10500
25	500
26	1500
27	5500
28	1500
29	1500
girth = 12	7428000 cycles
EMD	Number
2	1000
3	1500
4	20000
5	27000
6	21000
7	55500
8	141500
9	152000
10	212500
11	277500
12	338000
13	448000
14	571000
15	429000
16	526000
17	591000
18	460000
19	577500
20	401500
21	327500
22	422000
23	318500
24	252500
25	215000
26	158500
27	108000
28	130500
29	61500
30	49500
31	38000
32	20500
33	29500
34	13000
35	19000
36	13000*




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



SA method constructed 3x12 regular code with girth 8 less than one hour on multitread (Amd Ryzen 3950x) and less than 21 hours on 1 thread (Intel i7700k), file "12_3_42girth8upGirth6emd0seed11protograph_from_proto.txt_matrix629558.txt", "12_3_42girth8upGirth6emd0seed333protograph_from_proto.txt_matrix828012.txt",  attached to github. 





Summary: SA lifting method still one of the best published QC-LDPC protograph lifting method from "flat" model(pure girth and girth+EMD maximization) of Trapping sets breaking.


If you need to construct regular codes with maximal girth use kroneker based approach results from
Alireza Tasdighi, Emmanuel Boutillon, "Integer Ring Sieve for Constructing Compact QC-LDPC Codes with Girths 8, 10, and 12.", Submitted to IEEE Trans. on Information Theory, Feb. 2021. http://www-labsticc.univ-ubs.fr/~boutillon/ldpc/ldpc.htm
but alwayse remember max girth not mean good Trapping set and low weight codeword Spectrum (good EMD Spectrum) - performance/complexity especialy compare to irregular LDPC which allow trade-off between waterfall and error-floor.


With BR,
Vasiliy.

