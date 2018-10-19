## BOLTZMANN LEARNING CODE
Code dca.cpp (c++ language) implements a parallelized sthocastic gradient descent procedure adopted to solve the direct coupling analysis inverse Potts problem on RNA families.

### EXECUTION
Script run.sh must be executed. It requires two arguments: 
- first argument  must be the input file (any name), a multiple sequence alignment of homolougus RNA sequences in FASTA format (alphabeth `AUCGaucg-`)
- second argument must be the number of processes to be run in parallel
The script takes care of proper manipulation of the input file and of compiling/executing the main code dca.cpp in parallel (MPI protocol) 

### EXAMPLE
````
./run.sh MSA 20
````
where MSA is a reported input file example and 20 would be the required number of processors

### OUTPUT
- `conv_fi`, `conv_fij`: empirical frequency counts (one-site and two-sites respectively) vs model frequency counts (obtained via Monte Carlo sampling from the obtained Potts distribution)
- `h_vals`, `j_vals`: all parameters resulting from learning procedure (local fields and direct couplings respectively)
- scores: pairs covariance scores and respective nucleotides indexes (couplings Frobenius norm and subsequent APC correction)

