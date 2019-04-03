## BOLTZMANN LEARNING CODE
Code dca.cpp (c++ language) implements a parallelized sthocastic gradient descent procedure adopted to solve the direct coupling analysis inverse Potts problem on RNA families. It should be compiled with a C++11 supporting compiler. Recent GCC versions support C++11 by default, whereas older version
require `-std=c++11` or `-std=c++0x` flag.

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

Notice that the `run.sh` script preprocess the MSA first and produces a file named `seqs` with only the columns corresponding to the
non-gap positions in the first sequence of the alignment. An example `seqs.test` file already preprocessed is also included in this repository.
It can be analyzed using the command
````
mv seqs.test seqs
mpirun -np 20 dca.o 
````

### OUTPUT
- `conv_fi`, `conv_fij`: empirical frequency counts (one-site and two-sites respectively) vs model frequency counts (obtained via Monte Carlo sampling from the obtained Potts distribution)
- `h_vals`, `j_vals`: all parameters resulting from learning procedure (local fields and direct couplings respectively)
- scores: pairs covariance scores and respective nucleotides indexes (couplings Frobenius norm and subsequent APC correction)

