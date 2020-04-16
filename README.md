# Parallel-Mergesort
Utilize MPI to do mergesort in parallel
The "pmerge" function is what the actual change is between a regular mergesort and the parallel mergesort.

First, the entire non-duplicated random array a will be split in half.
Then, both halves are quicksorted (essentially 2 sorted arrays glued together as one array).
Halves are called A and B for theory purposes.
Those halves are passed to the pmerge function.

Pmerge does the work in stages:
1) Calculate sampled ranks by striping the work across the processors.
  - We want to calculate the rank position within B of a value in A.
  - However, we don't want to calculate every rank for each value (takes too much time)
  - So, we take sampled rank values and pick them based on striping the work among processors
  - The local sranka and srankb values are then MPI reduced into the global sranka and srankb
2) Determine what work needs to be done with a sequential merge function
  - When we only sampled certain ranks, we left lots of small problems left to be solved. That's where
    smerge comes in.
   - Partition between the srank values and smerge them
3) Reduce values together with MPI
  - Reduce the local smerge values to the array WIN
