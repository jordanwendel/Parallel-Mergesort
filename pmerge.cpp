/**
    Pmerge Project
    Authors: Jordan Wendel and Katelyn Murphy
    Sorts an array of n size in parallel
*/

// New compile and run commands for MPI!
// mpicxx -o blah file.cpp
// mpirun -q -np 4 blah

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include <cmath>

using namespace std;

// ***** GLOBAL VARIABLES ***** //
int my_rank; // my CPU number for this process
int p; // number of CPUs

void QuickSort (int * a, int first, int last)
{ // Quicksort algorithm, run time nlogn
    int i = first, j = last;
    int tmp;
    int pivot = a[(first + last) / 2];
    
    while (i <= j) {
        while (a[i] < pivot)
            i++;
        while (a[j] > pivot)
            j--;
        if (i <= j) {
            tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
            i++;
            j--;
        }
    }

    if (first < j)
        QuickSort(a, first, j);

    if (i < last)
        QuickSort(a, i, last);
}

void PrintArray (int * a, int size) 
{ // Print array values to reduce redundant code

    for (int x=0; x<size; x++) {
        cout << a[x] << " ";
    }
    cout << endl;
}

// Helper function to remove redundant code in smerge
void move (int * a1, int first1, int last1, int * a2, int first2)
{
    while (first1 <= last1)
        a2[first2++] = a1[first1++];
}

int rank (int * a, int first, int last, int valToFind)  
{ // Binary Search to find the rank
    int middle = 0;
    int count = last;

    while (first <= last) {
        // First compare valToFind with middle element
        middle = (last + first) / 2;

        // If target value is greater than middle, continue search in upper half of array
        if (a[middle] < valToFind) {
            first = middle + 1;
        }
        // If target value is less than middle, continue search in lower half of array
        else if (a[middle] > valToFind) {
            last = middle - 1;
        }
        // If target value matches middle, return its position
        else {
            return middle;
        }
    }
}

void smerge(int *a, int first1, int last1, int * b, int first2, int last2, int * c, int first3)
{
    int index1, index2, index3;

    index1 = first1;
    index2 = first2;
    index3 = first3;

    //Merge two lists together until one list runs out of items

    while ((index1 <= last1) && (index2 <= last2)) {
        if (a[index1] < b[index2])
            c[index3++] = a[index1++];
        else
            c[index3++] = b[index2++];
    }
    move (a, index1, last1, c, index3);

    move (b, index2, last2, c, index3);
}


void pmerge (int* a, int* aLeft, int* aRight, int n, int my_rank, int p)
{

// ************************ SRANKA and SRANKB ******************** //
    int local_start = my_rank;
    int q = n/2;
    int loggy = log(q) / log(2);

    // Size of stripe
    int num = ceil((double(q) / (loggy))); // Ceiling command will round down to nearest whole number
   
    // Create arrays
    int * local_sranka = new int [num];
    int * local_srankb = new int [num];
    int * SRANKA = new int [num];
    int * SRANKB = new int [num];

    // Fill arrays with 0 so that the local values can be reduced into the correct position
    for (int i = 0; i < num; i++) {
		local_sranka[i] = 0;
        local_srankb[i] = 0;
		SRANKA[i] = 0;
        SRANKB[i] = 0;
    }

    // Print both sides of array
    if (my_rank==0) {
        cout << "Left half of array: " << endl;
        PrintArray(aLeft, q);
		cout<<endl;

        cout << "Right half of array: " << endl;
        PrintArray(aRight, q);
    }

    // Calculate local sranks by striping across processors
    for (int x = local_start; x < num; x += p)
    {
        local_sranka[x] = rank(aRight, 0, q-1, aLeft[x*loggy]);
        local_srankb[x] = rank(aLeft, 0, q-1, aRight[x*loggy]);
    //     cout<< "local_sranka: "<< local_sranka[x] << " ";
    //     cout << "local_srankb: " << local_srankb[x] << " ";
    }

    // Reduce local sranks
    MPI_Allreduce(local_sranka, SRANKA, num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_srankb, SRANKB, num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   
    if (my_rank == 0){
        // Print SRANKA
        cout<< endl << "SRANKA all reduced: "<< endl;
        cout << endl;
        PrintArray(SRANKA, num);
        cout << endl;

        // Print SRANKB
        cout<< endl << "SRANKB all reduced: "<< endl;
        cout << endl;
        PrintArray(SRANKB, num);
        cout << endl;
    }

    // Broadcast SRANKA and SRANKB to every processor from root 0.
    MPI_Bcast(SRANKA, p, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(SRANKB, p, MPI_INT, 0, MPI_COMM_WORLD);

// ****************************** SHAPE GAME ****************************** //

    // Partition problem into size of 2 * num
    // Each srank element represents a block
    // These 2 arrays will hold unsampled ranks that make up the shapes
    int * aShape = new int[2*num];
    int * bShape = new int[2*num];
    
    /**
    Plan is to fill the first part of array with end points (dotted blocks),
     then fill the rest of the array with the sampled ranks (shapes),
     then sort the in the array.
    
    Determine endpoints of each shape of size loggy that falls between each SRANKA or SRANKB element.
     Each endpoint will act as the start of the next shape.
     Fill the first indexes from 0 to num with the endpoints.
    **/
    for (int i = 0; i < num; i++) {
        aShape[i] = loggy * i;
        bShape[i] = loggy * i;
    }

    /**
    aShape and bShape are subarrays being filled with the sampled ranks
     of B for aShape and A for bShape after the first num indexes between
     start point of A and endpoint of B (for aShape)
     and start point of B and endpoint of A (for bShape).
    **/
    for (int x = 0; x < num; x++) {
        aShape[num + x] = SRANKB[x];
        bShape[num + x] = SRANKA[x];
    }

    
    //Sorting the endpoints and the sampled ranks together to make organized shapes for the shape game
    QuickSort(aShape, 0, 2 * num);
    QuickSort(bShape, 0, 2 * num);

    if (my_rank == 0) {
        // Print both shape arrays
        cout << "This is aShape: " << " ";
        PrintArray(aShape, 2 * num);
        cout << "This is bShape: " << " ";
        PrintArray(bShape, 2 * num);
        cout << endl;
    }

    // Local array playing shape game
    int * corona_local = new int[n];
    int * WIN = new int[n];

    // Fill arrays with 0's to be reduced later
    for (int i = 0; i < n; i++) {
        corona_local[i] = 0;
        WIN[i] = 0;
    }

    
    // Smerge all the shapes but the last using striping
    for (int x = local_start; x < 2 * num - 1; x += p) {
        
        smerge (aLeft, aShape[x], aShape[x+1]-1, aRight, bShape[x], bShape[x+1]-1, corona_local, aShape[x] + bShape[x]);
        
        if (my_rank == 0) {
            cout << "Sorted smerged shape when x= " << x
                << endl;
            PrintArray(corona_local, n);
            cout << endl << endl;
        }
    }


    if (my_rank == 0) {
        //smerge the last shape
        smerge(aLeft, aShape[2*num-1], n/2-1, aRight, bShape[2*num-1], n/2-1, corona_local, aShape[2*num-1] + bShape[2*num-1]);

        cout<< "Last smerge shape: " <<endl;
        PrintArray(corona_local, n);
        cout<<endl;
    }

    MPI_Allreduce(corona_local, WIN, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   
    for (int x = 0; x < n; x++)
        a[x] = WIN[x];
}

void mergesort(int *a, int n, int first, int last)
{
    int middle;
    
    // Base case
    if (n <= 32) {
        QuickSort(a, 0, n-1);
    }

    else if (first < last)
    {
        middle = n/2; //pivot

        //recursively merge sort
        mergesort(&a[0], middle, my_rank, p);
        mergesort(&a[middle], middle, my_rank, p);
        pmerge(a, &a[0], &a[n/2], n, my_rank, p);
    }
}

int main (int argc, char * argv [])
{ // number of CPUs that we have
    int source; // rank of the sender
    int dest; // rank of destination
    int tag = 0; // message number
    MPI_Status status; // return status for receive

    // Start MPI
    MPI_Init(&argc, &argv);

    // Find out my rank!
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Find out the number of processes!
    MPI_Comm_size(MPI_COMM_WORLD, &p);

 
   

    // ***** CREATE ARRAY  ***** //

    
    //CHANGE SIZE TO WHATEVER VALUE FOR N
	int n = 128;

	int* a = new int[n];  // creates an array that will hold all the random letters

    if (my_rank == 0) {
        srand(time(NULL)); //seed

        // Prints unsorted array items
        cout << "Unsorted: ";

        // Assigns random numbers to array a without duplicates
        for (int i = 0; i < n; i++)
        {
            bool check; // Checks if number is a duplicate
            int x; // Temporary "check" variable that holds the value of the random number in case it needs to be changed

            // do loop that will check if the number is a repeat
            do {
                x = rand() % 500;
                check = true; // Not a repeat (base state)
                for (int j = 0; j < i; j++) { // If the random number is already in array
                    if (x == a[j]) {
                        check = false; // It is a repeat
                        break;
                    }
                }
            } while (!check); // Runs do loop while check is false

            a[i] = x;
            cout << a[i] << " ";
        }
        cout << endl;
        cout << endl;
    }

	MPI_Bcast(&a[0], n, MPI_INT, 0, MPI_COMM_WORLD); // Broadcasts a of size n from root 0
	
	mergesort(a, n, 0, n);

	if (my_rank == 0) {
        cout << "Not even a global pandemic could stop us, praise the good Lord, here is the sorted array: ";
        PrintArray(a, n);
    }

	delete [] a;
	MPI_Finalize();
	return 0;
}