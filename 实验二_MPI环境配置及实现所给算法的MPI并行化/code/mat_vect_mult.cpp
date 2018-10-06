/* File:     mat_vect_mult.c
 *
 * Purpose:  Implement serial matrix-vector multiplication using
 *           one-dimensional arrays to store the vectors and the
 *           matrix.
 *
 * Compile:  gcc -g -Wall -o mat_vect_mult mat_vect_mult.c
 * Run:      ./mat_vect_mult
 *
 * Input:    Dimensions of the matrix (m = number of rows, n
 *              = number of columns)
 *           n-dimensional vector x
 * Output:   Product vector y = Ax
 *
 * Errors:   if the number of user-input rows or column isn't
 *           positive, the program prints a message and quits.
 * Note:     Define DEBUG for verbose output
 *
 * IPP:      Section 3.4.9 (pp. 113 and ff.), Section 4.3 (pp. 159
 *           and ff.), and Section 5.9 (pp. 252 and ff.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpich/mpi.h>
#include<iostream>
#include<cstdlib>
#include<ctime>
#define random(a,b) (rand()%(b-a+1)+a)
#include<vector>

using namespace std;


/*-------------------------------------------------------------------
 * Function:   Read_matrix
 * Purpose:    Read the contents of the matrix from stdin
 * In args:    prompt:  description of matrix
 *             m:       number of rows
 *             n:       number of cols
 * Out arg:    A:       the matrix
 */
void Read_matrix(
	char    prompt[]  /* in  */,
	double  A[]       /* out */,
	int     m         /* in  */,
	int     n         /* in  */) {
	int i, j;

	printf("creating the matrix %s\n", prompt);
	for (i = 0; i < m; i++) {

		for (j = 0; j < n; j++) {

			//scanf("%lf", &A[i*n + j]);
			A[i*n + j] = random(1, 10);

		}

	}


}

/*-------------------------------------------------------------------
 * Function:   Read_matrix
 * Purpose:    Read a vector from stdin
 * In args:    prompt:  description of matrix
 *             n:       order of matrix
 * Out arg:    x:       the vector being read in
 */
void Read_vector(
	char    prompt[]  /* in  */,
	double  x[]       /* out */,
	int     n         /* in  */) {
	int i;

	printf("creating the vector %s\n", prompt);
	for (i = 0; i < n; i++)
		x[i] = random(1, 10);
	//scanf("%lf", &x[i]);
}  /* Read_vector */



/*-------------------------------------------------------------------*/
int main(int argc, char *argv[]) {

      int myid, numprocs;
      MPI_Init(NULL, NULL);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);



      double* A = NULL;
	double* x = NULL;
	//double* y = NULL;

	char ca[2] = "A";
	char cx[2] = "x";
	char cy[2] = "y";

      //int* m[1], n[1];
      int * m = new int[1];
      int * n = new int[1];
      
      //int* section_num[1];
      //int* last_sectionnum[1];
      int * section_num = new int[1];
      int * last_sectionnum = new int[1];

      double * result;
      double * local_result;
      int local_i = 0;

      double  start_time = 0;
      double  end_time = 0;
      double  diff_time = 0;
      double  total_time = 0;


      if(myid ==0){

      m[0] = 10000; //  atoi(argv[1]);
      n[0]=10000; //atoi(argv[2]);

      section_num[0] = m[0]/numprocs;
      last_sectionnum[0] = section_num[0] + m[0]%numprocs;

      result = new double[m[0]];

      }

      
      MPI_Bcast(m,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(section_num,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(last_sectionnum,1,MPI_INT,0,MPI_COMM_WORLD);
      
      

      A = new double[m[0]*n[0]];
      x = new double[n[0]];

	if (A == NULL || x == NULL) {
		fprintf(stderr, "Can't allocate storage\n");
		exit(-1);
      }
      if(myid != numprocs-1){
            local_result = new double[section_num[0]];
      }else{
            local_result = new double[last_sectionnum[0]];
      }

      

      if(myid ==0){
	Read_matrix( ca, A, m[0], n[0]);
      Read_vector(cx, x, n[0]);


      }



      MPI_Bcast(A,m[0]*n[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(x,n[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
      

      start_time = MPI_Wtime();

      if(myid !=numprocs-1 )
      {
            //cout<<"myid"<<myid<<endl;

            for(int i= myid*section_num[0]; i<myid*section_num[0]+section_num[0]; i++) {
                  local_result[local_i] = 0;
                  for(int j=0;j<n[0];j++){
                        local_result[local_i] = local_result[local_i] + A[i*n[0]+j]*x[j];
                        //cout<<"A[i*n[0]+j] "<<A[i*n[0]+j]<<" x[j] "<<x[j]<<" local_result[local_i] "<<local_result[local_i]<<endl;
                        }
                  local_i++;
                  }

      }else{

            //cout<<"myid"<<myid<<endl;

            for(int i= myid*section_num[0]; i<myid*section_num[0]+last_sectionnum[0]; i++) {
                  local_result[local_i] = 0;
                  for(int j=0;j<n[0];j++){
                        local_result[local_i] = local_result[local_i] + A[i*n[0]+j]*x[j];
                        //cout<<"A[i*n[0]+j] "<<A[i*n[0]+j]<<" x[j] "<<x[j]<<" local_result[local_i] "<<local_result[local_i]<<endl;                        
                        }
                  local_i++;
                  }

      }

      end_time = MPI_Wtime();
      diff_time = end_time - start_time;


      MPI_Gather(local_result, section_num[0], MPI_DOUBLE, result, section_num[0], MPI_DOUBLE, 0 , MPI_COMM_WORLD);


      if(myid == 0){
            cout<<"total time is "<<diff_time<<" s"<<endl;
      }


  

   free(A);
   free(x);
   free(local_result);
   

   MPI_Finalize();
   return 0;
}  /* main */


