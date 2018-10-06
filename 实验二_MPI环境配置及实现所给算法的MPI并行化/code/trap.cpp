#include <stdio.h>
#include <stdlib.h>
#include <mpich/mpi.h>
#include<iostream>
#include<cstdlib>
#include<ctime>
#include<vector>

using namespace std;

double f(double x) {
    double return_val;
 
    return_val = x*x;
    return return_val;
 }

double Trap(double left_endpt,double right_endpt,int trap_count,double base_len){

    double estimate,x;
    int i;

    estimate = ( f(left_endpt) + f(right_endpt) )/2.0;
    for(i = 0;i<=trap_count-1;i++){
        x = left_endpt + i*base_len;
        estimate += f(x);
    }
    estimate = estimate*base_len;

    return estimate;
}

int main(){

    int my_rank, comm_sz, n =100000000, local_n;
    double a=10.0,b=1000.0,h,local_a,local_b;
    double local_int,total_int;
    int source;

    double start_time;
    double end_time;
    double total_time;

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);

    start_time = MPI_Wtime();

    h = (b-a)/n;
    local_n = n/comm_sz;

    local_a = a+ my_rank*local_n*h;
    local_b = local_a + local_n*h;
    local_int = Trap(local_a,local_b,local_n,h);

    if(my_rank!=0){

        MPI_Send(&local_int,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);

    }else{
        total_int = local_int;
        for(source = 1; source < comm_sz; source++){
            MPI_Recv(&local_int,1,MPI_DOUBLE,source,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            total_int += local_int;
        }

    }

    end_time = MPI_Wtime();

    if(my_rank == 0){
        cout<<"Total time is "<<end_time - start_time<<" s"<<endl;
    }

    MPI_Finalize();
    return 0;

}