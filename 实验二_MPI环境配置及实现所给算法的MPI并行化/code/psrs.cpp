#include <stdio.h>
#include <stdlib.h>
#include <mpich/mpi.h>
#include<iostream>
#include<cstdlib>
#include<ctime>
#include<vector>
#include<fstream>
#include<algorithm>
#include<math.h>

using namespace std;

int dn = pow(2,20)*128;

#define data_num 40000000


void read(int * &buffer, long size) {
    
        ifstream infile;
        infile.open("../../data/data128M.edata", ios::binary);
        if (!infile.is_open()) {
            cout << "error\n";
            exit(0);
        }
        buffer = new int[size];
        infile.read((char*)buffer, sizeof(int)*size); 
        infile.close();
    
    
    }

int main(int argc, char *argv[]) {
        
    int myid, numprocs;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int * data; //store the global data
    int data_per_node = data_num/numprocs ;//the number of data for each node to process
    int * local_data = new int[data_per_node]; //store the local data, after partition

    int * local_sample = new int[numprocs]; //store the local sample
    int count = 0;

    int * root_sample = new int[numprocs*numprocs];

    int * pivot = new int[numprocs-1];

    int * local_partion_num = new int[numprocs];
    for(int i=0;i<numprocs;i++){
        local_partion_num[i] = 0;
    }
    int * local_partion_num_all = new int[numprocs*numprocs];
    // for(int i=0;i<numprocs;i++){
    //     local_partion_num_matrix[i]= new int[numprocs];
    //     // for(int j=0;j<numprocs;j++){
    //     //     local_partion_num_matrix[i][j] = 0;
    //     // }
    // }

    //store the number of data each compute node should process
    int* local_partion_num_sum = new int[1];
    

    //store the location of each pivot in each compute node
    int* local_partition_location = new int[numprocs];
    for(int i=0;i<numprocs;i++){
        local_partition_location[i] = 0;
    }
    int * local_partition_location_all = new int[numprocs*numprocs];


    



    int * after_partition_data;



    //read data
    if(myid ==0){
        
        data = new int [data_num];
        //read data
        read(data,data_num);

        // for(int i=0;i<data_num;i++){
        //     cout<<data[i]<<" ";
        //     if(i%10 == 0)
        //         cout<<endl;
        // }
    }


    double start_time = MPI_Wtime();
    


    //scatter the data to each compute node 
    MPI_Scatter(data,data_per_node,MPI_INT,local_data,data_per_node,MPI_INT,0,MPI_COMM_WORLD);

    //sort the local data for each node
    sort(local_data,local_data+data_per_node);

    //sample
    for(int i=0;i<numprocs;i++){
        //cout<<"local sample:"<< local_data[  i * (data_per_node / (numprocs + 1))    ] <<endl;
        local_sample[i] = local_data[  i * (data_per_node / (numprocs + 1))  ];
    }

    //gather the sample
    MPI_Gather(local_sample,numprocs,MPI_INT,root_sample,numprocs,MPI_INT,0,MPI_COMM_WORLD);

    //sort the sample in the root and pick the pivot
    if(myid == 0){

        //sort the local sample
        sort(root_sample,root_sample+numprocs*numprocs);

        //select the pivot
        for(int i=0;i<numprocs-1;i++){
            pivot[i] = root_sample[  (i+1)*(numprocs+1)  ];
        }

    }

    //broadcast the pivot to all compute node
    MPI_Bcast(pivot,numprocs-1,MPI_INT,0,MPI_COMM_WORLD);

    
    //traverse the local_data and record the number for each partition
    for(int i=0;i<data_per_node;i++){

        //record the number for each partition
        for(int j=0;j<numprocs-1;j++){

            //record the partition location
            // if(local_data[i]<=pivot[j] && local_data[i+1]>pivot[j]){
            //     //cout<<"local_data[i] "<<local_data[i]<<" local_data[i+1] "<< local_data[i+1] <<"  pivot "<<pivot[j]<<endl;
            //     // if(j == 0){
            //     //     local_partition_location[j] = i;
            //     // }
            //     // else{
            //     //     local_partition_location[j] = i - local_partition_location[j-1];
            //     //     if(j == numprocs-2){
            //     //         local_partition_location[j+1] = data_per_node-i;
            //     //     }
            //     // }
            //     local_partition_location[j+1] = i;                        
            // }

            if(local_data[i]<pivot[j]){
                local_partion_num[j] = local_partion_num[j]+1;
                break;
            }

            if( local_data[i]>=pivot[j] && j == numprocs -2  ){
                local_partion_num[j+1] = local_partion_num[j+1]+1;
                break;
            }

            

        }

        
    }

    int ** local_partition_data = new int*[numprocs] ;
    //do the copy
    //1.allocat the space
    for(int i=0;i<numprocs;i++){
        local_partition_data[i] = new int[local_partion_num[i]];
    }
    //2.do the copy
    int * local_count = new int[numprocs];
    for(int i=0;i<numprocs;i++){
        local_count[i] = 0;
    }
    for(int i=0;i<data_per_node;i++){
        
                for(int j=0;j<numprocs-1;j++){
        
                    if(local_data[i]<pivot[j]){
                         local_partition_data[j][ local_count[j]++ ] = local_data[i];
                        // local_partition_data[ local_count[j]++ ] = 100;
                        //  memcpy(local_partition_data[][ local_count[j]++ ],local_data[i],sizeof(int));
                        break;
                    }
        
                    if( local_data[i]>=pivot[j] && j == numprocs -2  ){
                        local_partition_data[j+1][ local_count[j+1]++ ] = local_data[i];
                        break;
                    }
        
                }
        
            }


    //reduce the num the each compute node
    for(int i=0;i<numprocs;i++){
        MPI_Reduce(local_partion_num+i,local_partion_num_sum,1,MPI_INT,MPI_SUM,i,MPI_COMM_WORLD);
    }
    // for(int i=0;i<numprocs;i++){
    //     MPI_Reduce(local_partion_num[myid]+i,local_partion_num_sum,1,MPI_INT,MPI_SUM,i,MPI_COMM_WORLD);
    // }

    //each compute node store all the partition situation
    for(int i=0;i<numprocs;i++){
        MPI_Gather(local_partion_num,numprocs,MPI_INT,local_partion_num_all,numprocs,MPI_INT,i,MPI_COMM_WORLD);
    }

    // for(int i=0;i<numprocs;i++){
    //     MPI_Gather(local_partition_location,numprocs,MPI_INT,local_partition_location_all,numprocs,MPI_INT,i,MPI_COMM_WORLD);
    // }

    //allocat the space for the local partition data
    after_partition_data = new int[local_partion_num_sum[0]];

    //do the partition
    for(int i=0;i<numprocs;i++){
    
    //1.construct the parameter
    int * reccount = new int[numprocs];
    for(int j=0;j<numprocs;j++){
        reccount[j] = local_partion_num_all[ i + j*numprocs ];
    }
    int * displs = new int[numprocs];
    for(int j=0;j<numprocs;j++){
        if(j == 0)
            displs[j] = 0;
        else
            displs[j] = displs[j-1] + reccount[j-1];
    }

    //2.send the data
    MPI_Gatherv(local_partition_data[i],local_partion_num[i],MPI_INT,after_partition_data,reccount,displs,MPI_INT,i,MPI_COMM_WORLD);
    }

    sort(after_partition_data,after_partition_data+local_partion_num_sum[0]);

    // cout<<"  myid"<<myid<<endl;
    // for(int i=0;i<local_partion_num_sum[0];i++){
    //     cout<<after_partition_data[i]<<" ";
    // }
    // cout<<endl;
    
   
    
    // for(int i=0;i<numprocs;i++){
    //     cout<<"  myid"<<myid<<" "<<local_partion_num[i]<<endl;
    // }


    // cout<<"  myid"<<myid<<endl;
    // for(int i=0;i<numprocs*numprocs;i++){
        
    //     cout<<local_partition_location_all[i]<<" ";
    // }
    // cout<<endl;




    //cout<<"  myid"<<myid<<" local_partion_num_sum:"<<local_partion_num_sum[0]<<endl;

    // cout<<"  myid"<<myid<<endl;
    // for(int i=0;i<numprocs*numprocs;i++){
        
    //     cout<<local_partion_num_all[i]<<" ";
    // }
    // cout<<endl;

    // for(int i=0;i<numprocs;i++){
    //     cout<<"  myid"<<myid<<" "<<local_partition_location[i];
    // }cout<<endl;

    // for(int i=0;i<numprocs;i++){
    //     cout<<"  myid"<<myid<<"  "<<i<<"  "<<local_partion_num[myid][i]<<endl;
    // }



    MPI_Barrier(MPI_COMM_WORLD);



    

    // if(myid == 0){
    // for(int i=0;i<numprocs*numprocs;i++){
    //         //cout<<data_per_node<<endl;
    //         cout<<"  myid"<<myid<<" "<<root_sample[i]<<endl;
    //     }
    // }

    //cout the pivot
    // for(int i=0;i<numprocs-1;i++){
    //     cout<<"  myid"<<myid<<" pivot:"<<pivot[i]<<endl;
    // }

    double end_time = MPI_Wtime();
    double diff_time = end_time - start_time;
    if(myid == 0){
        //cout<<"total time is "<<total_time[0]<<" s"<<endl;
        cout<<"total time is "<<diff_time<<" s"<<endl;
    }



    MPI_Finalize();
    return 0;
}