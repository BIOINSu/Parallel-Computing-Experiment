#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include<iostream>
#include<cstdlib>
#include<ctime>
#include<vector>
#include<fstream>
#include<algorithm>
#include<math.h>
#include <malloc.h>
#include<string>
#include<ctime>
#include<time.h>
#include<omp.h>
#include <unistd.h>


using namespace std;


#define num_threads 4


void PSRS(double *& array, long n)
{

	int id;
	int i = 0;
	int base = n / num_threads;
	double p[num_threads*num_threads];
	double pivot[num_threads - 1];

	int yu = n % num_threads;
	int last_section_num = base + yu;


	int * pivot_num = new int[num_threads];
	for (int r = 0; r < num_threads; r++) {
		pivot_num[r] = 0;
	}
	double** pivot2 = new double*[num_threads];
	int count[num_threads] = { 0 };



	omp_set_num_threads(num_threads);


#pragma omp parallel shared(base,n,i,p,pivot,array,count,pivot2) private(id)
	{
		id = omp_get_thread_num();

		if (id != num_threads - 1)
			sort(array + id * base, array + (id + 1)*base);
		else
			sort(array + id * base, array + id * base + last_section_num);

#pragma omp critical
		{
			for (int k = 0; k < num_threads; k++) {
				p[i++] = array[id *base + k*(base / num_threads + 1)];
			}

		}


#pragma omp barrier
#pragma omp master
		{
			sort(p, p + i);
			for (int m = 0; m < num_threads - 1; m++) {
				pivot[m] = p[(m + 1)*num_threads];
			}
		}

		//#pragma omp master		
		//{
		//			for (int m = 0; m < num_threads*num_threads; m++) {
		//				cout << "p[" + to_string(m) + "]:" + to_string(p[m]) << endl;
		//			}
		//
		//for (int m = 0; m < num_threads - 1; m++) {
		//	cout << "openmp thread id " + to_string(id) + " pivot[" + to_string(m) + "]:" + to_string(pivot[m]) + "\r\n";
		//}
		//}

#pragma omp barrier

		if (id != num_threads - 1)
			for (int k = 0; k<base; k++)
			{
				for (int m = 0; m<num_threads - 1; m++)
				{
					if (array[id*base + k] <= pivot[m])
					{
						pivot_num[id] = pivot_num[id] + 1;
						break;
					}

					if (array[id*base + k] > pivot[m] && m == num_threads - 2) {
						pivot_num[id] = pivot_num[id] + 1;
						break;
					}
				}
			}
		else
			for (int k = 0; k<last_section_num; k++)
			{
				for (int m = 0; m<num_threads - 1; m++)
				{
					if (array[id*base + k] <= pivot[m])
					{
						pivot_num[id] = pivot_num[id] + 1;
						break;
					}

					if (array[id*base + k] > pivot[m] && m == num_threads - 2) {
						pivot_num[id] = pivot_num[id] + 1;
						break;
					}
				}
			}

		//#pragma omp master
		//cout << " the number is ok " << endl;
		//for (int m = 0; m < num_threads; m++) {
		//	cout << " pivot_num[" + to_string(m) + "]:" + to_string(pivot_num[m]) << endl;
		//}


#pragma omp barrier
		pivot2[id] = (double *)malloc(pivot_num[id] * sizeof(double));



#pragma omp barrier
		if (id != num_threads - 1)
			for (int k = 0; k<base; k++)
			{
				for (int m = 0; m<num_threads - 1; m++)
				{
					if (array[id*base + k] <= pivot[m])
					{
						//if(count[id]<pivot_num[id])
						pivot2[id][count[id]++] = array[id*base + k];
						break;
					}

					if (array[id*base + k] > pivot[m] && m == num_threads - 2) {
						//if (count[id]<pivot_num[id])
						pivot2[id][count[id]++] = array[id*base + k];
					}
				}
			}
		else
			for (int k = 0; k<last_section_num; k++)
			{
				for (int m = 0; m<num_threads - 1; m++)
				{
					if (array[id*base + k] <= pivot[m])
					{
						//if (count[id]<pivot_num[id])
						pivot2[id][count[id]++] = array[id*base + k];
						break;
					}

					if (array[id*base + k] > pivot[m] && m == num_threads - 2) {
						//if (count[id]<pivot_num[id])
						pivot2[id][count[id]++] = array[id*base + k];
					}
				}
			}


	}

#pragma omp master
	//cout << " the partiton is ok " << endl;


#pragma omp parallel shared(pivot2,pivot_num,count) private(id)
	{
		int id = omp_get_thread_num();

		sort(pivot2[id], pivot2[id] + pivot_num[id]);


		int start_position = 0;
		for (int i = 0; i < id + 1; i++) {

			if (i == 0)
				start_position = 0;
			else
				start_position = start_position + pivot_num[i - 1];
		}


		for (int i = 0; i < pivot_num[id]; i++) {
			array[start_position + i] = pivot2[id][i];
		}


	}

#pragma omp master
	//cout << " the sort is ok " << endl;

	for (int f = 0; f < num_threads; f++) {
		free(pivot2[f]);
	}
	delete[] pivot2;


	i = 0;


}



void read(double * &buffer, long size) {

	ifstream infile;
	infile.open("../../40G_data.txt", ios::binary);
	if (!infile.is_open()) {
		cout << "error\n";
		exit(0);
	}

	//buffer = new int[size];
	buffer = (double*)malloc(size*sizeof(double));


	infile.read((char*)buffer, sizeof(double)*size);
	infile.close();


}

double* ReadData(string fileName, long length, long offset) {
	
		double *Array = (double*)malloc(length * sizeof(double));
	
		ifstream file;
		file.open(fileName.c_str(), ios::in | ios::binary);
	
		file.seekg(offset * sizeof(double));
		file.read((char*)Array, length * sizeof(double));
		if (!file.is_open()) {
			cout << "open file " << fileName << " fail!" << endl;
			exit(-1);
		}
		file.close();
		cout << "Read Data Complete!" << endl;
		return Array;
	}
	

int main(int argc, char *argv[]) {



	long  data_num = pow(2, 30)*2 ;
	//string filepath = ;


	int myid, numprocs;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	double * data = NULL; //store the global data
	long data_per_node = data_num / numprocs;//the number of data for each node to process

												 //int * local_data = new int[data_per_node]; //store the local data, after partition
	double * local_data = (double*)malloc(data_per_node*sizeof(double));

	//int * local_sample = new int[numprocs]; //store the local sample
	double * local_sample = (double*)malloc(numprocs*sizeof(double));
	int count = 0;

	//int * root_sample = new int[numprocs*numprocs];
	double * root_sample = (double*)malloc(numprocs*numprocs*sizeof(double));

	double * pivot = new double[numprocs - 1];

	int * local_partion_num = new int[numprocs];
	for (int i = 0; i<numprocs; i++) {
		local_partion_num[i] = 0;
	}
	

	//store the number of data each compute node should process
	int* local_partion_num_sum = new int[1];


	//store the location of each pivot in each compute node
	long* local_partition_location = new long[numprocs];
	for (int i = 0; i<numprocs; i++) {
		local_partition_location[i] = 0;
	}
	long * local_partition_location_all = new long[numprocs*numprocs];


	double * after_partition_data;

	//cout<<"can allocate 1"<<endl;


	//read data
	if(numprocs == 1)
	if (myid == 0) {

		//read data
		read(data, data_num);

		double begin, end, time;
		begin = omp_get_wtime();

		sort(data,data+data_num);

		end = omp_get_wtime();
		time = end - begin;
		printf("The running time is %lfs\n", time);

		//sleep(5000);

		//exit(0);
	}

	local_data = ReadData("../../40G_data.txt", data_per_node, myid*(data_num/numprocs));

	double start_time = MPI_Wtime();

	//cout<<"can allocate 2"<<endl;


	//scatter the data to each compute node 
	//MPI_Scatter(data, data_per_node, MPI_DOUBLE, local_data, data_per_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//cout<<"can allocate 2.1"<<endl;

	free(data);




	//sort the local data for each node
	PSRS(local_data, data_per_node);

	//cout<<"can allocate 3"<<endl;

	//sample
	for (int i = 0; i<numprocs; i++) {
		//cout<<"myid "+ to_string(myid) +" local sample:" + to_string( local_data[  i * (data_per_node / (numprocs + 1))    ] ) +"\r\n";
		local_sample[i] = local_data[i * (data_per_node / (numprocs + 1))];
	}

	//gather the sample
	MPI_Gather(local_sample, numprocs, MPI_DOUBLE, root_sample, numprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//sort the sample in the root and pick the pivot
	if (myid == 0) {

		//sort the local sample
		sort(root_sample, root_sample + numprocs*numprocs);

		//select the pivot
		for (int i = 0; i<numprocs - 1; i++) {
			pivot[i] = root_sample[(i + 1)*(numprocs + 1)];

			//cout << "myid " + to_string(myid) + " pivot:" + to_string(pivot[i]) + "\r\n";
		}

	}

	//cout<<"can allocate 4"<<endl;

	//broadcast the pivot to all compute node
	MPI_Bcast(pivot, numprocs - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);




	//traverse the local_data and record the number for each partition
	for (int i = 0; i<data_per_node; i++) {

		//record the number for each partition
		for (int j = 0; j<numprocs - 1; j++) {

			if (local_data[i]<pivot[j]) {
				local_partion_num[j] = local_partion_num[j] + 1;
				break;
			}

			if (local_data[i] >= pivot[j] && j == numprocs - 2) {
				local_partion_num[j + 1] = local_partion_num[j + 1] + 1;
				break;
			}



		}


	}

	//cout<<"can allocate 5"<<endl;


	double ** local_partition_data =  (double**)malloc(numprocs*sizeof(double*));   // new double*[numprocs];

	//cout<<"can allocate 5.1"<<endl;
	//do the copy
	//1.allocat the space
	for (int i = 0; i<numprocs; i++) {
		//local_partition_data[i] = new int[local_partion_num[i]];
		local_partition_data[i] = (double*)malloc(local_partion_num[i] * sizeof(double));
		//cout<<"myid "<<myid<<" local_partion_num ["<<i<<"] "<<local_partion_num[i]<<endl;
	}

	//cout<<"can allocate 5.2"<<endl;
	//2.do the copy
	int * local_count = new int[numprocs];
	for (int i = 0; i<numprocs; i++) {
		local_count[i] = 0;
	}

	//cout<<"can allocate 5.3"<<endl;
	for (int i = 0; i<data_per_node; i++) {

		for (int j = 0; j<numprocs - 1; j++) {

			if (local_data[i]<pivot[j]) {
				local_partition_data[j][local_count[j]++] = local_data[i];
				// local_partition_data[ local_count[j]++ ] = 100;
				//  memcpy(local_partition_data[][ local_count[j]++ ],local_data[i],sizeof(int));
				break;
			}

			if (local_data[i] >= pivot[j] && j == numprocs - 2) {
				local_partition_data[j + 1][local_count[j + 1]++] = local_data[i];
				break;
			}

		}

	}

	//cout<<"can allocate 5.4"<<endl;
	free(local_data);

	//cout<<"can allocate 5.5"<<endl;
	//reduce the num the each compute node
	for (int i = 0; i<numprocs; i++) {
		MPI_Reduce(&local_partion_num[i], local_partion_num_sum, 1, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
	}

	

	//cout<<"myid "<<myid<<" local_partion_num_sum "<<local_partion_num_sum[0]<<endl;

	//cout<<"can allocate 5.6"<<endl;

	//exit(0);

	//each compute node store all the partition situation
	//int * local_partion_num_all = new int[numprocs*numprocs];
	int * local_partion_num_all = (int*)malloc(sizeof(int)*numprocs*numprocs);
	// for (int i = 0; i<numprocs; i++) {
	// 	MPI_Gather(&local_partion_num, numprocs, MPI_INT, local_partion_num_all, numprocs, MPI_INT, i, MPI_COMM_WORLD);
	// }
	//cout<<"numprocs "<<numprocs<<endl;
	MPI_Allgather(local_partion_num,numprocs,MPI_INT,local_partion_num_all,numprocs,MPI_INT,MPI_COMM_WORLD);
	//cout<<"myid "<<myid<<" local_partion_num_all "<<local_partion_num_all[0]<<endl;

	//cout<<"can allocate 5.7"<<endl;
	//allocat the space for the local partition data
	//after_partition_data = new double[local_partion_num_sum[0]];
	after_partition_data = (double*)malloc(sizeof(double)*local_partion_num_sum[0]);
	//cout<<"myid "<<myid<<" local_partion_num_sum[0]  "<<local_partion_num_sum[0]<<endl;
	
	//cout<<"can allocate 6"<<endl;

	// if(myid == 0)
	// for(int i=0;i<numprocs*numprocs;i++){
	// 	cout<<local_partion_num_all[i]<<" ";
	// 	if(i!=0&&i%3==0)
	// 		cout<<endl;
	// }

	int * reccount =  (int*)malloc(sizeof(int)*numprocs);          //new int[numprocs];
	int * displs =  (int*)malloc(sizeof(int)*numprocs);                //new int[numprocs];
	// for(int i = 0; i<numprocs; i++) {

	// 	for (int j = 0; j<numprocs; j++) {
	// 		reccount[j] = local_partion_num_all[i + j*numprocs];
	// 		if(myid == 0)
	// 			cout<<"myid "<<myid<<" "<<j<<" recount"<<reccount[j]<<endl;
	// 			//cout<<endl;
	// 	}
		
	// 	for (int j = 0; j<numprocs; j++) {
	// 		if (j == 0){
	// 			displs[j] = 0;
	// 			if(myid == 0)
	// 				cout<<"myid "<<myid<<" "<<j<<" displs"<<displs[j]<<endl;
	// 				//cout<<endl;
	// 		}
	// 		else{
	// 			displs[j] = displs[j - 1] + reccount[j - 1];
	// 			if(myid == 0)
	// 				cout<<"myid "<<myid<<" "<<j<<" displs"<<displs[j]<<endl;
	// 				//cout<<endl;
	// 		}
				
	// 	}

	// }

	for (int j = 0; j<numprocs; j++) {
		reccount[j] = local_partion_num_all[myid + j*numprocs];
		//if(myid == 0)
			//cout<<"myid "<<myid<<" "<<j<<" recount"<<reccount[j]<<endl;
			//cout<<endl;
	}

	for (int j = 0; j<numprocs; j++) {
		if (j == 0){
			displs[j] = 0;
			//if(myid == 0)
				//cout<<"myid "<<myid<<" "<<j<<" displs"<<displs[j]<<endl;
				//cout<<endl;
		}
		else{
			displs[j] = displs[j - 1] + reccount[j - 1];
			//if(myid == 0)
				//cout<<"myid "<<myid<<" "<<j<<" displs"<<displs[j]<<endl;
				//cout<<endl;
		}
			
	}


	//do the partition
	for (int i = 0; i<numprocs; i++) {

		//cout<<"can allocate 6.1"<<endl;
		//1.construct the parameter

		//cout<<"can allocate 6.2"<<endl;
		//2.send the data
		// if(local_partition_data[i] == NULL)
		// 	cout<<"local_partition_data[i] == NULL"<<endl;
		// if(local_partion_num[i] == NULL)
		// 	cout<<"local_partion_num[i] == NULL"<<endl;
		// if(reccount == NULL)
		// 	cout<<"reccount == NULL"<<endl;
		// if(displs == NULL)
		// 	cout<<"displs == NULL"<<endl;
		// if(after_partition_data == NULL)
		// 	cout<<"after_partition_data == NULL"<<endl;
		MPI_Gatherv(local_partition_data[i], local_partion_num[i], MPI_DOUBLE, after_partition_data, reccount, displs, MPI_DOUBLE, i, MPI_COMM_WORLD);
	}

	
		
	// cout<<"can allocate 6.1"<<endl;
	// //1.construct the parameter
	// int * reccount = new int[numprocs];
	// for (int j = 0; j<numprocs; j++) {
	// 	reccount[j] = local_partion_num_all[myid + j*numprocs];
	// }
	// int * displs = new int[numprocs];
	// for (int j = 0; j<numprocs; j++) {
	// 	if (j == 0)
	// 		displs[j] = 0;
	// 	else
	// 		displs[j] = displs[j - 1] + reccount[j - 1];
	// }

	// cout<<"can allocate 6.2"<<endl;
	// //2.send the data
	// MPI_Gatherv(local_partition_data[myid], local_partion_num[myid], MPI_DOUBLE, after_partition_data, reccount, displs, MPI_DOUBLE, myid, MPI_COMM_WORLD);





	//cout<<"can allocate 7.1"<<endl;
	free(local_partition_data);

	//cout<<"can allocate 7.2"<<endl;
	PSRS(after_partition_data, local_partion_num_sum[0]);

	//cout<<"can allocate 7.3"<<endl;


	//cout the pivot
	// for (int i = 0; i<numprocs - 1; i++) {
	// 	cout << "myid" << myid << " pivot:" << pivot[i] << endl;
	// }


	double end_time = MPI_Wtime();
	double diff_time = end_time - start_time;

	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0) {
		cout << "total time is " << diff_time << " s" << endl;
	}



	MPI_Finalize();
	return 0;
}
