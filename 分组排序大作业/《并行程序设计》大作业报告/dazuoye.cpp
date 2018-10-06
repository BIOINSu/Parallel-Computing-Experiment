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
#include<cstring>
#include<ctime>
#include<time.h>
#include<omp.h>
#include <unistd.h>
#include <limits.h>


using namespace std;


#define num_threads 8


void PSRS(int64_t *& array, long n)
{

	int id;
	int i = 0;
	int base = n / num_threads;
	int64_t p[num_threads*num_threads];
	int64_t pivot[num_threads - 1];

	int yu = n % num_threads;
	int last_section_num = base + yu;


	int * pivot_num = new int[num_threads];
	for (int r = 0; r < num_threads; r++) {
		pivot_num[r] = 0;
	}
	int64_t** pivot2 = new int64_t*[num_threads];
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

#pragma omp barrier
		pivot2[id] = (int64_t *)malloc(pivot_num[id] * sizeof(int64_t));

#pragma omp barrier
		if (id != num_threads - 1)
			for (int k = 0; k<base; k++)
			{
				for (int m = 0; m<num_threads - 1; m++)
				{
					if (array[id*base + k] <= pivot[m])
					{
						pivot2[id][count[id]++] = array[id*base + k];
						break;
					}

					if (array[id*base + k] > pivot[m] && m == num_threads - 2) {
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
						pivot2[id][count[id]++] = array[id*base + k];
						break;
					}

					if (array[id*base + k] > pivot[m] && m == num_threads - 2) {
						pivot2[id][count[id]++] = array[id*base + k];
					}
				}
			}


	}

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

	for (int f = 0; f < num_threads; f++) {
		free(pivot2[f]);
	}
	delete[] pivot2;

	i = 0;

}




int64_t* ReadData(string fileName, long length, long offset) {
	
		int64_t *Array = (int64_t*)malloc(length * sizeof(int64_t));
	
		ifstream file;
		file.open(fileName.c_str(), ios::in | ios::binary);
	
		file.seekg(offset * sizeof(int64_t));
		file.read((char*)Array, length * sizeof(int64_t));
		if (!file.is_open()) {
			cout << "open file " << fileName << " fail!" << endl;
			exit(-1);
		}
		file.close();
		cout << "Read Data Complete!" << endl;
		return Array;
	}

void WriteData(int64_t * &Array,string fileName, long length, long offset) {
	
	
		ofstream file;
		file.open(fileName.c_str(), ios::out | ios::binary);

		if (!file.is_open()) {
			cout << "open file " << fileName << " fail!" << endl;
			exit(-1);
		}
	
		file.seekp(offset * sizeof(int64_t));
		file.write((char*)Array, length * sizeof(int64_t));
		
		file.close();
		cout << "Write Data Complete!" << endl;	
	}
	

void MPI_PSRS(int loops) {

	long  data_num = pow(2, 20)*512 ;
	double total_time = 0;

	int myid, numprocs;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	for(int loop = 0 ; loop<loops;loop++)
	{


	int64_t * data = NULL; //store the global data
	long data_per_node = data_num / numprocs;//the number of data for each node to process

	 //store the local data, after partition
	int64_t * local_data = (int64_t*)malloc(data_per_node*sizeof(int64_t));

	//int * local_sample = new int[numprocs]; //store the local sample
	int64_t * local_sample = (int64_t*)malloc(numprocs*sizeof(int64_t));
	int count = 0;

	//int * root_sample = new int[numprocs*numprocs];
	int64_t * root_sample = (int64_t*)malloc(numprocs*numprocs*sizeof(int64_t));

	int64_t * pivot = new int64_t[numprocs - 1];

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


	int64_t * after_partition_data;

	//read data
	local_data = ReadData("../../dazuoye/100G_File.txt", data_per_node, loop*data_num + myid*(data_num/numprocs));

	double start_time = MPI_Wtime();

	//cout<<"can allocate 2"<<endl;

	//sort the local data for each node
	PSRS(local_data, data_per_node);
	//sort(local_data,local_data+data_per_node);
	
	//sample
	for (int i = 0; i<numprocs; i++) {
		//cout<<"myid "+ to_string(myid) +" local sample:" + to_string( local_data[  i * (data_per_node / (numprocs + 1))    ] ) +"\r\n";
		local_sample[i] = local_data[i * (data_per_node / (numprocs + 1))];
	}

	//gather the sample
	MPI_Gather(local_sample, numprocs, MPI_INT64_T, root_sample, numprocs, MPI_INT64_T, 0, MPI_COMM_WORLD);

	//sort the sample in the root and pick the pivot
	if (myid == 0) {

		//sort the local sample
		sort(root_sample, root_sample + numprocs*numprocs);

		//select the pivot
		for (int i = 0; i<numprocs - 1; i++) {
			pivot[i] = root_sample[(i + 1)*(numprocs + 1)];
		}

	}

	//broadcast the pivot to all compute node
	MPI_Bcast(pivot, numprocs - 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

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

	int64_t ** local_partition_data =  (int64_t**)malloc(numprocs*sizeof(int64_t*));   

	//do the copy
	//1.allocat the space
	for (int i = 0; i<numprocs; i++) {
		local_partition_data[i] = (int64_t*)malloc(local_partion_num[i] * sizeof(int64_t));
	}

	//2.do the copy
	int * local_count = new int[numprocs];
	for (int i = 0; i<numprocs; i++) {
		local_count[i] = 0;
	}

	for (int i = 0; i<data_per_node; i++) {

		for (int j = 0; j<numprocs - 1; j++) {

			if (local_data[i]<pivot[j]) {
				local_partition_data[j][local_count[j]++] = local_data[i];
				break;
			}

			if (local_data[i] >= pivot[j] && j == numprocs - 2) {
				local_partition_data[j + 1][local_count[j + 1]++] = local_data[i];
				break;
			}

		}

	}

	free(local_data);

	//cout<<"can allocate 5.5"<<endl;
	//reduce the num the each compute node
	for (int i = 0; i<numprocs; i++) {
		MPI_Reduce(&local_partion_num[i], local_partion_num_sum, 1, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
	}

	//each compute node store all the partition situation
	int * local_partion_num_all = (int*)malloc(sizeof(int)*numprocs*numprocs);
	
	
	MPI_Allgather(local_partion_num,numprocs,MPI_INT,local_partion_num_all,numprocs,MPI_INT,MPI_COMM_WORLD);
	
	
	//allocat the space for the local partition data
	after_partition_data = (int64_t*)malloc(sizeof(int64_t)*local_partion_num_sum[0]);	

	int * reccount =  (int*)malloc(sizeof(int)*numprocs);          //new int[numprocs];
	int * displs =  (int*)malloc(sizeof(int)*numprocs);                //new int[numprocs];
	

	for (int j = 0; j<numprocs; j++) {
		reccount[j] = local_partion_num_all[myid + j*numprocs];
	}

	for (int j = 0; j<numprocs; j++) {
		if (j == 0){
			displs[j] = 0;
			
		}
		else{
			displs[j] = displs[j - 1] + reccount[j - 1];
			
		}
			
	}


	//do the partition
	for (int i = 0; i<numprocs; i++) {
		MPI_Gatherv(local_partition_data[i], local_partion_num[i], MPI_INT64_T, after_partition_data, reccount, displs, MPI_INT64_T, i, MPI_COMM_WORLD);
	}


	//cout<<"can allocate 7.1"<<endl;
	free(local_partition_data);

	//cout<<"can allocate 7.2"<<endl;
	PSRS(after_partition_data, local_partion_num_sum[0]);
	//sort(after_partition_data,after_partition_data+local_partion_num_sum[0]);


	//calculate the offset    local_partion_num_all
	long offset = 0;
	for(int i=0;i<myid;i++){
		for(int j=0;j<numprocs;j++){
			offset = offset + local_partion_num_all[i+j*numprocs];
		}
	}


	WriteData(after_partition_data,to_string((long long)loop)+".txt", local_partion_num_sum[0], offset); //"../../dazuoye/100G_File.txt", data_per_node, myid*(data_num/numprocs)


	ofstream resultfile;
	resultfile.open((to_string( ((long long)loop))+"result"+to_string( ((long long)myid))).c_str(), ios::out | ios::binary);
	int64_t finalcount = 1;
	//iterate the after partition data
	for(int i=0;i<local_partion_num_sum[0];i++){
		
		if(after_partition_data[i] != after_partition_data[i+1]){
			resultfile<<after_partition_data[i]<<finalcount;
			finalcount = 1;
		}else{
			finalcount++;
		}
	}

	
	//close the resultfile stream
	resultfile.close();

	//free all the allocated source
	free(local_sample);
	free(root_sample);
	free(local_partion_num_all);
	free(after_partition_data);
	free(reccount);
	free(displs);

	double end_time = MPI_Wtime();
	double diff_time = end_time - start_time;

	total_time = total_time + diff_time;

	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0) {
		cout << "total time is " << diff_time << " s" << endl;
	}

	}

	if (myid == 0) {
		cout << "All time is " << total_time << " s" << endl;
	}

	MPI_Finalize();
}

void m_mergesort(int loop){

	
	//declare the final result file stream
	ofstream ofile;
	ofile.open("result.txt",ios::out | ios::binary ); 

	//declare the array to store 
	int64_t* data = new int64_t[loop];
	int64_t* read_buff = new int64_t[1];

	//declare the array to store the number of each strand
	long* num_count = new long[loop];
	for(int i=0;i<loop;i++){
		num_count[i] = 0;
	}

	//declare the var to store all the sum of the data num
	long sum_num = 0;

	//declare the input file stream
	ifstream*  infile = new ifstream[loop];
	//open the file stream
	for(int i=0;i<loop;i++){

		infile[i].open( to_string( ((long long)i) )+".txt" , ios::in | ios::binary);

		//seekp move the pointer to the right place
		infile[i].seekg(num_count[i] * sizeof(int64_t));
		
		//read the first data of each strand
		infile[i].read((char*)read_buff, 1 * sizeof(int64_t));
		data[i] = read_buff[0];

		//add the count
		num_count[i]++;
	}

	int64_t min;
	int min_index;
	int64_t* write_buff = new int64_t[1];
	//begin the loop
	while(sum_num != pow(2, 20)*loop){

		//find the min value in data array
		min = data[0];
		min_index = 0;
		for(int i=0;i<loop;i++){
			if(data[i]<min){
				min = data[i];
				min_index = i;
			}
		}
		write_buff[0] = min;

		//write the min value into the result file
		ofile.seekp(sum_num * sizeof(int64_t));

		ofile.write((char*)write_buff, 1 * sizeof(int64_t));
		sum_num++;


		//if the file is end do not read data anymore
		if(num_count[min_index]==pow(2, 10)){
			
			data[min_index] = INT_MAX;

		}else{
		//read the data into buffer
		//seekp move the pointer to the right place
		infile[min_index].seekg(num_count[min_index] * sizeof(int64_t));

		//read the first data of each strand
		infile[min_index].read((char*)read_buff, 1 * sizeof(int64_t));
		data[min_index] = read_buff[0];

		//add the count
		num_count[min_index]++;
		}

	}

	//close the file stream
	for(int i=0;i<loop;i++){
		infile[i].close();	
	}
	ofile.close();

}

int main(){

	int loop = 4;

	//sort the data and store them 
	 MPI_PSRS(loop);

	//do the m_mergesort and calculate the number
	m_mergesort(loop);

}
