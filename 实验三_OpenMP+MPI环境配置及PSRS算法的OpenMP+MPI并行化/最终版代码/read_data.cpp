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


using namespace std;



int main(int argc, char *argv[]) {


	double * buffer;

	ifstream infile;
	infile.open("../../40G_data.txt", ios::binary);
	if (!infile.is_open()) {
		cout << "error\n";
		exit(0);
	}

	long size = pow(2,20) * 128;

	//buffer = new int[size];
	buffer = (double*)malloc(size*sizeof(double));


	infile.read((char*)buffer, sizeof(double)*size);
	infile.close();

	
	for(int i=0;i<100;i++)
		cout<<buffer[i]<<endl;

}
