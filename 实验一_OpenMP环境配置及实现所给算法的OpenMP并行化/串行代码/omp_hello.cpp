#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

void hello();

int main(int argc,char **argv)
{
	int thread_count=atoi(argv[1]);

#	pragma omp parallel num_threads(thread_count)
	hello();

	return 0;
}

void hello()
{
	int my_rank=omp_get_thread_num();
	int thread_count=omp_get_num_threads();

	printf("hello from rank %d of %d \n",my_rank,thread_count);

}