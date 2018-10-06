//该代码用于生成排序所用的数据，倒序

#include <iostream>
#include <fstream>
#include <string>
#define BLOCK 1024*1024

using namespace std;
int main()
{
	cout << "Enter the Data Length(MB): " << endl;
	int length;
	cin >> length;
	
	int* Array = (int*)malloc(length * BLOCK * sizeof(int));

	for (int i = 0; i < length * BLOCK; i++) {
		Array[i] = 100 - i;
	}

	string fileName = "data" + to_string(length) + "M.txt";
	ofstream newFile;
	newFile.open(fileName, ios::out | ios::binary);
	
	if (!newFile) {
		cout << "Create file " + fileName + " fail!" << endl;
	}
	newFile.write((char*)Array, length * BLOCK * sizeof(int));

	newFile.close();
	free(Array);

	return 0;
}
