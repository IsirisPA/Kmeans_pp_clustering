#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include "cube.hpp"
#include "utilities.hpp"

using namespace std;

int main(int argc, char const **argv){


	double R = 0.0;
	int query_lines = 0, d = 784, k = 0, N = 0,  M = 0, probes;
	int inputnumberoflines =0;
	char pathname[100], queryPathname[100], outputPathName[100];



	if(argc == 7){ // only filenames specified
				   // in this case the non specified command line arguments are given their default values
		k = 14;
		M = 10;
		N = 1;
		R = 10000.0;
		probes = 2;

		for(int i = 1; i < argc; i++){

			if (strcmp(argv[i], "-d") == 0){ 

				strcpy(pathname, argv[++i]);
			}
			else if (strcmp(argv[i], "-q") == 0){ 

				strcpy(queryPathname, argv[++i]);
			}
			else if (strcmp(argv[i], "-o") == 0){

				strcpy(outputPathName, argv[++i]);
				//we need it for writing the result of our programm!	
			}
			else{

				cout << "Incorrect command line arguments" << endl;
				return -1;
			}
		}
	}
	else if(argc == 17){ // all command line arguments specified

		for (int i = 1; i < argc; i++){

			if (strcmp(argv[i], "-d") == 0){

				strcpy(pathname, argv[++i]);
			}
			else if (strcmp(argv[i], "-q") == 0){
				
				strcpy(queryPathname, argv[++i]);	
			}
			else if (strcmp(argv[i], "-o") == 0){
				
				strcpy(outputPathName, argv[++i]);	
			}
			else if (strcmp(argv[i], "-k") == 0){

				k = atoi(argv[++i]);	
			}
			else if (strcmp(argv[i], "-M") == 0){

				M = atoi(argv[++i]);				
			}
			else if (strcmp(argv[i], "-probes") == 0){

				probes = atoi(argv[++i]);		
			}
			else if (strcmp(argv[i], "-N") == 0){

				N = atoi(argv[++i]);			
			}
			else if (strcmp(argv[i], "-R") == 0){
				
				R = atof(argv[++i]);		
			}
			else{

				cout << "Incorrect command line arguments" << endl;
				return -1;
			}
		}
	}
	else{

		cout << "Incorrect number of command line arguments" << endl;
		return -1;
	}

	/*cout << "R = " << R << endl;
	cout << "k = " << k << endl;
	cout << "N = " << N << endl;
	cout << "M = " << M << endl;
	cout << "input file = " << pathname << endl;
	cout << "query file = " << queryPathname << endl;
	cout << "output file = " << outputPathName << endl;*/

	vector<Point> dataset;
	get_dataset(dataset, pathname, &inputnumberoflines); 
	
	int m = int(pow(2, (32/k))/2) - 1;
	int Mh = int(pow(2, (32/k)));// Mh to M pou theloun oi synarthseis ths class H
	int w = 7500;
	int table_size = inputnumberoflines/16;


	Cube cube(dataset, m, Mh, w, d, k); // na kanoume k = dtonos
	cube.query_cube(outputPathName, queryPathname, N, R,  probes,  M);// M = posa points tha exetastoun se kathe bucket

	char ask[10];

	cout << "Continue with another query? (Please press yes or no)" << endl;
	cin >> ask;
	while(strcmp(ask, "yes") == 0){
		cout << "Please type new query file name" << endl;
		cin >> queryPathname;
		cube.query_cube(outputPathName, queryPathname, N, R,  probes,  M);
		cout << "Continue with another query? (Please press yes or no)" << endl;
		cin >> ask;
	}
	
	return 0;
}