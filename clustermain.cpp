#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include "clusterkmeans.hpp"

using namespace std;

int main(int argc, char const **argv){


	double R = 0.0;
	int query_lines = 0, d = 784, l = 0, K = 0, k = 0, kappa = 0,  M = 0, probes = 0;
	int inputnumberoflines =0;
	char pathname[100], configPathname[100], outputPathName[100];
	char method[50];
	int complete = 0;



	if(argc == 9){ // only filenames specified
				   // in this case the non specified command line arguments are given their default values
	
		for(int i = 1; i < argc; i++){

			if (strcmp(argv[i], "-i") == 0){ 

				strcpy(pathname, argv[++i]);

			}
			else if (strcmp(argv[i], "-c") == 0){ 

				strcpy(configPathname, argv[++i]);
				//query_lines = calculateInputSizefile(queryPathname);
			}
			else if (strcmp(argv[i], "-o") == 0){

				strcpy(outputPathName, argv[++i]);
				//we need it later for writing the result of our programm!
				
			}
			else if (strcmp(argv[i], "-m") == 0){
				strcpy(method, argv[++i]);
			}
			else{

				cout << "Incorrect command line arguments" << endl;
				return -1;
			}

		}

	}
	else if(argc == 10){ // all command line arguments specified

		for (int i = 1; i < argc; i++){

			if (strcmp(argv[i], "-i") == 0){

				strcpy(pathname, argv[++i]);

			}
			else if (strcmp(argv[i], "-c") == 0){
				
				strcpy(configPathname, argv[++i]);
				//query_lines = calculateInputSizefile(queryPathname);
				
			}
			else if (strcmp(argv[i], "-o") == 0){
				
				strcpy(outputPathName, argv[++i]);
				
			}
			else if (strcmp(argv[i], "-complete") == 0){

				complete = 1;
			}
			else if (strcmp(argv[i], "-m") == 0){

				strcpy(method, argv[++i]);
			}
		}
	}
	else{

		cout << "Incorrect number of command line arguments" << endl;
		return -1;
	}

	vector<Point> dataset;
	get_dataset(dataset, pathname, &inputnumberoflines);
	
	ifstream file (configPathname);
	char line[256];
	int count = 0;

	if(file.is_open()){
		cout << "file is open" << endl;
		while(file){
			file.getline(line, 256);
			strtok(line, " ");
			if (count == 0)
			{
				K = stoi(strtok(NULL, " ")); //K for medoids
			}
			else if(count == 1){
				l = stoi(strtok(NULL, " ")); // L for hash tables
			}
			else if(count == 2){
				k = stoi(strtok(NULL, " ")); // k for hash functions lsh
			}
			else if(count == 3){
				M = stoi(strtok(NULL, " ")); // M of hypercube
			}
			else if(count == 4){
				kappa = stoi(strtok(NULL, " ")); // k of hypercube
			}
			else if(count == 5){
				probes = stoi(strtok(NULL, " ")); //probes of hypercube
			}
			count ++;
		}
		
	}
	else{
		cout << "File doesn't open" << endl;
	}
	file.close();
	int m = int(pow(2, (32/k))/2) -1;
	int Mh = int(pow(2, (32/k)));// Mh to M pou theloun oi synarthseis ths class H
	int w = 7500;
	int table_size = inputnumberoflines/16;	
	
	//call cluster query here
	if(strcmp(method, "Lloyds") == 0){
		Lloyd lloyd(dataset, K);
		lloyd.query_cluster(outputPathName, lloyd.Lloydcluster, method, complete);
	}
	else if(strcmp(method, "Range_Search_LSH") == 0){
		Rev_LSH revLSH(dataset, K, m, Mh,  d,  w,  k,  l, table_size, 20);
		revLSH.query_cluster(outputPathName, revLSH.ClusterLsh, method, complete);
	}
	else if(strcmp(method, "Range_Search_Hypercube") == 0){
		Rev_Cube revCube(dataset, K, m, M, w, d, kappa, 20, 10);
		revCube.query_cluster(outputPathName, revCube.ClusterCube, method, complete);
	}
	else{
		cout << "WRONG METHOD!!! " << endl;
	}
	
	return 0;
}