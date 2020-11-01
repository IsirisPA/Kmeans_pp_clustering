#include <iostream>
#include <random>
#include <bits/stdc++.h>
#include <ctime> 
#include <algorithm>
#include <string>
#include "lsh.hpp"
#include "utilities.hpp"


LSH::LSH(vector<Point> dataset, int m, int M, int d, int w, int k, int l, int table_size){

	this->d = d;
	this->w = w;
	this->m = m;
	this->M = M;
	this->k = k;
	this->l = l;	
	this->dataset = dataset;
	this->table_size = table_size;

	for (int i = 0; i < l; i++)
	{
		G g(d, w, m, M, k);
		gi.push_back(g);
	}

	this->hash_tables = new vector<Point>* [l];
	for(int i = 0; i < l ;i++){
			this->hash_tables[i] = new vector<Point> [table_size];
	}

   //unsigned int position = 0; //key value 
   for(int i = 0; i < l; i++){
   		cout << "Creating Hash table " << i << " ....." << endl;
		for(int j = 0; j < dataset.size(); j++){

			int position = (gi.at(i).gfunction(dataset.at(j)))%table_size;
			hash_tables[i][position].push_back(dataset.at(j));
		}
		cout << "Hash table ready!" << endl;
	}

}

void LSH::print_LSH(){

	for(int i = 0; i < l; i++){ //for each hash table

		cout << "Printing hash table at " << i << endl;

		for(int j = 0; j < table_size; j++){ //for each bucket

			cout << "Bucket at " << j << endl;
			//cout << "size" << hash_tables[i][j].size() << endl;

			for(int k = 0; k < hash_tables[i][j].size(); k++){
				cout << "Point id at " << k << " " ;
				cout << hash_tables[i][j].at(k).get_id() << endl;
			}
		}
	}
}


void LSH::RangeSearch(vector<Point>& neighbours, Point q, double r){
	
	unsigned int position = 0;
	double distance = 0.0;
	
	vector<int> range_point = q.getpoint();

	int bucket_size = 0;

	for(int i=0; i < l; i++){


		position = (gi.at(i).gfunction(q))%table_size; 

		bucket_size = hash_tables[i][position].size();
		for(int j=0; j < bucket_size; j++){

			distance = manhattanDistance(range_point, hash_tables[i][position].at(j).getpoint());
			//searching for NN within r
			if (distance < r){
				neighbours.push_back(hash_tables[i][position].at(j));	//store NNs in a vector
			}
			if (neighbours.size() > 20*l){
				return ;
			}
		}
	}
}


void LSH::ANN(vector<NN>& neighbours, Point q, int N){

	unsigned int position = 0;
	double distance = 0.0;
	double bestDistance = numeric_limits<double>::infinity();
	vector<int> range_point = q.getpoint();

	NN x;

	vector <NN> kNN; // krataei ton geitona kai tin apostasi toy apo to query

	int bucket_size = 0;

	for(int i=0; i < l; i++){

		position = (gi.at(i).gfunction(q))%table_size;

		bucket_size = hash_tables[i][position].size();

		for(int j=0; j < bucket_size; j++){

			distance = manhattanDistance(range_point, hash_tables[i][position].at(j).getpoint());
			if (distance < bestDistance){
				bestDistance = distance;
				x.point = hash_tables[i][position].at(j);
				x.distance = distance;
				kNN.push_back(x);
			}
			
			if (kNN.size() > 10*l){		//optional upper bound, can be omitted

				sort(kNN.begin(), kNN.end(),myfunction);	// sort gia na einai se ayksousa seira na pairnoyme ta N me tin mikroteri apostash
				// an size > N epilegoyme mono ta N me tin mikroteri apostasi
				if(kNN.size() >= N){

					for(int t = 0; t < N; t++){
						neighbours.push_back(kNN.at(t));
						//cout << "id = " <<kNN.at(t).point.get_id() << " distance= " << kNN.at(t).distance << endl;
					}

				}
				else{		// alliws epistrefoyme osa vroume

					for(int t = 0; t < kNN.size(); t++){
						neighbours.push_back(kNN.at(t));
						//cout << "id = " <<kNN.at(t).point.get_id() << " distance= " << kNN.at(t).distance << endl;
					}
				}

				
				return ;
			}
		}
	}

	// auto to kommati codika tha ektelesti mono gia knn.size < 10*l
	sort(kNN.begin(), kNN.end(),myfunction); // sort gia na einai se ayksousa seira na pairnoyme ta N me tin mikroteri apostash

	if(kNN.size() >= N){
		// an size > N epilegoyme mono ta N me tin mikroteri apostasi
		for(int t = 0; t < N; t++){
			neighbours.push_back(kNN.at(t));
			//cout << "id = " <<kNN.at(t).point.get_id() << " distance= " << kNN.at(t).distance << endl;
		}

	}
	else{	// alliws epistrefoyme osa vroume

		for(int t = 0; t < kNN.size(); t++){
			neighbours.push_back(kNN.at(t));
			//cout << "id = " <<kNN.at(t).point.get_id() << " distance= " << kNN.at(t).distance << endl;
		}

	}

}


int estimateW(vector<Point>& dataset, int size){
	double sum = 0.0;
	int w = 0;
	Point p;
	for (int i = 0; i < size; i++)
	{
		sum += ExactNN(dataset, dataset.at(i).getpoint(), p, size);
		cout << sum << endl;

	}
	w = int(sum/size);
	return w;
}


void LSH::query_lsh(char* output_file, char* query_file, int N, double R){
	ifstream qfile;
	ofstream ofile;

	vector<Point> queries;
	int number_of_queries = 0;

	get_dataset(queries , query_file, &number_of_queries);
	qfile.open(query_file, ios::binary);
	ofile.open(output_file, ios::out);
	if (ofile.is_open()){cout << "is open\n";}
	else {cout << "is closed\n";}
	int query_size = queries.size();
	Point x;

	for (int i = 0; i < query_size; i++)
	{
	    ofile << "Query: " << i + 1 << endl;
	    ofile.flush();
		vector<NN> neighbours;
		vector<Point> Rneighbours;
		double distanceTrue;

		auto t1 = chrono::steady_clock::now();

		ANN(neighbours, queries.at(i), N);
		auto t2 = chrono::steady_clock::now();
		auto TimeResult = t2 - t1;
		auto time_spanLsh = (chrono::duration<double, milli>(TimeResult).count()) / 1000;
		
		t1 = chrono::steady_clock::now();
		distanceTrue = ExactNN(dataset, queries.at(i), x, dataset.size());
		t2 = chrono::steady_clock::now();
		TimeResult = t2 - t1;
		auto time_spanTrue = (chrono::duration<double, milli>(TimeResult).count()) / 1000;
		for (int i = 0; i < neighbours.size(); i++)
		{
			ofile << "Nearest neighbor-" << i << ": " << neighbours.at(i).point.get_id() << endl;
			ofile << "distanceLSH: " << neighbours.at(i).distance << endl;
			ofile << "distanceTrue: " << distanceTrue << endl;
		}

		ofile << "tLSH: " << time_spanLsh << endl;
		ofile << "tTrue: " << time_spanTrue << endl;
		
		RangeSearch(Rneighbours, queries.at(i), R);
		ofile << "R-near neighbours:" << endl;
		for (int i = 0; i < Rneighbours.size(); i++)
		{
			ofile << Rneighbours.at(i).get_id() << endl;
		}

		
	}
	qfile.close();
	ofile.close();	
}

/*
Query: image_number_in_query_set
Nearest neighbor-1: image_number_in_data_set
distanceLSH: <double> [Î® distanceHypercube Î±Î½Ï„Î¯ÏÏ„Î¿Î¹Ï‡Î±]
distanceTrue: <double>
...
Nearest neighbor-N: image_number_in_data_set
distanceLSH: <double> [Î® distanceHypercube Î±Î½Ï„Î¯ÏÏ„Î¿Î¹Ï‡Î±]
distanceTrue: <double>
tLSH: <double>
tTrue: <double>
R-near neighbors:
image_number_A
image_number_B
. . .
image_number_Z 
*/
	