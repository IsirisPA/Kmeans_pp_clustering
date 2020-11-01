#include <iostream>
#include <random>
#include "cube.hpp"
#include <ctime> 
#include <algorithm>
#include <chrono>


F::F(int m, int M, int w, int d, int dtonos){
	this->d = d;
	this->w = w;
	this->m = m;
	this->M = M;
	this->dtonos = dtonos;
	H hi(d, w, m, M);
	this->hi = hi;

}

int F::ffunction(Point p){

	unordered_map<int, int>::iterator flag;
	random_device rand_dev;									// uniform random generation
	mt19937	generator(rand_dev());
	uniform_int_distribution<int> distr(0, 1);
	int hvalue = 0; 
	hvalue = hi.hfunction(p); 
	flag = fi.find(hvalue); // if hvalue already mapped to {0,1} there is no need to map it again by flipping a coin
	if (flag ==	fi.end()){ // if value not mapped to {0,1}
		fi.insert({hvalue, distr(generator)}); //by flipping a coin
	}

	return fi[hvalue];
}


Cube::Cube(vector<Point> dataset,int m, int M, int w, int d, int dtonos){
	this->d = d;
	this->w = w;
	this->m = m;
	this->M = M;
	this->dtonos = dtonos;
	this->dataset = dataset;


	for(int i =0; i < dtonos; i++){

		F f(m, M, w, d, dtonos);
		fresults.push_back(f);
	}
			
	for(int i = 0; i < dataset.size(); i++){ // for each point
		
		long int key = 0;

		for(int j = 0; j < dtonos; j++){ // construct key

			key = key*10 + fresults.at(j).ffunction(dataset.at(i));
		}

		unordered_map<long int, vector<Point>>::iterator flag;
		flag = cube.find(key);

		if (flag ==	cube.end()){
			vector<Point> bucket;
			bucket.push_back(dataset.at(i));
			cube.insert({key, bucket});
		}
		else{		
			cube[key].push_back(dataset.at(i));
			keys.push_back(key);
		
		}

	}

}

void Cube::print_cube(){

	for (auto& x: cube) {

    	cout << "Key = " << x.first << " ";

    	for(int i = 0; i < x.second.size(); i++)
    		cout <<  x.second.at(i).get_id() << " ";

    	cout << endl;
  	}
}

void Cube::RangeSearch(vector<Point>& neighbours, Point q, double r, int probes, int M){
	
	double distance = 0.0;
	
	vector<int> range_point = q.getpoint();

	long int key = 0;

	for(int i = 0; i < dtonos; i++){ // construct key

		key = key*10 + fresults.at(i).ffunction(q);

	}
	
	int bucket_size = cube[key].size();
	for(int j=0; j < bucket_size; j++){

		distance = manhattanDistance(range_point, cube[key].at(j).getpoint());

				
		if (distance < r){

			neighbours.push_back(cube[key].at(j));
		}
		if (j > M){			// se kathe geitoniko probe koitaei to polu M geitones
			break;
		}
	}
	if(probes > 0){// we want to visit at least one neighbouring vertex.

		int humming_dist = 1;
		int count = 0; // how many neighbouring vertices we've already visited
		int keys_size = keys.size();
		int key_length =  floor(log10(key) + 1); // number of digits of the key

		while(count != probes && humming_dist < key_length ){

			for(int i = 0; i < keys_size; i++){

				if(hamming_distance(key, keys.at(i)) == humming_dist && count != probes){// we found a neighbouring vertex

					int bucket_size = cube[keys.at(i)].size();
					for(int j=0; j < bucket_size; j++){

						distance = manhattanDistance(range_point, cube[keys.at(i)].at(j).getpoint());

						if (distance < r){
							neighbours.push_back(cube[keys.at(i)].at(j));
						}
						if (j > M){	  // se kathe geitoniko probe koitaei to polu M geitones
							break;
						}
					}
					count++;
				}
			}
			humming_dist++;
		}		
	}	
}

void Cube::ANN(vector<NN>& neighbours, Point q, int N, int probes, int M){

	double distance = 0.0;
	double bestDistance = numeric_limits<double>::infinity();
	vector<int> range_point = q.getpoint();

	long int key = 0;

	for(int i = 0; i < dtonos; i++){ // construct key

		key = key*10 + fresults.at(i).ffunction(q);

	}

	NN x;

	vector <NN> kNN; 

	int bucket_size = cube[key].size();
	for(int j=0; j < bucket_size; j++){

		distance = manhattanDistance(range_point, cube[key].at(j).getpoint());
		if (distance < bestDistance){
				bestDistance = distance;
				x.point = cube[key].at(j);
				x.distance = distance;
				kNN.push_back(x);
		}			
		
		if (j > M){

			break;
		}
	}
	if(probes > 0){// we want to visit at least one neighbouring vertex.

		int humming_dist = 1;
		int count = 0; // how many neighbouring vertices we've already visited
		int keys_size = keys.size();
		int key_length =  floor(log10(key) + 1); // number of digits of the key

		while(count != probes && humming_dist < key_length ){

			for(int i = 0; i < keys_size; i++){

				if(hamming_distance(key, keys.at(i)) == humming_dist && count != probes){// we found a neighbouring vertex

					int bucket_size = cube[keys.at(i)].size();
					for(int j=0; j < bucket_size; j++){

						distance = manhattanDistance(range_point, cube[keys.at(i)].at(j).getpoint());
						if (distance < bestDistance){
							bestDistance = distance;
							x.point = cube[keys.at(i)].at(j);
							x.distance = distance;
							kNN.push_back(x);
						}
			
						if (j > M){
							break;
						}
					}
					count++;
				}
			}
			humming_dist++;
		}
	}

	sort(kNN.begin(), kNN.end(),myfunction); // sort gia na einai se ayksousa seira na pairnoyme ta N me tin mikroteri apostash
	// an size > N epilegoyme mono ta N me tin mikroteri apostasi
	if(kNN.size() >= N){

		for(int t = 0; t < N; t++){
			neighbours.push_back(kNN.at(t));
			//cout << "id = " <<kNN.at(t).point.get_id() << " distance= " << kNN.at(t).distance << endl;
		}
	}
	else{
		// alliws epistrefoyme osa vroume
		for(int t = 0; t < kNN.size(); t++){
			neighbours.push_back(kNN.at(t));
			//cout << "id = " <<kNN.at(t).point.get_id() << " distance= " << kNN.at(t).distance << endl;
		}
	}
}



int hamming_distance(long x, long y)
{
    int dist = 0;
    // Count the number of bits set
    for (unsigned val = x ^ y; val > 0; val = val >> 1)
    {
        // If A bit is set, so increment the count
        if (val & 1)
            dist++;
        // Clear (delete) val's lowest-order bit
    }
    // Return the number of differing bits
    return dist;
}

void Cube::query_cube(char* output_file, char* query_file, int N, double R, int probes, int M){
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

		ANN(neighbours, queries.at(i), N, probes, M);
		auto t2 = chrono::steady_clock::now();
		auto TimeResult = t2 - t1;
		auto time_spanCube = (chrono::duration<double, milli>(TimeResult).count()) / 1000;
		
		t1 = chrono::steady_clock::now();
		distanceTrue = ExactNN(dataset, queries.at(i), x, dataset.size());
		t2 = chrono::steady_clock::now();
		TimeResult = t2 - t1;
		auto time_spanTrue = (chrono::duration<double, milli>(TimeResult).count()) / 1000;
		for (int i = 0; i < neighbours.size(); i++)
		{
			ofile << "Nearest neighbor-" << i << ": " << neighbours.at(i).point.get_id() << endl;
			ofile << "distanceHypercube: " << neighbours.at(i).distance << endl;
			ofile << "distanceTrue: " << distanceTrue << endl;
		}

		ofile << "tHypercube: " << time_spanCube << endl;
		ofile << "tTrue: " << time_spanTrue << endl;
		
		RangeSearch(Rneighbours, queries.at(i), R, probes, M);
		ofile << "R-near neighbours:" << endl;
		for (int i = 0; i < Rneighbours.size(); i++)
		{
			ofile << Rneighbours.at(i).get_id() << endl;
		}

		
	}
	qfile.close();
	ofile.close();	
}

