#include <iostream>
#include <random>
#include <bits/stdc++.h>
#include <ctime> 
#include <algorithm>
#include <cmath>
#include "utilities.hpp"


bool myfunction (NN i, NN j) { return (i.distance<j.distance); }


double manhattanDistance(vector<int> vector1, vector<int> vector2){
	double totalDistance = 0;
	if (vector1.size() != vector2.size())
	{
		cout << "CAUTION!!!, NOT SAME SIZE OF VECTORS" << endl;
		return 0;
	}
	for (int i = 1; i < vector1.size(); i++)
	{
		totalDistance += fabs(vector1.at(i) - vector2.at(i));
	}
	return totalDistance;
}

/*int modFunction(int a, int mod_b)
{
    int result = a % mod_b;
    if (result < 0)
        result += mod_b;
    return result;
}*/

//h sinartisi auti mas epitrepei na min au3anetai o ari8mos ka8ws to ipsonoume se dinami
int H::modular_exponentiation(int vasi, int ekthetis, int modwith)
{
    unsigned int result = 0;
    long number = 0;
    // elegxos vasis
    if (vasi == 0)
        return 0;
    if (ekthetis == 0)
        return 1;

    if (ekthetis % 2 == 0) // an o ekthetis einai artios
    {
        number = modular_exponentiation(vasi, ekthetis / 2, modwith);
        number = (number * number) % modwith;
    }
    else // an o ekthetis einai perittos
    {
        number = vasi % modwith;
        number = (number * modular_exponentiation(vasi, ekthetis - 1, modwith) % modwith) % modwith;
    }
    result = (number + modwith) % modwith;
    return result;
}



H::H(int d, int w, int m, int M){
	vector<int> s;
	this->d = d;
	this->w = w;
	this->m = m;
	this->M = M;
	random_device rand_dev;									// s uniform random generation
	mt19937	generator(rand_dev());
	uniform_int_distribution<int> distr(0, w-1);
	//generate s
	for (int i = 0; i < d; i++)
	{
		s.push_back(distr(generator));
		//cout << s.at(i)<<endl;
	}
   this->s = s;
  // cout << this->s.size() << endl;

}


//hash function
unsigned int H::hfunction(Point x){

	vector<int> a;
	vector<int> xi;
	xi = x.Point::getpoint();
	unsigned int hvalue = 0;
	
   
    for (unsigned int j = 0; j < d; j++){
    	a.push_back( floor((double)(xi.at(j) - s.at(j))/w)) ;
    	//cout <<  "x = " << xi.at(j) << "s = " << s.at(j) << "a = " << a.at(j) << endl;
    }


    for (unsigned int i = 0; i < d; i++)
    {
    	unsigned int temp = a.at(d-i-1)%M;
    	if( temp  <  0)
    		temp += M;
    	hvalue += (temp * modular_exponentiation(m, i ,M))%M;


    }

    hvalue %= M;

 
    return hvalue;
	
}
// function for bruteforce Nearest Neighbors
double ExactNN(vector<Point> dataset, Point q, Point& neighbour,int size){
	double minimunDistance = manhattanDistance(dataset.at(0).getpoint(), q.getpoint());
	double distance = 0;
	for(int i=1; i < size; i++){
		distance = manhattanDistance(dataset.at(i).getpoint(), q.getpoint());
		
		if(distance < minimunDistance && distance !=0){
			minimunDistance = distance;
			neighbour = dataset.at(i);
		}
	}

	return minimunDistance;
}


G::G(int d, int w, int m, int M, int k){

	vector<H> hi;
	this->d = d;
	this->w = w;
	this->m = m;
	this->M = M;
	this->k = k;
	//generate h
	for (int i = 0; i < k; i++)
	{
		H h(d, w, m, M);
		hi.push_back(h);
	}
   this->hi = hi;


}

unsigned int G::gfunction(Point x){

	vector<unsigned int> hvalues;
	unsigned int gvalue = 0;
	vector<int> p;
	p = x.Point::getpoint();
	for(int i = 0; i < k; i++){
		H obj;
		obj = hi.at(i);
		hvalues.push_back(obj.hfunction(p));
	}
	for(int i = 0; i < k; i++){
		gvalue |= hvalues.at(i)<<(i*32/k);
	}

	return gvalue;
}


int ReverseInt (int i){ 

// Function reversing the binary words stored in the file in little-endian format
// since the data are processed in a big-endian CPU

    unsigned char ch1, ch2, ch3, ch4;

    ch1 = i & 255;
    ch2 = (i >> 8) & 255;
    ch3 = (i >> 16) & 255;
    ch4 = (i >> 24) & 255;

    return((int)ch1 << 24) + ((int) ch2 << 16) + ((int)ch3 << 8) + ch4;
}

void get_dataset(vector<Point>& dataset, char* input_file, int* number_of_points){


	 ifstream file (input_file,ios::binary);

	 if (file.is_open()){

	 	int magic_number = 0; // not used anywhere in the program
	    int number_of_images = 0;
        int n_rows = 0;
        int n_cols = 0;


        file.read((char*)&magic_number,sizeof(magic_number));
        magic_number= ReverseInt(magic_number);

        file.read((char*)&number_of_images,sizeof(number_of_images));
        number_of_images= ReverseInt(number_of_images);

        file.read((char*)&n_rows,sizeof(n_rows));
        n_rows= ReverseInt(n_rows);

        file.read((char*)&n_cols,sizeof(n_cols));
        n_cols= ReverseInt(n_cols);

        *number_of_points = number_of_images;

        for(int i=0;i<number_of_images;++i){ //for every image in the file

        	vector<int> image; // stores the overall pixels (n_rows x n_col) of each image

            for(int r=0;r<n_rows;++r){

                for(int c=0;c<n_cols;++c){

                	unsigned char temp=0;
                    file.read((char*)&temp,sizeof(temp));
                    image.push_back((int)temp);
                }
            }

            dataset.push_back(Point(image)); 
            // a class Point object is constructed, initialized with the above vector 
            // and stored in the  vector representing our dataset 

        }

      
	 }else{ cout << "Error while reading file" << endl; }

}