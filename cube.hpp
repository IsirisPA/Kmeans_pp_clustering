#ifndef CUBE_HPP_
#define CUBE_HPP_
#include <iostream>
#include <unordered_map>
#include "utilities.hpp"

class F{
	int m, M, d, w, dtonos; 			//orismata tou h
	H hi;
	unordered_map<int, int> fi;
	public:
		F(int, int, int, int, int);
		~F(){}
		unordered_map<int, int> getfi(){return this->fi;}
		int getm(){ return this->m;}
		int getd(){ return this->d;}
		int getM(){ return this->M;}
		int getw(){ return this->w;}
		int getdtonos(){ return this->dtonos;}
		int ffunction(Point );

};

class Cube{
	int m, M, d, w, dtonos; 
	vector<Point> dataset;
	vector<F> fresults;
	vector<long int> keys;
	public:
		unordered_map<long int, vector<Point>> cube;// is defined as public so that we can have direct access at the points stored in the HYper Cube by the clustering algorithm
		Cube(vector<Point>, int,  int, int, int, int);
		~Cube(){}
		int getm(){ return this->m;}
		int getd(){ return this->d;}
		int getM(){ return this->M;}
		int getw(){ return this->w;}
		int getdtonos(){ return this->dtonos;}
		vector<F> get_fresults(){ return this->fresults; }
		void print_cube();
		void RangeSearch(vector<Point>&, Point, double, int, int);
		void ANN(vector<NN>&, Point, int, int, int);
		void query_cube(char* , char*, int, double, int, int);
};

int hamming_distance(long int  , long int );
#endif