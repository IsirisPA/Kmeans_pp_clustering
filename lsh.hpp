#ifndef LSH_H_
#define LSH_H_
#include <iostream>
#include <vector>
#include <unordered_map>
#include "utilities.hpp"

using namespace std;



class LSH{
	int m, M, d, w, k, l, table_size;
	vector<G> gi;
	vector<Point> dataset;
	public:
		vector<Point> ** hash_tables;
		LSH(vector<Point>,int , int , int , int, int, int, int);
		~LSH(){
			for(int i = 0; i < l; i++)
				delete [] hash_tables[i];
			delete [] hash_tables;

		}
		int getm(){ return this->m;}
		int getd(){ return this->d;}
		int getM(){ return this->M;}
		int getw(){ return this->w;}
		int getk(){ return this->k;}
		int getl(){ return this->l;}
		vector<G> getgi(){return this->gi;}
		int get_table_size(){ return this->table_size;}
		void print_LSH();
		void RangeSearch(vector<Point>&, Point, double);
		void ANN(vector<NN>&, Point, int);
		void query_lsh(char* , char*, int, double);
};
		


int estimateW(vector<Point>&, int);

#endif