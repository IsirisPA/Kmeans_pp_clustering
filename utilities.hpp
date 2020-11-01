#ifndef UTILITIES_H_
#define UTILITIES_H_
#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <string>
using namespace std;


class Point
{
	vector<int> p;
	int flag;
	int id;

public:
	Point(){}
	Point(vector <int> p){ 
		static int object_count;
		this->p = p;
		this->id = ++object_count;
		this->flag = -1;
	}
	~Point(){this->p.clear();}
	vector<int> getpoint(){ return this->p;}
	int get_id(){return this->id;}
	int get_flag(){return this-> flag;}
	void set_flag(int flag){this-> flag = flag;}
	void print_point(){ 
		for(int i=0; i< p.size(); i++)
			cout << p.at(i) << " ";
		cout << endl;
	}
};

class H{
	//vector<int> a; 
	vector<int> s;
	int m, M, d, w;
	int modular_exponentiation(int ,  int , int );
	public:
		H(){}
		H(int, int , int , int );
		~H(){this->s.clear();}
		vector<int> getSi(){ return this->s;}
		int getm(){ return this->m;}
		int getd(){ return this->d;}
		int getM(){ return this->M;}
		int getw(){ return this->w;}
		unsigned int hfunction(Point );
};

class G{

	vector<H> hi;
	int m, M, d, w, k;
	public:
		G(int, int , int , int , int );
		~G(){this->hi.clear();}
		vector<H> getHi(){ return this->hi;}
		int getm(){ return this->m;}
		int getd(){ return this->d;}
		int getM(){ return this->M;}
		int getw(){ return this->w;}
		int getk(){ return this->k;}
		unsigned int gfunction(Point );

};
// krataei ton geitona kai tin apostasi tou apo to query
//used at ANN <-- LSH
struct NN{
	Point point;
	int distance;
};

double manhattanDistance(vector<int> , vector<int> );
void get_dataset(vector<Point>& , char* , int* );
bool myfunction (NN ,NN );
double ExactNN(vector<Point>, Point, Point&, int);

#endif