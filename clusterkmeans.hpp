#ifndef CLUSTERKMEANS_HPP_
#define CLUSTERKMEANS_HPP_
#include <iostream>
#include <random>
#include "lsh.hpp"
#include "cube.hpp"

using namespace std;

class Cluster{
	protected:
		vector<Point> dataset;
		vector<Point> centroids;
		double clusterTime;
		int k;
	public:
		Cluster(){}
		Cluster(vector<Point>, int); // Initializes centroids with kmeans++ method.
		~Cluster(){}
		double getClusterTime(){return this->clusterTime;}
		vector<double> Silhouette(unordered_map<int, vector<Point>>);
		void query_cluster(char*, unordered_map<int, vector<Point>>, char*, int);
};

class Lloyd: public Cluster{	
	public:
		unordered_map<int, vector<Point>> Lloydcluster;
		Lloyd(vector<Point>, int);
		~Lloyd(){}
		
};

class Rev_LSH: public Cluster{
	LSH lsh;
	public:
		unordered_map<int, vector<Point>> ClusterLsh;
		Rev_LSH(vector<Point>, int, int, int, int, int, int, int, int, int);
		~Rev_LSH(){}
		
};

class Rev_Cube: public Cluster{
	Cube hypercube;
	public:
		unordered_map<int, vector<Point>> ClusterCube;
		Rev_Cube(vector<Point>, int,  int, int, int, int, int, int, int);
		~Rev_Cube(){}
};

bool compare (int, int );

#endif