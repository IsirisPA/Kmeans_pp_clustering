#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include "clusterkmeans.hpp"

bool compare (int i, int j) { return (i < j); }


Cluster::Cluster(vector<Point> dataset, int k){
	//start clock now
	auto t1 = chrono::steady_clock::now();
	int t=0; // number of centroids found
	random_device rand_dev;									// t uniform random generation
	mt19937	generator(rand_dev());
	uniform_real_distribution<> dis(0.0, dataset.size());
	int rand_selection = dis(generator);
	centroids.push_back(dataset.at(rand_selection)); // randomly pick first centroid
	t++;

	this-> k = k; // number of centroids to be initialized
	this->dataset = dataset;

	vector<double> partial_sums;
	vector<double> distances;
	vector<int> centroid_position;
	centroid_position.push_back(rand_selection);

	while(t != k){

		vector<double> probabilities;
		int max_dist = -1;
		for(int i = 0; i < dataset.size(); i++){

			// max distance metaksi 2 centroids
			int distance = manhattanDistance(dataset.at(i).getpoint(), centroids.at(t-1).getpoint()); 
			if(max_dist < distance){
				max_dist = distance;
			}	
		}

		for (int i = 0; i < dataset.size(); i++)
		{
			unsigned int sum = 0;
			int centroid_check = 0;

			for (int j = 0; j < centroids.size(); j++)
			{
				//check if point is centroid do not calculate the distance(manhattan)
				if(dataset.at(i).get_id() == centroids.at(j).get_id()){
					centroid_check = 1;
					break;
				}
			}
			if (centroid_check == 1)
				continue;

			// kanei assigned to point sto centroid pou apexei tin mikroteri apostasi 
			int min_dist = manhattanDistance(dataset.at(i).getpoint(), centroids.at(0).getpoint());
			for (int j = 1; j < centroids.size(); j++)
			{
				int distance = manhattanDistance(dataset.at(i).getpoint(), centroids.at(j).getpoint());

				if(min_dist > distance){
					min_dist = distance;
				}
				
			}
			// upologismos merikon athroismaton
			if(partial_sums.size() == 0){
				sum +=  (min_dist*min_dist)/max_dist;
				partial_sums.push_back(sum);
			}
			else{
				sum = partial_sums.at(partial_sums.size() - 1);
				sum +=  (min_dist*min_dist)/max_dist;
				partial_sums.push_back(sum);
			}
			
		}
	
		// epilegei omoiomorfa ena simeio anamesa sto diastima 0 - P(n-t)
		random_device rand_dev;									
		mt19937	generator(rand_dev());
		uniform_real_distribution<> dis(0.0, partial_sums.at(partial_sums.size() - 1));
		double x = dis(generator);
		//cout << "x:" << x << endl;
		int position = 0;

		for (int i = 0; i < partial_sums.size() - 1; i++)
		{ // vriskoyme anamesa se poia partia sum vrisketai kai dialegei to megalutero
			if(x > partial_sums.at(i) && x <= partial_sums.at(i + 1)){
				position = i + 1;
				if(centroid_position.size() == 1){
					//an exoyme mono 1 centroid 
					if (position < centroid_position.at(0)){ 
						break;
					}
					else{
						position += 1;
						break;
					}


				}else if(centroid_position.size() > 1){
					// an exoume >1 centroids
					if (position < centroid_position.at(0)){
						break;
					}
					else{
						// vriskoume anamesa se poia centroid vrisketai kai posa proigountai
						// kai prosthetoyme ton arithmo ton centroids pou proigountai sto position 
						for (int j = 0; j < centroid_position.size() - 1; j++)
						{
							if (position > centroid_position.at(j) && position <= centroid_position.at(j + 1))
							{
								position += (j + 1);
								break;
							}
						}
					}		
				}
				else{
					cout << "Error no centroids" << endl;
				}					
			}
		}

		centroids.push_back(dataset.at(position));
		centroid_position.push_back(position);
		sort(centroid_position.begin(), centroid_position.end(),compare);
		t++;
		partial_sums.clear();
	}
	auto t2 = chrono::steady_clock::now();
	auto TimeResult = t2 - t1;
	auto time_spanCluster = (chrono::duration<double, milli>(TimeResult).count()) / 1000;

	this ->clusterTime = time_spanCluster;
}

//elegxos an oi pinakes 2 centroid einai akrivws idioi
int configuration(vector<Point> old_centroids, vector<Point> new_centroids){

	if(old_centroids.size() == new_centroids.size()){

		for(int i = 0; i < old_centroids.size(); i++){

			for(int z = 0; z < new_centroids.at(i).getpoint().size(); z++){
				// an diaferoun esto kai se mia suntetagmeni einai diaforetika opote return 0
				if(old_centroids.at(i).getpoint().at(z) != new_centroids.at(i).getpoint().at(z))
					return 0;
			}
		}
		return 1;
	}
}


Lloyd::Lloyd(vector<Point> dataset, int k): Cluster(dataset, k){
	int centroid_size = centroids.size();
	int dataset_size = dataset.size();
	int distance = 0; 
	int position = 0;
	int config_call = 0;
	// Lloyd's assigment
	
	//starting clock
	auto t1 = chrono::steady_clock::now();

	while(config_call == 0) {
		// gia kathe point  
		for (int i = 0; i < dataset_size; i++)
		{	
			int centroid_check = 0;
			// paralupoume ta point poy einai centroids sto dataset
			for (int j = 0; j < centroids.size(); j++)
			{		
				if(dataset.at(i).get_id() == centroids.at(j).get_id()){
					centroid_check = 1;
					break;
				}
			}
			if (centroid_check == 1)
				continue;

			int minimum = manhattanDistance(dataset.at(i).getpoint(), centroids.at(0).getpoint()); 
			position = 0;
			// vriskoyme to pio kontino centroid
			for (int j = 1; j < centroid_size; j++)
			{
				distance = manhattanDistance(dataset.at(i).getpoint(), centroids.at(j).getpoint());
				if (minimum > distance){
					minimum = distance;
					position = j;
				}			
			}
			unordered_map<int, vector<Point>>::iterator flag;
			flag = Lloydcluster.find(centroids.at(position).get_id());
			//to kanoume assign sto cluster toy pio kontinou centroid pou exoyme vrei
			if (flag ==	Lloydcluster.end()){
				vector<Point> bucket;
				bucket.push_back(dataset.at(i));
				Lloydcluster.insert({centroids.at(position).get_id(), bucket});
			}
			else{
				Lloydcluster[centroids.at(position).get_id()].push_back(dataset.at(i));
				
			}
		}
		//endl of Lloyd's assigment

		//update with kmediods for Lloyd
		vector<Point> new_centroids;
		for (int i = 0; i < centroid_size; i++)
		{	//gia kathe centroid 
			
			vector<int> kmedoid;
			vector<Point> bucket = Lloydcluster[centroids.at(i).get_id()];
			if(bucket.size() == 0){
				for (int z = 0; z < centroids.at(i).getpoint().size(); z++)
				{
					kmedoid.push_back(0);
				}
				new_centroids.push_back(Point(kmedoid));
				continue;
			}
			kmedoid = bucket.at(0).getpoint();
			for (int j = 1; j < bucket.size(); j++)
			{	//gia kathe point sto cluster tou centroid
				for (int z = 0; z < bucket.at(j).getpoint().size(); z++)	//ypologismos athroismatos dianusmatwn 
				{
					
					kmedoid.at(z) += bucket.at(j).getpoint().at(z);
				}
			}
			for (int z = 0; z < kmedoid.size(); z++)	//ypologismos mesou orou athroismatos dianusmatwn 
			{
				kmedoid.at(z) /= bucket.size();
			}

			new_centroids.push_back(Point(kmedoid));
		}


		if(configuration(centroids,new_centroids) == 1){
			// optimal centroids found, no change in new centroids
			
			config_call = 1;
		} 
		else{ // update centroids
			for(int i = 0; i < new_centroids.size(); i++){
				centroids.at(i) = new_centroids.at(i);
			}
			Lloydcluster.clear();
		}
	}
	//end clock
	auto t2 = chrono::steady_clock::now();
	auto TimeResult = t2 - t1;
	auto time_spanCluster = (chrono::duration<double, milli>(TimeResult).count()) / 1000;

	this->clusterTime += time_spanCluster;
}


Rev_LSH::Rev_LSH(vector<Point> dataset, int kappa, int m, int M, int d, int w, int k, int l, int table_size, int max_iter = 20) : 
	Cluster(dataset, kappa), lsh(dataset, m, M, d, w, k, l, table_size)
{

	int centroid_size = centroids.size();
	int distance = 0;
	int mindistance;
	int radius = 0;
	int config_call = 0;
	
	unordered_map<int, vector<Point>> assigned_point;
	vector<Point> second_to_last_centroids;
	for(int i = 0; i < centroids.size(); i++){
		second_to_last_centroids.push_back(centroids.at(i));
	}
	int count = 1;
	
	//clock starts
	auto t1 = chrono::steady_clock::now();

	while(config_call == 0){

		for (int j = 0; j < lsh.getl(); j++) // for each hash table in LSH 
		{
			for (int i = 0; i < centroids.size(); i++) // for each centroid
			{
				// hash the centroid to eash LSH hash table and find the bucket it would have normally enter if it was a non centroid point
				int position = (lsh.getgi().at(j).gfunction(centroids.at(i)))%table_size; 

				for (int w = 0; w < lsh.hash_tables[j][position].size(); w++) // for each point in that bucket
				{	
					int ClosestCentroidFlag = 0; // flag = 0 has not been assigned a centroid, flag = 1 has already been visited
				 	    						 // but can be also revisited by another centroid in order to be assigned at that centroid
												 // in case the later centroid is the closest to the point 
					mindistance = manhattanDistance(centroids.at(i).getpoint(), lsh.hash_tables[j][position].at(w).getpoint());
					int min_pos = 0;	
					for (int t = 0; t < centroids.size(); t++) // for each centroid find the one closest to the point
					{
						if(t != i){ //don't compare each centroid with itself
							int positionT = (lsh.getgi().at(j).gfunction(centroids.at(t)))%table_size;
							if(position == positionT){
								distance = manhattanDistance(centroids.at(t).getpoint(), lsh.hash_tables[j][positionT].at(w).getpoint());
								if (mindistance > distance){
									mindistance = distance;
									min_pos = t;
								}
							}	
						}
					}
					// insert barrier
					if(lsh.hash_tables[j][position].at(w).get_flag() == -1){
						lsh.hash_tables[j][position].at(w).set_flag(min_pos);
						
					}

					if(lsh.hash_tables[j][position].at(w).get_flag() != -1){

						unordered_map<int, vector<Point>>::iterator Rev_LSHfind;
							
						Rev_LSHfind = ClusterLsh.find(centroids.at(i).get_id());

						if (Rev_LSHfind ==	ClusterLsh.end()){ //first assignment of a point to the centroid
							
							vector<Point> bucket;
							bucket.push_back(lsh.hash_tables[j][position].at(w));
							ClusterLsh.insert({centroids.at(i).get_id(), bucket});
							}
						else {
							ClusterLsh[centroids.at(i).get_id()].push_back(lsh.hash_tables[j][position].at(w));								
						}

					} 
				}
			}		
		}

		// update centroids for reverse assignement
		vector<Point> new_centroids;
		for (int i = 0; i < centroid_size; i++)
		{	//for every centroid 
				
			vector<int> kmedoid;
			vector<Point> bucket = ClusterLsh[centroids.at(i).get_id()];
			if(bucket.size() == 0){  // initialize bucket if empty
				for (int z = 0; z < centroids.at(i).getpoint().size(); z++)
				{
					kmedoid.push_back(0);
				}
				new_centroids.push_back(Point(kmedoid));
				continue;
			}
			kmedoid = bucket.at(0).getpoint();
			for (int j = 1; j < bucket.size(); j++)
			{	//for every  point in centroid's cluster
				
				for (int z = 0; z < bucket.at(j).getpoint().size(); z++)	//compute sum of all the vectors of the Pointd
				{
						
					kmedoid.at(z) += bucket.at(j).getpoint().at(z);
				}
			}
			
			
			for (int z = 0; z < kmedoid.size(); z++)//divide vector of sum of all the vectors of the Points by the number of the points in the cluster
			{
				kmedoid.at(z) /= bucket.size();
			}

			new_centroids.push_back(Point(kmedoid));
		}
		cout << endl;

		/*for(int i = 0; i < new_centroids.size(); i++){ //print previous vs current centroids
			cout << "old_centroid\n";
			centroids.at(i).print_point();
			cout << "new_centroid\n";
			new_centroids.at(i).print_point();
		}*/
		
		/*if(count % 2 == 0 ){		//alternative stop like max iterator

			for(int i = 0; i < centroids.size(); i++){
				second_to_last_centroids.at(i) = centroids.at(i);
			}
		}*/
	

		if(configuration(centroids,new_centroids) == 1 || configuration(new_centroids, second_to_last_centroids) == 1 || count > max_iter){
			// optimal centroids found, no change in new c entroids
			config_call = 1;
			
		} 
		else{ // update centroids

			for(int i = 0; i < new_centroids.size(); i++){
				centroids.at(i) = new_centroids.at(i);
			}
			ClusterLsh.clear();
		}
		count++;
	}

	//assignment

	int min_dist = manhattanDistance(centroids.at(0).getpoint(), centroids.at(1).getpoint());


	for(int i = 0; i < centroid_size; i++){ // finds the nearest centroids ie the once with minimum manhattan Distance


		for(int j = 0; j < centroid_size; j++){

			if(i != j){ // avoid computing distance of a cetroid with itself, since manhattan distance in this case = 0 
				distance = manhattanDistance(centroids.at(i).getpoint(), centroids.at(j).getpoint());
				if(min_dist > distance && distance > 0)
					min_dist = distance;
			}else
				continue;
		}

	}

	radius = min_dist/2;

	//ClusterLsh.clear();
	int previous_size = 0;



	while(1){

	
		for (int j = 0; j < lsh.getl(); j++) // for each hash table in LSH 
		{
			for (int i = 0; i < centroids.size(); i++) // for each centroid
			{
				// hash the centroid to eash LSH hash table and find the bucket it would have normally enter if it was a non centroid point
				int position = (lsh.getgi().at(j).gfunction(centroids.at(i)))%table_size; 

				for (int w = 0; w < lsh.hash_tables[j][position].size(); w++) // for each point in that bucket
				{	
					int ClosestCentroidFlag = 0; // flag = 0 has not been assigned a centroid, flag = 1 has already been visited
				 	    						 // but can be also revisited by another centroid in order to be assigned at that centroid
												 // in case the later centroid is the closest to the point 
					mindistance = manhattanDistance(centroids.at(i).getpoint(), lsh.hash_tables[j][position].at(w).getpoint());
					
					int min_pos = 0;	
					for (int t = 0; t < centroids.size(); t++) // for each centroid find the one closest to the point
					{
						if(t != i){ //don't compare each centroid with itself
							int positionT = (lsh.getgi().at(j).gfunction(centroids.at(t)))%table_size;
							if(position == positionT){
								distance = manhattanDistance(centroids.at(t).getpoint(), lsh.hash_tables[j][positionT].at(w).getpoint());
								if (mindistance > distance){
									mindistance = distance;
									min_pos = t;
								}
							}	
						}
					}

					if(lsh.hash_tables[j][position].at(w).get_flag() == -1){
							lsh.hash_tables[j][position].at(w).set_flag(min_pos);
							
					}
					// now checking query balls and insert point to closest cluster / centroid
					if(mindistance <= radius && mindistance > radius/2){

						if(lsh.hash_tables[j][position].at(w).get_flag() != -1){

							unordered_map<int, vector<Point>>::iterator Rev_LSHfind;
								
							Rev_LSHfind = ClusterLsh.find(centroids.at(i).get_id());

							if (Rev_LSHfind ==	ClusterLsh.end()){ //first assignment of a point to the centroid
								
								vector<Point> bucket;
								bucket.push_back(lsh.hash_tables[j][position].at(w));
								ClusterLsh.insert({centroids.at(i).get_id(), bucket});
								}
							else {
								ClusterLsh[centroids.at(i).get_id()].push_back(lsh.hash_tables[j][position].at(w));								
							}		
						} 
					}
				}
			}		
		}

		int current_size = 0;
		for (int i = 0; i < centroids.size(); i++)
		{
			current_size += ClusterLsh[centroids.at(i).get_id()].size();
		}

		if(previous_size <= current_size - max_iter){	//max iterator

			previous_size = current_size;
			radius *= 2;

		}else{

			break;
		}
	}//end while 

	int position = 0;
	//ta stoixeia poy emeinan ektos ball-queries
	for (int i = 0; i < lsh.getl(); i++)
	{
		for (int j = 0; j < lsh.getk(); j++)
		{
			for (int z = 0; z < lsh.hash_tables[i][j].size(); z++)
			{
				if (lsh.hash_tables[i][j].at(z).get_flag() == -1)
    			{
    				position = 0;
    				int minimum = manhattanDistance(lsh.hash_tables[i][j].at(z).getpoint(), centroids.at(0).getpoint());
    				for (int c = 1; c < centroid_size; c++)
					{
						distance = manhattanDistance(lsh.hash_tables[i][j].at(z).getpoint(), centroids.at(c).getpoint());
						if (minimum > distance){
							minimum = distance;
							position = c;
						}			
					}
					ClusterLsh[centroids.at(position).get_id()].push_back(lsh.hash_tables[i][j].at(z));
					lsh.hash_tables[i][j].at(z).set_flag(position);
    			}	
			}
		}
	}

	//end clock
	auto t2 = chrono::steady_clock::now();
	auto TimeResult = t2 - t1;
	auto time_spanCluster = (chrono::duration<double, milli>(TimeResult).count()) / 1000;

	this->clusterTime += time_spanCluster;

	/*for(int i = 0; i < centroids.size(); i++){	//printing functiom

		vector<Point> v = ClusterLsh[centroids.at(i).get_id()]; 

		for(int j = 0; j < v.size(); j++){
			//v.at(j).getpoint();
			cout << v.at(j).get_id() << " ";
		}
		cout << endl;
	} */
}

Rev_Cube::Rev_Cube(vector<Point> dataset, int kappa, int m, int M,  int w, int d, int dtonos, int max_iter = 20, int tolerance = 10) : 
	Cluster(dataset, kappa), hypercube(dataset, m, M, w, d, dtonos)
{

	int centroid_size = centroids.size();
	int distance = 0;
	int mindistance;
	int radius = 0;
	int config_call = 0;
	
	vector<Point> second_to_last_centroids;
	for(int i = 0; i < centroids.size(); i++){
		second_to_last_centroids.push_back(centroids.at(i));
	}
	int count = 1;

	//clock start
	auto t1 = chrono::steady_clock::now();
	while(config_call == 0){
		for(int i = 0; i < centroids.size(); i++){ // for each centroid

			long int key = 0;

			for(int j = 0; j < dtonos; j++){ // construct key

				key = key*10 + hypercube.get_fresults().at(j).ffunction(centroids.at(i));
			}

			unordered_map<long int, vector<Point>>::iterator flag;

			flag = hypercube.cube.find(key);

			if (flag !=	hypercube.cube.end()){ // if centroid is  hashed at an existing key in the cube
					
				for(int w = 0; w < hypercube.cube[key].size(); w++){
					
					mindistance = manhattanDistance(centroids.at(i).getpoint(),  hypercube.cube[key].at(w).getpoint());
					int min_pos = 0;	
					for (int t = 0; t < centroids.size(); t++) // for each centroid find the one closest to the point
					{
						if(t != i){ //don't compare each centroid with itself

							long int keyT = 0;
							for(int s = 0; s < dtonos; s++){

								keyT = keyT*10 + hypercube.get_fresults().at(s).ffunction(centroids.at(t));

							}
							if(key == keyT){
								distance = manhattanDistance(centroids.at(t).getpoint(), hypercube.cube[keyT].at(w).getpoint());
								if (mindistance > distance){
									mindistance = distance;
									min_pos = t;
								}
							}	
						}
					}


					if(hypercube.cube[key].at(w).get_flag() == -1){
						
						hypercube.cube[key].at(w).set_flag(min_pos);
						
					}

					if(hypercube.cube[key].at(w).get_flag() != -1){

						unordered_map<int, vector<Point>>::iterator Rev_Cubefind;
								
						Rev_Cubefind = ClusterCube.find(centroids.at(i).get_id());

						if (Rev_Cubefind ==	ClusterCube.end()){ //first assignment of a point to the centroid
							
							
							vector<Point> bucket;
							bucket.push_back(hypercube.cube[key].at(w));
							ClusterCube.insert({centroids.at(i).get_id(), bucket});
						}
						else {
						
							ClusterCube[centroids.at(i).get_id()].push_back(hypercube.cube[key].at(w));								
						}
					} 
				}
			}
		}

		vector<Point> new_centroids;
		for (int i = 0; i < centroid_size; i++)
		{	//for every centroid 
				
			vector<int> kmedoid;
			vector<Point> bucket = ClusterCube[centroids.at(i).get_id()];
			if(bucket.size() == 0){  // initialize bucket if empty
				for (int z = 0; z < centroids.at(i).getpoint().size(); z++)
				{
					kmedoid.push_back(0);
				}
				new_centroids.push_back(Point(kmedoid));
				continue;
			}
			kmedoid = bucket.at(0).getpoint();
			for (int j = 1; j < bucket.size(); j++)
			{	//for every  point in centroid's cluster
				
				for (int z = 0; z < bucket.at(j).getpoint().size(); z++)	//compute sum of all the vectors of the Pointd
				{
						
					kmedoid.at(z) += bucket.at(j).getpoint().at(z);
				}
			}
			
			for (int z = 0; z < kmedoid.size(); z++)//divide vector of sum of all the vectors of the Points by the number of the points in the cluster
			{
				kmedoid.at(z) /= bucket.size();
			}

			new_centroids.push_back(Point(kmedoid));
		}
		cout << endl;

		/*for(int i = 0; i < new_centroids.size(); i++){ //update centroids
			cout << "old_centroid\n";
			//cout << "count " << count << endl;
			centroids.at(i).print_point();
			cout << "new_centroid\n";
			//cout << "count " << count << endl;
			new_centroids.at(i).print_point();

		}*/
		
		if(count % 2 == 0 ){

			for(int i = 0; i < centroids.size(); i++){
				second_to_last_centroids.at(i) = centroids.at(i);
			}
		}

		if(configuration(centroids,new_centroids) == 1 || configuration(new_centroids, second_to_last_centroids) == 1 || count > max_iter){
			config_call = 1;
			
		} 
		else{ // update centroids

			for(int i = 0; i < new_centroids.size(); i++){
				centroids.at(i) = new_centroids.at(i);
			}
			ClusterCube.clear();
		}
		count++;
	}

	int min_dist = manhattanDistance(centroids.at(0).getpoint(), centroids.at(1).getpoint());


	for(int i = 0; i < centroid_size; i++){ // finds the nearest centroids ie the once with minimum manhattan Distance


		for(int j = 0; j < centroid_size; j++){

			if(i != j){ // avoid computing distance of a cetroid with itself, since manhattan distance in this case = 0 
				distance = manhattanDistance(centroids.at(i).getpoint(), centroids.at(j).getpoint());
				//cout << "distance " << distance << endl;
				if(min_dist > distance && distance > 0)
					min_dist = distance;
			}
			else
				continue;
		}
	}

	//cout << "min_dist " << min_dist << endl;

	radius = min_dist/2;

	int previous_size = 0;

	while(1){
		for(int i = 0; i < centroids.size(); i++){ // for each centroid

			long int key = 0;

			for(int j = 0; j < dtonos; j++){ // construct key

				key = key*10 + hypercube.get_fresults().at(j).ffunction(centroids.at(i));
			}

			unordered_map<long int, vector<Point>>::iterator flag;

			flag = hypercube.cube.find(key);

			if (flag !=	hypercube.cube.end()){ // if centroid is  hashed at an existing key in the cube
			
			
				for(int w = 0; w < hypercube.cube[key].size(); w++){
					

					mindistance = manhattanDistance(centroids.at(i).getpoint(),  hypercube.cube[key].at(w).getpoint());
					int min_pos = 0;	
					for (int t = 0; t < centroids.size(); t++) // for each centroid find the one closest to the point
					{
						if(t != i){ //don't compare each centroid with itself

							long int keyT = 0;
							for(int s = 0; s < dtonos; s++){

								keyT = keyT*10 + hypercube.get_fresults().at(s).ffunction(centroids.at(t));
							}
							if(key == keyT){ // an vroume duo centroid einai stin idia korifi
								distance = manhattanDistance(centroids.at(t).getpoint(), hypercube.cube[keyT].at(w).getpoint());
								if (mindistance > distance){
									mindistance = distance;
									min_pos = t;
								}
							}	
						}
					}


					if(hypercube.cube[key].at(w).get_flag() == -1){
						
						hypercube.cube[key].at(w).set_flag(min_pos);
						
					}
					if(mindistance <= radius && mindistance > radius/2){

						if(hypercube.cube[key].at(w).get_flag() != -1){

							unordered_map<int, vector<Point>>::iterator Rev_Cubefind;
									
							Rev_Cubefind = ClusterCube.find(centroids.at(i).get_id());

							if (Rev_Cubefind ==	ClusterCube.end()){ //first assignment of a point to the centroid
								
								vector<Point> bucket;
								bucket.push_back(hypercube.cube[key].at(w));
								ClusterCube.insert({centroids.at(i).get_id(), bucket});
							}
							else {
							
								ClusterCube[centroids.at(i).get_id()].push_back(hypercube.cube[key].at(w));								
							}
						} 
					}
				}			
			}
		}

		int current_size = 0;
		for (int i = 0; i < centroids.size(); i++)
		{
			current_size += ClusterCube[centroids.at(i).get_id()].size();
		}

		//cout << "previous_size " << previous_size << " current_size " << current_size << endl;

		if(previous_size <= current_size - tolerance){	// terminantion function

			previous_size = current_size;
			radius *= 2;
		}
		else{
			break;
		}
	}

	// ta stoixeia poy emeinan ektos ball query ginetai manhattan distance me ta centroids
	// kai to kanoyme insert sto plisiestero apo ta centroids (lines 804 - 823)
	int minimum = 0;
	int position = 0;
	for (auto& x: hypercube.cube) {

    	for(int i = 0; i < x.second.size(); i++){
    		if (x.second.at(i).get_flag() == -1)
    		{
    			position = 0;
    			minimum = manhattanDistance(x.second.at(i).getpoint(), centroids.at(0).getpoint());
    			for (int j = 1; j < centroid_size; j++)
				{
					distance = manhattanDistance(x.second.at(i).getpoint(), centroids.at(j).getpoint());
					if (minimum > distance){
						minimum = distance;
						position = j;
					}			
				}
				ClusterCube[centroids.at(position).get_id()].push_back(x.second.at(i));
				x.second.at(i).set_flag(position);
    		}
    	}	
  	}	
	//end clock here!
  	auto t2 = chrono::steady_clock::now();
	auto TimeResult = t2 - t1;
	auto time_spanCluster = (chrono::duration<double, milli>(TimeResult).count()) / 1000;

	this->clusterTime += time_spanCluster;
}



vector<double> Cluster::Silhouette(unordered_map<int, vector<Point>> clusters){
	vector<double> ai;		// a(i)
	vector<double> bi;		// b(i)
	vector<double> si;		// s(i)
	double distance;
	double mindistance;
	int position = 0;    //the id for nearest cluster b(i)

	for (auto& x: clusters) {

    	Point q = x.second.at(0);
    	int sum = 0;
    	for(int i = 1; i < x.second.size(); i++){
    		sum += manhattanDistance(q.getpoint(), x.second.at(i).getpoint());	//athroisma apostaseon
    	}
    	ai.push_back(sum / x.second.size()-1); 	// (ai's avg)
    	// edo tha vriskoume ton kontinero cluster apo ton torino (lines 851 -863)
    	if (x.first == centroids.at(0).get_id())
    	{
    		mindistance = manhattanDistance(centroids.at(1).getpoint(), q.getpoint());
    		position = 1;
    		for (int i = 2; i < centroids.size(); i++)
    		{
    			distance = manhattanDistance(centroids.at(i).getpoint(), q.getpoint());
    			if(mindistance > distance){
    				mindistance = distance;
    				position = i;
    			}
    		}
    	}
    	else{
    		//same with if but not for the first centroid
			mindistance = manhattanDistance(centroids.at(0).getpoint(), q.getpoint());
			position = 0;
    		for (int i = 1; i < centroids.size(); i++)
    		{
    			if(x.first != centroids.at(i).get_id()){
    				distance = manhattanDistance(centroids.at(i).getpoint(), q.getpoint());
    				
    				if(mindistance > distance){
    					mindistance = distance;
    					position = i;
    				}
    			}
    		}  		
  		}
  		sum = 0;
  		// bi avg
  		for (int i = 0; i < clusters[centroids.at(position).get_id()].size(); i++)
  		{
  			sum += manhattanDistance(q.getpoint(), clusters[centroids.at(position).get_id()].at(i).getpoint());
  		}
  		bi.push_back(sum / clusters[centroids.at(position).get_id()].size()-1);
    }
    // silhouette mathematical fraction type here
    for (int i = 0; i < ai.size(); i++)
    {
    	if (ai.at(i) > bi.at(i))
    	{
    		si.push_back(double(bi.at(i) - ai.at(i))/ai.at(i));
    	}
    	else{
    		si.push_back(double(bi.at(i) - ai.at(i))/bi.at(i));
    	}
    }
    double totalsi = 0.0;
    for (int i = 0; i < si.size(); i++)
    {
    	totalsi += si.at(i);
    }
    totalsi /= si.size();
    si.push_back(totalsi);	// store fractions here

    return si;
}

void Cluster::query_cluster(char* output_file, unordered_map<int, vector<Point>> clusters, char* method, int complete){
	ofstream ofile;

	ofile.open(output_file, ios::out);
	if (ofile.is_open()){cout << "is open\n";}
	else {cout << "is closed\n";}
	ofile << "Algorithm: " << method << endl;
	int count = 1;
	for (auto& x: clusters) {
		ofile << "CLUSTER-" << count++ << " ";
		ofile << "{size: " << x.second.size() << ", centroid: "; 
    	for(int i = 0; i < centroids.size(); i++){
  			if (x.first == centroids.at(i).get_id())
  			{
  				for(int j=0; j< centroids.size(); j++)
					ofile << centroids.at(i).getpoint().at(j) << " ";
				ofile << "}" << endl;
  			}	
    	}	
  	}	
  	ofile << "clustering_time: " << getClusterTime() << endl;
  	vector <double> SilhouetteResults = Silhouette(clusters);
  	ofile << "Silhouette: [ ";
  	for (int i = 0; i < SilhouetteResults.size(); i++)
	{
		if (i == SilhouetteResults.size() - 1){
			ofile << SilhouetteResults.at(i) << " " ;
		}
		else{
			ofile << SilhouetteResults.at(i) << ", " ;
		}
	}
	ofile << "]" << endl;
	if (complete == 1)
	{
		count = 1;
		for (auto& x: clusters){
			ofile << "CLUSTER-" << count++ << " ";
			ofile << "{ " << x.first << ", "; 
	    	for(int i = 0; i < x.second.size(); i++){
	    		if (i == x.second.size() - 1){
					ofile << x.second.at(i).get_id() << " ";
				}
				else{
					ofile << x.second.at(i).get_id() << ", ";
				}	
	    	}
	    	ofile << "}" << endl;	
  		}
	}
}
/*
Algorithm: Lloyds OR Range Search LSH OR Range Search Hypercube
CLUSTER-1 {size: <int>, centroid: πίνακας με τις συντεταγμένες του centroid}
. . . . . . .
CLUSTER-Κ {size: <int>, centroid: πίνακας με τις συντεταγμένες του centroid}
clustering_time: <double> //in seconds
Silhouette: [s1,...,si,...,sΚ, stotal]
/* si=average s(p) of points in cluster i, stotal=average s(p) of points in
dataset */
/* Optionally with command line parameter –complete 
CLUSTER-1 {centroid, image_numberA, ..., image_numberX}
. . . . . . .
CLUSTER-Κ {centroid, image_numberR, ..., image_numberZ}
*/