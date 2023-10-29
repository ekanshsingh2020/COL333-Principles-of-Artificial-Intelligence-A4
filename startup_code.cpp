#include <bits/stdc++.h>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

#define N 4

unordered_map<string, int> idx;
// unordered_map<string, int> n_values;
vector<vector<int>> data;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

public:

	string Node_Name;		// Variable name
	vector<int> Children;	// Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;			// Number of categories a variable represented by this node can take
	vector<string> values;	// Categories of possible values
	vector<float> CPT;		// conditional probability table as a 1-d array . Look for BIF format to understand its meaning
	vector<int> num_samples; // number of samples for each combination of values of parents
	vector<int> total_samples;

	vector<int> num_samples_copy; // number of samples for each combination of values of parents
	vector<int> total_samples_copy;

// public:
	// Constructor- a node is initialised with its name and its categories
	Graph_Node(string name, int n, vector<string> vals){

		Node_Name = name;
		nvalues = n;
		values = vals;

		num_samples.clear();
		total_samples.clear();

		num_samples_copy.clear();
		total_samples_copy.clear();

		// (nvalues+1) * 2^(number of parents * 2) {2 bits for each parent since max number of values can be 4} 
		// CPT.resize((nvalues+1) * (1 << ((int)Parents.size() * 2)));
		num_samples.resize((nvalues) * (1 << ((int)Parents.size() * 2)));
		total_samples.resize(1 << ((int)Parents.size() * 2));

		num_samples_copy.resize((nvalues) * (1 << ((int)Parents.size() * 2)));
		total_samples_copy.resize(1 << ((int)Parents.size() * 2));

		return;
	}

	string get_name(){
		return Node_Name;
	}

	vector<int> get_children(){
		return Children;
	}

	vector<string> get_Parents(){
		return Parents;
	}

	vector<float> get_CPT(){
		return CPT;
	}

	int get_nvalues(){
		return nvalues;
	}

	vector<string> get_values(){
		return values;
	}

	void set_CPT(vector<float> new_CPT){
		CPT.clear();
		CPT = new_CPT;
	}
	
	void set_Parents(vector<string> Parent_Nodes){
		Parents.clear();
		Parents = Parent_Nodes;
	}

	// add another node in a graph as a child of this node
	int add_child(int new_child_index){
		for (int i = 0; i < Children.size(); i++){
			if (Children[i] == new_child_index)
				return 0;
		}
		Children.push_back(new_child_index);
		return 1;
	}
};

// The whole network represted as a list of nodes
class Network{

public:

	list<Graph_Node> Pres_Graph;

	int addNode(Graph_Node node){

		Pres_Graph.push_back(node);
		return 0;
	}

	int netSize(){
		return Pres_Graph.size();
	}

	// get the index of node with a given name
	int get_index(string val_name){
		
		list<Graph_Node>::iterator listIt;
		int count = 0;
		for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
		{
			if (listIt->get_name().compare(val_name) == 0)
				return count;
			count++;
		}
		return -1;
	}
	
	// get the node at nth index
	list<Graph_Node>::iterator get_nth_node(int n){
		list<Graph_Node>::iterator listIt;
		int count = 0;
		for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++){
			if (count == n){
				cout << "returning\n";
				// cout << listIt->Node_Name << endl;
				string str = listIt->Node_Name;
				return listIt;
			}
			count++;
		}

		cout << "node not found\n";
		assert(false);
		return listIt;
	}

	// get the iterator of a node with a given name
	list<Graph_Node>::iterator search_node(string val_name){
		list<Graph_Node>::iterator listIt;
		for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
		{
			if (listIt->get_name().compare(val_name) == 0)
				return listIt;
		}

		cout << "node not found\n";
		return listIt;
	}

};

Network read_network(){

	Network Alarm;
	string line;
	int find = 0;
	ifstream myfile("alarm.bif");
	string temp;
	string name;
	vector<string> values;

	if (myfile.is_open()){
		while (!myfile.eof()){
			stringstream ss;
			getline(myfile, line);
			ss.str(line);
			ss >> temp;
			if (temp.compare("variable") == 0){
				ss >> name;
				getline(myfile, line);

				stringstream ss2;
				ss2.str(line);
				for (int i = 0; i < 4; i++){
					ss2 >> temp;
				}
				values.clear();
				while (temp.compare("};") != 0){
					values.push_back(temp);
					ss2 >> temp;
				}

				cout << name << endl;
				Graph_Node new_node(name, values.size(), values);
				int pos = Alarm.addNode(new_node);
			}

			else if (temp.compare("probability") == 0){

				ss >> temp;
				ss >> temp;

				list<Graph_Node>::iterator listIt;
				list<Graph_Node>::iterator listIt1;
				listIt = Alarm.search_node(temp);
				int index = Alarm.get_index(temp);
				ss >> temp;
				values.clear();
				while (temp.compare(")") != 0){
					listIt1 = Alarm.search_node(temp);
					listIt1->add_child(index);
					values.push_back(temp);
					ss >> temp;
				}
				listIt->set_Parents(values);
				getline(myfile, line);
				stringstream ss2;

				ss2.str(line);
				ss2 >> temp;

				ss2 >> temp;

				vector<float> curr_CPT;
				string::size_type sz;

				while (temp.compare(";") != 0){
					curr_CPT.push_back(atof(temp.c_str()));
					ss2 >> temp;
				}

				listIt->set_CPT(curr_CPT);
			}
			else
			{
			}
		}
		if (find)
			myfile.close();
	}
	return Alarm;
}

// We will write EM algorithm.

void init(Network &Alarm){

	// for (auto it : Alarm.Pres_Graph)
	// 	n_values[it.get_name()] = it.get_nvalues();

	int num_nodes = Alarm.netSize();

	// data is nodes x samples
	data.resize(num_nodes);

	// initialize idx map
	int i = 0;
	for (auto it : Alarm.Pres_Graph){
		idx[it.Node_Name] = i;
		i++;
	}

	// open file for reading
	ifstream myfile("records.dat");
	string line;
	string temp;
	if (myfile.is_open()){
		while (!myfile.eof()){ // TODO: Check this!!!
			getline(myfile, line);
			// if(myfile.eof()) break;
			stringstream ss;
			ss.str(line);
			auto it = Alarm.Pres_Graph.begin();
			bool question_mark = false;
			for (int j = 0; j < num_nodes; j++){
				assert(it != Alarm.Pres_Graph.end());
				ss >> temp;
				// remove first and last character
				// temp = temp.substr(1, temp.size() - 2);
				if(temp.compare("\"?\"") == 0){
					question_mark = true;
					data[j].push_back(-1);
				}
				else{
					// check the index of this value in the list of values of this node
					int index = find(it->values.begin(), it->values.end(), temp) - it->values.begin();
					data[j].push_back(index);
				}
				it++;
			}
			if(!question_mark){
				//this means that this row is completely known, so use this to initiliaze the parameters
				auto it = Alarm.Pres_Graph.begin();
				for (int j = 0; j < num_nodes; j++){
					assert(it != Alarm.Pres_Graph.end());
					// update the num_samples table for this node
					int index = data[j].back();
					// get data of all its parents
					int val = 0;
					int i = ((int)it->get_Parents().size() - 1) * 2;
					for (auto parent : it->get_Parents()){
						int parent_idx = idx[parent];
						int par_val = data[parent_idx].back();
						int lsb = par_val & 1, msb = par_val & 2;
						val += (lsb << i) + (msb << (i+1));
						i -= 2;
					}

					i = ((int)it->get_Parents().size()) * 2;
					
					it->total_samples[val]++;
					
					int lsb = index & 1, msb = index & 2;
					val += (lsb << i) + (msb << (i+1));
					it->num_samples[val]++;
					it++;
				}
			}
		}
		myfile.close();
	}
	else
		cout << "Unable to open file\n";

	return;
}

void E_step(Network &Alarm){

	// need to calculate the probability of each node given its parents
	// for each node
	int num_nodes = Alarm.netSize();
	list<Graph_Node>::iterator it = Alarm.Pres_Graph.begin();
	
	vector<vector<int>> data_copy = data;


	for(int i = 0 ; i < num_nodes ; i++){

		assert(it != Alarm.Pres_Graph.end());
		// cout << i << endl;
		// copy num_samples and total_samples to num_samples_copy and total_samples_copy
		for(int j = 0 ; j < it->num_samples.size() ; j ++){
			it->num_samples_copy[j] = it->num_samples[j];
		}
		for(int j = 0 ; j < it->total_samples.size() ; j ++){
			it->total_samples_copy[j] = it->total_samples[j];
		}

		// for each value of the node
		for(int j = 0 ; j < data_copy[i].size() ; j ++){
			if(data_copy[i][j] == -1){

				cout << i << " " << j << " " << it->Node_Name << endl;

				// I need to estimate this probability given its parents in this column
				// get data of all its parents
				float prob_x_given_parents = 1.0;
				int val = 0;
				int ix = ((int)it->Parents.size() - 1) * 2;

				for(string& par : it -> Parents){
					int par_idx = idx[par];
					int par_val = data_copy[par_idx][j];
					int lsb = par_val & 1, msb = par_val & 2;
					val += (lsb << ix) + (msb << (ix+1));
					ix -= 2;
				}


				cout << i << " " << j << " " << it->Node_Name << endl;

				(it->total_samples[val])++;

				cout << i << " " << j << " " << it->Node_Name << endl;

				// with laplace smoothing
				int t_samples = it -> total_samples[val];
				t_samples += it -> nvalues;

				cout << i << " " << j << " " << it->Node_Name << endl;

				float prob_x[it->nvalues];

				ix = ((int)it->Parents.size()) * 2;

				cout << i << " " << j << " " << it->Node_Name << endl;

				for(int k = 0 ; k < it->nvalues ; k ++){
					int lsb = k & 1, msb = k & 2;
					int val1 = val + (lsb << ix) + (msb << (ix+1));
					int n_samples = it -> num_samples[val1];
					n_samples += 1;
					prob_x[k] = (float)n_samples / (float)t_samples;

					cout << i << " " << j << " " << it->Node_Name << endl;

					// markov blanket
					for(int p = 0 ; p < (int)it->Children.size() ; p ++){

						cout << i << " " << j << " " << it->Node_Name << endl;

						int child = it->Children[p];
						assert(child >= 0);
						int child_idx = child;
						int child_val = data_copy[child_idx][j];
						int lsb = child_val & 1, msb = child_val & 2;
						// calculate parents val
						int val = 0;
						cout << i << " " << j << " " << it->Node_Name << endl;
						cout << child << endl;
						list<Graph_Node>::iterator it_child = Alarm.get_nth_node(child);
						cout << i << " " << j << " " << it->Node_Name << endl;
						cout << it_child->num_samples.size() << endl;
						// assert(it_child != nullptr);
						cout << it_child->Node_Name << endl;
						int iy = ((int)it_child->Parents.size() - 1) * 2;
						for (string &parent : it_child->Parents){
							int parent_idx = idx[parent];
							int par_val = data_copy[parent_idx][j];
							int lsb = par_val & 1, msb = par_val & 2;
							val += (lsb << iy) + (msb << (iy+1));
							iy -= 2;
						}
						cout << i << " " << j << " " << it->Node_Name << endl;

						int t_samples_child = it_child -> total_samples[val];
						t_samples_child += it_child -> nvalues;

						iy = ((int)it_child->Parents.size()) * 2;
						
						val += (lsb << iy) + (msb << (iy+1));
						int n_samples_child = it_child -> num_samples[val];
						n_samples_child += 1;

						prob_x[k] *= (float)n_samples_child / (float)t_samples_child;
						cout << i << " " << j << " " << it->Node_Name << endl;
					}
				}

				// now sample from this distribution
				float r = (float)rand() / (float)RAND_MAX;
				float sum = 0.0;
				int k = 0;
				for(k = 0 ; k < it->nvalues ; k ++){
					sum += prob_x[k];
					if(sum > r)
						break;
				}
				data_copy[i][j] = k;
				int lsb = k & 1, msb = k & 2;
				val += (lsb << i) + (msb << (i+1));
				(it->num_samples[val])++;
				cout << i << " " << j << " " << it->Node_Name << endl;
			}
		}
		it++;
	}

	return;
}

void M_step(Network &Alarm){

	// write to the solved_alarm.bif file
	ofstream myfile("solved_alarm.bif");
	if (myfile.is_open()){
		// open alarm.bif file
		// ifstream myfile2("alarm.bif");
		// string line;
		// string temp;
		// int find = 0;
		// if (myfile2.is_open()){
		// 	while (!myfile2.eof()){
		// 		getline(myfile2, line);
		// 		stringstream ss;
		// 		ss.str(line);
		// 		ss >> temp;
		// 		if (temp.compare("variable") == 0){
		// 			// just paste these 4 lines into myfile
		// 			myfile << line << "\n";
		// 			getline(myfile2, line);
		// 			myfile << line << "\n";
		// 			getline(myfile2, line);
		// 			myfile << line << "\n";
		// 			getline(myfile2, line);
		// 			myfile << line << "\n";
		// 		}
		// 		else if (temp.compare("probability") == 0){
		// 			// write this line to output file
		// 			myfile << line << "\n";
		// 			// get the name of the node
		// 			ss >> temp;
		// 			ss >> temp;
		// 			// get the node
		// 			// remove first and last character
		// 			temp = temp.substr(1, temp.size() - 2);
		// 			// get the node
		// 			list<Graph_Node>::iterator listIt;
		// 			listIt = Alarm.search_node(temp); //TODO: Should be linear
					
		// write all probabilities
		for(auto node : Alarm.Pres_Graph){
			myfile << node.get_name() << endl;
			int num_parents = node.get_Parents().size();
			int num_values = node.get_nvalues();
			int num_combinations = 1 << (2 * num_parents);
			for(int i = 0 ; i < num_combinations ; i ++){
				int val = 0;
				int j = (num_parents - 1) * 2;
				for(auto par : node.get_Parents()){
					int par_idx = idx[par];
					int par_val = i >> j;
					int lsb = par_val & 1, msb = par_val & 2;
					val += (lsb << par_idx) + (msb << (par_idx + 1));
					j -= 2;
				}
				int n_samples = node.num_samples_copy[val];
				int t_samples = node.total_samples_copy[i];
				n_samples += 1;
				t_samples += node.nvalues;
				float prob = (float)n_samples / (float)t_samples;
				myfile << prob << " ";
			}
			myfile << ";" << endl;
		}
		myfile.close();
	}
	else
		cout << "Unable to open file\n";

	return;

}

int main(){

	// start_time
	auto start_time = chrono::high_resolution_clock::now();
	
	Network Alarm;
	Alarm = read_network();

	int max_num_parents = 0;
	int max_nvalues = 0;

	for (auto it : Alarm.Pres_Graph){
		max_num_parents = max(max_num_parents, (int)it.get_Parents().size());
		max_nvalues = max(max_nvalues, it.get_nvalues());
	}

	// Initialize all the parameters
	init(Alarm);

	// do while current time - start_time < 2 mins
	while (1){

		// E step
		E_step(Alarm);

		// current_time
		auto current_time = chrono::high_resolution_clock::now();

		// if current_time - start_time > 0.5 mins
		if (chrono::duration_cast<chrono::minutes>(current_time - start_time).count() > 0.5)
			break;

		// M step
		M_step(Alarm);

		// current_time
		current_time = chrono::high_resolution_clock::now();

		// if current_time - start_time > 2 mins
		if (chrono::duration_cast<chrono::minutes>(current_time - start_time).count() > 0.5)
			break;

	}

	return 0;
}