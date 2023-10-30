#include <bits/stdc++.h>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

#define N 4

unordered_map<string, int> idx;
unordered_map<string, int> n_values;
vector<vector<int>> dat;

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

		num_samples.assign((nvalues) * (1 << (4 * 2 + 1)), 0);
		total_samples.assign(1 << (4 * 2 + 1), 0);

		num_samples_copy.assign((nvalues) * (1 << (4 * 2 + 1)), 0);
		total_samples_copy.assign(1 << (4 * 2 + 1), 0);

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
				return listIt;
			}
			count++;
		}

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

	for (auto it : Alarm.Pres_Graph)
		n_values[it.Node_Name] = it.nvalues;

	int num_nodes = Alarm.netSize();

	// dat is nodes x samples
	dat.clear();
	dat.resize(num_nodes);

	// initialize idx map
	int i_node_name = 0;
	idx.clear();
	for (auto it : Alarm.Pres_Graph){
		idx[it.Node_Name] = i_node_name;
		i_node_name++;
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
			auto alarm_pres_graph_iterator = Alarm.Pres_Graph.begin();
			bool question_mark = false;
			for (int j_node = 0; j_node < num_nodes; j_node++){
				assert(alarm_pres_graph_iterator != Alarm.Pres_Graph.end());
				ss >> temp;
				// remove first and last character
				// temp = temp.substr(1, temp.size() - 2);
				if(temp.compare("\"?\"") == 0){
					question_mark = true;
					dat[j_node].push_back(-1);
				}
				else{
					// check the index of this value in the list of values of this node
					int index = find(alarm_pres_graph_iterator->values.begin(), alarm_pres_graph_iterator->values.end(), temp) - alarm_pres_graph_iterator->values.begin();
					dat[j_node].push_back(index);
				}
				alarm_pres_graph_iterator++;
			}
			if(!question_mark){
				//this means that this row is completely known, so use this to initiliaze the parameters
				auto alarm_pres_graph_iterator_2 = Alarm.Pres_Graph.begin();
				for (int j_node = 0; j_node < num_nodes; j_node++){
					assert(alarm_pres_graph_iterator_2 != Alarm.Pres_Graph.end());
					// update the num_samples table for this node
					int index = dat[j_node][(int)dat[j_node].size() - 1];
					// get dat of all its parents
					int val = 0;
					int ix = ((int)alarm_pres_graph_iterator_2->Parents.size() - 1) * 2;
					for (auto parent : alarm_pres_graph_iterator_2->Parents){
						int parent_idx = idx[parent];
						int par_val = dat[parent_idx][(int)dat[parent_idx].size() - 1];
						int lsb = (par_val & 1);
						int msb = ((par_val >> 1) & 1);
						val += (lsb << ix) + (msb << (ix+1));
						ix -= 2;
					}

					ix = ((int)alarm_pres_graph_iterator_2->Parents.size()) * 2;
					
					alarm_pres_graph_iterator_2->total_samples[val]++;
					
					int lsb = (index & 1); 
					int msb = ((index >> 1) & 1);
					val += (lsb << ix) + (msb << (ix+1));
					alarm_pres_graph_iterator_2->num_samples[val]++;
					alarm_pres_graph_iterator_2++;
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
	list<Graph_Node>::iterator alarm_pres_graph_it = Alarm.Pres_Graph.begin();
	
	vector<vector<int>> dat_copy;
	dat_copy.clear();

	// copy dat
	for(int i_node = 0 ; i_node < num_nodes ; i_node++){
		dat_copy.push_back(dat[i_node]);
	}

	for(int i_node = 0 ; i_node < num_nodes ; i_node++){

		assert(alarm_pres_graph_it != Alarm.Pres_Graph.end());

		// for(int j = 0 ; j < (int)alarm_pres_graph_it->num_samples.size() ; j ++){
		// 	alarm_pres_graph_it->num_samples_copy[j] = alarm_pres_graph_it->num_samples[j];
		// }
		// for(int j = 0 ; j < (int)alarm_pres_graph_it->total_samples.size() ; j ++){
		// 	alarm_pres_graph_it->total_samples_copy[j] = alarm_pres_graph_it->total_samples[j];
		// }

		// for each value of the node
		for(int j_dat = 0 ; j_dat < (int)dat_copy[i_node].size() ; j_dat ++){

			if(dat_copy[i_node][j_dat] == -1){

				// I need to estimate this probability given its parents in this column
				// get dat of all its parents
				int val = 0;
				int ix = ((int)(alarm_pres_graph_it->Parents.size()) - 1) * 2;
				for(int i_par = 0 ; i_par < (int)alarm_pres_graph_it->Parents.size() ; i_par ++){
					string par = alarm_pres_graph_it->Parents[i_par];
					int par_idx = idx[par];
					int par_val = dat_copy[par_idx][j_dat];
					int lsb = (par_val & 1);
					int msb = ((par_val >> 1) & 1);
					val += ((lsb << ix) + (msb << (ix+1)));
					ix -= 2;
				}

				// alarm_pres_graph_it->total_samples[val]++;
				// with laplace smoothing
				int t_samples = alarm_pres_graph_it -> total_samples[val];
				t_samples += alarm_pres_graph_it -> nvalues;

				float prob_x[alarm_pres_graph_it->nvalues];

				// initialize with 1.0
				for(int k = 0 ; k < alarm_pres_graph_it->nvalues ; k ++){
					prob_x[k] = 1.0;
				}

				ix = ((int)(alarm_pres_graph_it->Parents.size())) * 2;

				for(int k = 0 ; k < alarm_pres_graph_it->nvalues ; k ++){
					int lsb = (k & 1);
					int msb = ((k >> 1) & 1);
					int val1 = val + (lsb << ix) + (msb << (ix+1));
					int n_samples = alarm_pres_graph_it -> num_samples[val1];
					n_samples += 1;
					prob_x[k] = (float)n_samples / (float)t_samples;

					// markov blanket
					for(int mb = 0 ; mb < (int)alarm_pres_graph_it->Children.size() ; mb ++){
						int child = alarm_pres_graph_it->Children[mb];
						int child_val = dat_copy[child][j_dat];
						int lsb = (child_val & 1);
						int msb = ((child_val >> 1) & 1);
						// calculate parents val
						int val_child = 0;
						auto it_child = Alarm.get_nth_node(child);
						int iy = ((int)it_child->Parents.size() - 1) * 2;
						for (auto parent : it_child->Parents){
							int parent_idx = idx[parent];
							int par_val = dat_copy[parent_idx][j_dat];
							int lsb = (par_val & 1);
							int msb = ((par_val >> 1) & 1);
							val_child += (lsb << iy) + (msb << (iy+1));
							iy -= 2;
						}

						int t_samples_child = it_child -> total_samples[val_child];
						t_samples_child += it_child -> nvalues;

						iy = ((int)it_child->Parents.size()) * 2;
						
						val_child += (lsb << iy) + (msb << (iy+1));

						int n_samples_child = it_child -> num_samples[val_child];
						n_samples_child += 1;

						prob_x[k] *= (float)n_samples_child / (float)t_samples_child;
					}
				}

				// now sample from this distribution
				float r = (float)rand() / (float)RAND_MAX;
				float sum = 0.0;
				int k_val = 0;
				for(k_val = 0 ; k_val < alarm_pres_graph_it->nvalues ; k_val ++){
					sum += prob_x[k_val];
					if(sum > r)
						break;
				}
				dat_copy[i_node][j_dat] = k_val;
				int lsb = (k_val & 1);
				int msb = ((k_val >> 1) & 1);
				(alarm_pres_graph_it->total_samples[val])++;
				ix = ((int)(alarm_pres_graph_it->Parents.size())) * 2;
				val += (lsb << ix) + (msb << (ix+1));
				(alarm_pres_graph_it->num_samples[val])++;
			}
		}
		alarm_pres_graph_it++;
	}

	return;
}

void generate_combination(Network &Alarm, Graph_Node &node, vector<int> &ans, int par_idx, int curr){
	
	int num_combinations = 1 << (2 * 4 + 1);
	if(par_idx >= (int)node.Parents.size()) {
		ans.push_back(curr);
		return;
	}

	for(int i = 0 ; i < n_values[node.Parents[par_idx]] ; i ++){
		int lsb = (i & 1);
		int msb = ((i >> 1) & 1);
		int val = curr + (lsb << (2 * par_idx)) + (msb << (2 * par_idx + 1));
		generate_combination(Alarm, node, ans, par_idx + 1, val);
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
			myfile << node.Node_Name << endl;
			int num_parents = node.Parents.size();
			// int num_parents = 4;
			// int num_values = node.get_nvalues();

			vector<int> combs;
			combs.clear();
			generate_combination(Alarm, node, combs, 0, 0);

			for(int i = 0 ; i < (int)combs.size() ; i ++){
				int val = combs[i];
				float t_samples = node.total_samples[val];
				// t_samples += 0.35*node.nvalues;
				t_samples += 1.0;
				vector<float> need_to_normalize;
				for(int k = 0 ; k < node.nvalues ; k++){
					int lsb = (k & 1);
					int msb = ((k >> 1) & 1);
					int val1 = val + (lsb << (2 * num_parents)) + (msb << (2 * num_parents + 1));
					int n_samples = node.num_samples[val1];
					n_samples += 1*0.35;
					float prob = (float)n_samples / (float)t_samples;
					if(prob > 1.0){
						cout << n_samples << " " << t_samples << endl;
					}
					assert(prob >= 0.0 && prob <= 1.0);
					need_to_normalize.push_back(prob);
				}
				float sum = 0.0;
				for(int k = 0 ; k < node.nvalues ; k ++){
					sum += need_to_normalize[k];
				}
				for(int k = 0 ; k < node.nvalues ; k ++){
					need_to_normalize[k] /= sum;
				}
				for(int k = 0 ; k < node.nvalues ; k ++){
					myfile << need_to_normalize[k] << " ";
				}
			}

			// int num_combinations = 1 << (2 * 4 + 1);


			// for(int i = 0 ; i < num_combinations ; i ++){
			// 	int val = i;
			// 	int t_samples = node.total_samples[val];
			// 	t_samples += node.nvalues;
			// 	for(int k = 0 ; k < node.get_nvalues() ; k++){
			// 		int lsb = (k & 1);
			// 		int msb = ((k >> 1) & 1);
			// 		int val1 = val + (lsb << (2 * num_parents)) + (msb << (2 * num_parents + 1));
			// 		int n_samples = node.num_samples[val1];
			// 		n_samples += 1;
			// 		float prob = (float)n_samples / (float)t_samples;
			// 		myfile << prob << " ";
			// 	}
			// }
			myfile << ";" << endl;
		}
		myfile.close();
	}
	else
		cout << "Unable to open file\n";

	return;

}

int main(){

	// seed for rand()
	srand(time(NULL));

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
	init(Alarm); // works fine

	M_step(Alarm);

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