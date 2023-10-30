#include <bits/stdc++.h>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

#define N 4
#define smoothing_factor 1

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

// public:
	// Constructor- a node is initialised with its name and its categories
	Graph_Node(string name, int n, vector<string> vals){

		Node_Name = name;
		nvalues = n;
		values = vals;

		num_samples.clear();
		total_samples.clear();

		num_samples.assign((nvalues) * (1 << (4 * 2 + 1)), 1);
		total_samples.assign(1 << (4 * 2 + 1), nvalues);

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

vector<pair<int, int>> question_mark_indices;

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
		while (!myfile.eof()){
			getline(myfile, line);
			stringstream ss;
			ss.str(line);
			auto alarm_pres_graph_iterator = Alarm.Pres_Graph.begin();
			bool question_mark = false;
			for (int j_node = 0; j_node < num_nodes; j_node++){
				assert(alarm_pres_graph_iterator != Alarm.Pres_Graph.end());
				ss >> temp;
				if(temp.compare("\"?\"") == 0){
					// question_mark = true;
					// dat[j_node].push_back(-1);
					// put a random value
					int index = rand() % alarm_pres_graph_iterator->nvalues;
					dat[j_node].push_back(index);
					question_mark_indices.push_back(make_pair(j_node, (int)dat[j_node].size() - 1));
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

// vector<vector<int>> dat_copy;

void E_step(Network &Alarm){

	// need to calculate the probability of each node given its parents
	// for each node
	int num_nodes = Alarm.netSize();

	// randomly choose a question_mark_index
	// int rand_index = rand() % (int)question_mark_indices.size();
	for(int rand_index = 0 ; rand_index < (int)question_mark_indices.size() ; rand_index++){

		int i_node = question_mark_indices[rand_index].first;
		int j_dat = question_mark_indices[rand_index].second;

		int val = 0;
		auto it = Alarm.get_nth_node(i_node);

		int ix = ((int)(it->Parents.size()) - 1) * 2;
		for(int i_par = 0 ; i_par < (int)it->Parents.size() ; i_par ++){
			string par = it->Parents[i_par];
			int par_idx = idx[par];
			int par_val = dat[par_idx][j_dat];
			int lsb = (par_val & 1);
			int msb = ((par_val >> 1) & 1);
			val += ((lsb << ix) + (msb << (ix+1)));
			ix -= 2;
		}

		// with laplace smoothing
		float t_samples = it -> total_samples[val];
		t_samples += (float)it -> nvalues * smoothing_factor;

		float prob_x[it->nvalues];

		// initialize with 1.0
		for(int k = 0 ; k < it->nvalues ; k ++){
			prob_x[k] = 1.0;
		}

		ix = ((int)(it->Parents.size())) * 2;

		for(int k = 0 ; k < it->nvalues ; k ++){
			int lsb = (k & 1);
			int msb = ((k >> 1) & 1);
			int val1 = val + (lsb << ix) + (msb << (ix+1));
			float n_samples = it -> num_samples[val1];
			n_samples += (float)1*smoothing_factor;
			prob_x[k] = (float)n_samples / (float)t_samples;
			assert(prob_x[k] >= 0.0 && prob_x[k] <= 1.0);

			// // markov blanket
			// for(int mb = 0 ; mb < (int)it->Children.size() ; mb ++){
			// 	int child = it->Children[mb];
			// 	int child_val = dat[child][j_dat];
			// 	int lsb_child = (child_val & 1);
			// 	int msb_child = ((child_val >> 1) & 1);
			// 	// calculate parents val
			// 	int val_child = 0;
			// 	auto it_child = Alarm.get_nth_node(child);
			// 	int iy = ((int)it_child->Parents.size() - 1) * 2;
			// 	for (auto parent : it_child->Parents){
			// 		int parent_idx = idx[parent];
			// 		int par_val = dat[parent_idx][j_dat];
			// 		if(parent_idx == i_node){
			// 			par_val = k;
			// 		}
			// 		int lsb_par = (par_val & 1);
			// 		int msb_par = ((par_val >> 1) & 1);
			// 		val_child += (lsb_par << iy) + (msb_par << (iy+1));
			// 		iy -= 2;
			// 	}

			// 	float t_samples_child = it_child -> total_samples[val_child];
			// 	t_samples_child += (float)it_child -> nvalues * smoothing_factor;

			// 	iy = ((int)it_child->Parents.size()) * 2;
				
			// 	val_child += (lsb_child << iy) + (msb_child << (iy+1));

			// 	float n_samples_child = it_child -> num_samples[val_child];
			// 	n_samples_child += (float)1*smoothing_factor;

			// 	if(n_samples_child > t_samples_child){
			// 		cout << n_samples_child << " " << t_samples_child << endl;
			// 	}

			// 	prob_x[k] *= (float)n_samples_child / (float)t_samples_child;
			// 	assert((float)n_samples_child / (float)t_samples_child <= 1.0);
			// 	assert(prob_x[k] >= 0.0 && prob_x[k] <= 1.0);
			// }
			
		}

		// float maxm = -1.0;
		// vector<int> maxm_idx;
		// for(int k = 0 ; k < it->nvalues ; k ++){
		// 	if(prob_x[k] > maxm){
		// 		maxm = prob_x[k];
		// 		maxm_idx.clear();
		// 		maxm_idx.push_back(k);
		// 	}
		// 	else if(prob_x[k] == maxm){
		// 		maxm_idx.push_back(k);
		// 	}
		// }

		// int k_val = maxm_idx[rand() % (int)maxm_idx.size()];
		// dat[i_node][j_dat] = k_val;
		float r = (float)rand() / (float)RAND_MAX;
		float sum = 0.0;
		int k_val = 0;
		for(k_val = 0 ; k_val < it->nvalues ; k_val ++){
			sum += prob_x[k_val];
			if(sum > r)
				break;
		}

		dat[i_node][j_dat] = k_val;

		// update the num_samples table for this node
		int index = dat[i_node][j_dat];
		// get dat of all its parents
		val = 0;
		ix = ((int)it->Parents.size() - 1) * 2;
		for (auto parent : it->Parents){
			int parent_idx = idx[parent];
			int par_val = dat[parent_idx][j_dat];
			int lsb = (par_val & 1);
			int msb = ((par_val >> 1) & 1);
			val += (lsb << ix) + (msb << (ix+1));
			ix -= 2;
		}

		it->total_samples[val]++;
		ix = ((int)it->Parents.size()) * 2;
		int lsb = (index & 1);
		int msb = ((index >> 1) & 1);
		val += (lsb << ix) + (msb << (ix+1));
		it->num_samples[val]++;
	}

	return;
}

void generate_combination(Network &Alarm, Graph_Node &node, vector<int> &ans, int par_idx, int curr, int ix){
	
	if(par_idx >= (int)node.Parents.size()) {
		ans.push_back(curr);
		return;
	}

	for(int i = 0 ; i < n_values[node.Parents[par_idx]] ; i ++){
		int lsb = (i & 1);
		int msb = ((i >> 1) & 1);
		int val = curr + (lsb << ix) + (msb << (ix + 1));
		generate_combination(Alarm, node, ans, par_idx + 1, val, ix - 2);
	}

	return;
}

void write_to_file(Network &Alarm){

	ofstream myfile("solved_alarm.bif");
	if (myfile.is_open()){

		ifstream myfile_read("alarm.bif");
		string line;

		auto it = Alarm.Pres_Graph.begin();

		if (myfile_read.is_open()){
			while (!myfile_read.eof()){
				getline(myfile_read, line);
				stringstream ss;
				ss.str(line);
				string temp;
				// if line does not begin with probability, just write as it is
				ss >> temp;
				if (temp.compare("probability") != 0){
					myfile << line << endl;
					continue;
				}
				else{
					myfile << line << endl;
					// get the name of the node
					ss >> temp;
					ss >> temp;

					assert(it->Node_Name.compare(temp) == 0);

					// now write the table
					getline(myfile_read, line);
					getline(myfile_read, line);

					myfile << "\ttable ";
					vector<int> combs;
					combs.clear();
					int ix = ((int)it->Parents.size() - 1) * 2;
					generate_combination(Alarm, *it, combs, 0, 0, ix);

					ix = ((int)it->Parents.size()) * 2;

					vector<vector<float>> need_to_normalize(it->nvalues, vector<float>(combs.size(), 0.0));

					for(int k = 0 ; k < it->nvalues ; k ++){
						int lsb = (k & 1);
						int msb = ((k >> 1) & 1);
						for(int comb = 0 ; comb < (int)combs.size() ; comb ++){
							int val = combs[comb];
							int val1 = val + (lsb << ix) + (msb << (ix + 1));
							float n_samples = it->num_samples[val1];
							float t_samples = it->total_samples[val];
							n_samples += (float)1*smoothing_factor;
							t_samples += (float)it->nvalues*smoothing_factor;
							float prob = (float)n_samples / (float)t_samples;
							if(prob > 1.0){
								cout << "Error\n";
								cout << "n_samples = " << n_samples << endl;
								cout << "t_samples = " << t_samples << endl;
							}
							assert(prob >= 0.0 && prob <= 1.0);
							need_to_normalize[k][comb] = prob;
						}
					}

					// now normalize
					for(int comb = 0 ; comb < (int)combs.size() ; comb ++){
						float sum = 0.0;
						for(int k = 0 ; k < it->nvalues ; k ++){
							sum += need_to_normalize[k][comb];
						}
						for(int k = 0 ; k < it->nvalues ; k ++){
							need_to_normalize[k][comb] /= sum;
						}
					}

					// now write to file
					for(int k = 0 ; k < it->nvalues ; k ++){
						for(int comb = 0 ; comb < (int)combs.size() ; comb ++){
							myfile << need_to_normalize[k][comb] << " ";
						}
					}

					myfile << ";" << endl;
					myfile << "}";
					it++;
					if(it != Alarm.Pres_Graph.end())
						myfile << endl;
				}
			}
			myfile_read.close();
		}
		else{
			cout << "Unable to open file\n";
			return;
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
	init(Alarm);

	int time_in_seconds = 15;

	while (1){

		// E step
		E_step(Alarm);
		// current_time
		auto current_time = chrono::high_resolution_clock::now();

		if (chrono::duration_cast<chrono::seconds>(current_time - start_time).count() > (time_t)(time_in_seconds-10))
			break;

		// M step
		// M_step(Alarm);

		// // current_time
		// current_time = chrono::high_resolution_clock::now();

		// if (chrono::duration_cast<chrono::seconds>(current_time - start_time).count() > (time_t)(time_in_seconds-10))
		// 	break;

	}

	// write to file
	write_to_file(Alarm);

	return 0;
}