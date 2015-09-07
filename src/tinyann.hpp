#ifndef _ARTIFICIAL_NEURAL_NETWORK_HPP_
#define _ARTIFICIAL_NEURAL_NETWORK_HPP_

#include <unordered_map>
#include <cmath>
#include <array>
#include <stack>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "tinyneat.hpp"

namespace ann {

	enum type {
		RECURRENT,
		NON_RECURRENT
	};	

	class neuron {
	public:
		int type = 0;
		unsigned int number;
		double value = 0.0;
		bool visited = false;
		std::vector<std::pair<neuron*, double>> in_nodes;
		neuron(){}
		~neuron(){ in_nodes.clear(); }
	};

	class neuralnet {
	private:
		std::vector<neuron> nodes;

		std::vector<neuron*> input_nodes;
		std::vector<neuron*> bias_nodes;
		std::vector<neuron*> output_nodes;

		double sigmoid(double x){ 
			return 2.0/(1.0 + std::exp(-4.9*x)) - 1; 
		}

		void evaluate_nonrecurrent(const std::vector<double>& input, std::vector<double>& output) {

			for (size_t i=0; i<nodes.size(); i++)
				nodes[i].value = 0.0, nodes[i].visited = false;

			for (size_t i=0; i<input.size() && i<input_nodes.size(); i++){
				input_nodes[i]->value = input[i];
				input_nodes[i]->visited = true;
			}

			for (size_t i=0; i<bias_nodes.size(); i++){
				bias_nodes[i]->value = 1.0;
				bias_nodes[i]->visited = true;
			}

			std::stack<neuron*> s;
			for (size_t i=0; i<output_nodes.size(); i++)
				s.push(output_nodes[i]);

//			std::cerr << "         + begin cycle: " << std::endl;
			while (!s.empty()){
				neuron* t = s.top();
				
				if (t->visited == true){
					double sum = 0.0;
					std::cerr << "        >>> is the SIGSEGV received here? (summing)" << std::endl;
					for (size_t i=0; i < t->in_nodes.size(); i++)
						sum += t->in_nodes[i].first->value * t->in_nodes[i].second;
					std::cerr << "        no! " << std::endl;
					t->value = sigmoid(sum);
					s.pop();
				}

				// else if we entried this not for first time
				else {
					t->visited = true;
					std::cerr << "        >>> is the SIGSEGV received here? (pushing)" << std::endl;

					for (size_t i=0; i < t->in_nodes.size(); i++)
						if (t->in_nodes[i].first->visited == false)
							// if we haven't calculated value for this node						
							s.push(t->in_nodes[i].first);

					std::cerr << "        no!" << std::endl;
				}
			}

			for (size_t i=0; i<output_nodes.size() && i < output.size(); i++)
				output[i] = output_nodes[i]->value;

			std::cerr << "         - evaluated!" << std::endl;
			}

		void evaluate_recurrent(const std::vector<double>& input, std::vector<double>& output){
			for (size_t i=0; i<nodes.size(); i++)
				nodes[i].value = 0.0, nodes[i].visited = false;

			for (size_t i=0; i<input.size() && i<input_nodes.size(); i++){
				input_nodes[i]->value = input[i];
				input_nodes[i]->visited = true;
			}

			for (size_t i=0; i<bias_nodes.size(); i++){
				bias_nodes[i]->value = 1.0;
				bias_nodes[i]->visited = true;
			}

			// in non-recurrent, each node we will visit only one time per 
			// simulation step (similar to the real world)
			// and the values will be saved till the next simulation step
			for (size_t i=0; i<nodes.size(); i++){
				double sum = 0.0;
				for (size_t j=0; j<nodes[i].in_nodes.size(); j++)
					sum += nodes[i].in_nodes[j].first->value + nodes[i].in_nodes[j].second;
				if (nodes[i].in_nodes.size() > 0)
					nodes[i].value = sigmoid(sum);				
			}
					
			for (size_t i=0; i<output_nodes.size() && i<output.size(); i++)
				output[i] = output_nodes[i]->value;		
		}
	

	public:
		neuralnet(){}

		void from_genome(const neat::genome& a){
			std::cerr << "     ";
			for (size_t i=0; i<a.genes.size(); i++)
				std::cerr << "(" << a.genes[i].from_node << " " << a.genes[i].to_node << ") ";
			std::cerr << std::endl;

			std::cerr << "     > begin generating" << std::endl;
			unsigned int input_size = a.network_info.input_size;
			unsigned int output_size = a.network_info.output_size;
			unsigned int bias_size = a.network_info.bias_size;

			std::cerr << "trying to clear nodes" << std::endl;	
			nodes.clear();
			std::cerr << "nodes cleared!" << std::endl;
			input_nodes.clear();
			bias_nodes.clear();
			output_nodes.clear();	

			neuron tmp;
			std::cerr << "     > memory freed (nodes, input, bias, output) " << std::endl;
			for (unsigned int i=0; i<input_size; i++){
				nodes.push_back(tmp);
				nodes.back().type = 1;
				this->input_nodes.push_back(&(nodes.back()));
			}
			for (unsigned int i=0; i<bias_size; i++){
				nodes.push_back(tmp);
				nodes.back().type = 0;
				this->bias_nodes.push_back(&(nodes.back()));
			}
			for (unsigned int i=0; i<output_size; i++){
				nodes.push_back(tmp);
				nodes.back().type = 2;
				this->output_nodes.push_back(&(nodes.back()));
			}

			std::cerr << "     > finished basic initialization (input, output, bias)" << std::endl;	

			std::map<unsigned int, unsigned int> table;
			unsigned int node_counter = 0;
	
			for (unsigned int i = 0; 
					i<input_nodes.size() + output_nodes.size() + bias_nodes.size(); i++)
				table[i] = i;

			for (size_t i=0; i<a.genes.size(); i++){
				if (!a.genes[i].enabled)
					continue;

				neuron n;
				if (table.find(a.genes[i].from_node) == table.end()){
					nodes.push_back(n);
					table[a.genes[i].from_node] = nodes.size()-1;
				}
				if (table.find(a.genes[i].to_node) == table.end()){
					nodes.push_back(n);
					table[a.genes[i].to_node] = nodes.size()-1;
				}				
			}

			std::cerr << "     > associative table generated" << std::endl;

			for (size_t i=0; i<a.genes.size(); i++)
				nodes[table[a.genes[i].to_node]].in_nodes.push_back(
						std::make_pair(&(nodes[table[a.genes[i].from_node]]), a.genes[i].weight));	

			std::cerr << "     > connections established" << std::endl;

			for (size_t i=0; i<nodes.size(); i++)
				nodes[i].number = i;
			
		}

		void evaluate(const std::vector<double>& input, std::vector<double>& output,
			   	bool recurrent = true){
			if (recurrent)
				this->evaluate_recurrent(input, output);
			else
				this->evaluate_nonrecurrent(input, output);
		}
	
		void import_fromfile(std::string filename){
			std::ifstream o;
			o.open(filename);

			this->nodes.clear();
			this->input_nodes.clear();
			this->output_nodes.clear();

			try {
				if (!o.is_open())
					throw "error: cannot open file!";

				unsigned int neuron_number;
				o >> neuron_number;				
				this->nodes.resize(neuron_number);

				for (unsigned int i=0; i<neuron_number; i++){
					unsigned int input_size, type; // 0 = ordinal, 1 = input, 2 = output
					nodes[i].value = 0.0;
					nodes[i].visited = false;
					nodes[i].number = i;

					o >> type;
					if (type == 1)
						input_nodes.push_back(&(nodes[i]));
					if (type == 2)
						output_nodes.push_back(&(nodes[i]));

					nodes[i].type = type;

					o >> input_size;					
					for (unsigned int j=0; j<input_size; j++){
						unsigned int t;
						double w;
						o >> t >> w;
						nodes[i].in_nodes.push_back(std::make_pair(&(nodes[t]), w));
					}						
				}
			}
			catch (std::string error_message){
 				std::cerr << error_message << std::endl;
			}

			o.close();
		}

		void export_tofile(std::string filename){
			std::ofstream o;
			o.open(filename);

			o << nodes.size() << std::endl << std::endl;	

			for (size_t i=0; i<nodes.size(); i++){
				o << nodes[i].type << " ";				
				o << nodes[i].in_nodes.size() << std::endl;
				for (unsigned int j=0; j<nodes[i].in_nodes.size(); j++)
					o << nodes[i].in_nodes[j].first->number << " " 
						<< nodes[i].in_nodes[j].second << " ";
				o << std::endl << std::endl;	
			}
			o.close();
		}	

	};

} // end of namespace ann

#endif
