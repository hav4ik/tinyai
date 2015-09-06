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


	typedef struct __neuron {
	public:		
		int type = 0; // 0 - ordinal, 1 - input, 2 - output
		unsigned int number;
		double value = 0.0;
		bool visited = false;
		std::vector<std::pair<__neuron*, double>> in_nodes;
	} neuron;

	class neuralnet {
	private:
		std::vector<neuron> nodes;

		std::vector<neuron*> input_nodes;
		std::vector<neuron*> bias_nodes;
		std::vector<neuron*> output_nodes;

		double sigmoid(double x){ 
			return 2.0/(1.0 + std::exp(-4.9*x)) - 1; 
		}

		std::vector<double>
		evaluate_nonrecurrent(std::vector<double>& input) {

			std::vector<double> answer(output_nodes.size(), 0.0);
			for (size_t i=0; i<nodes.size(); i++)
				nodes[i].value = 0.0, nodes[i].visited = false;

			for (size_t i=0; i<input.size(); i++)
				input_nodes[i]->value = input[i],
				input_nodes[i]->visited = true;			

			std::stack<neuron*> s;
			for (size_t i=0; i<output_nodes.size(); i++)
				s.push(output_nodes[i]);

			while (!s.empty()){
				neuron* t = s.top();
				
				if (t->visited == true){
					double sum = 0.0;
					for (size_t i=0; i < t->in_nodes.size(); i++)
						sum += t->in_nodes[i].first->value * t->in_nodes[i].second;
					t->value = sigmoid(sum);
					s.pop();
				}

				// else if we entried this not for first time
				else {
					t->visited = true;
					for (size_t i=0; i < t->in_nodes.size(); i++)
						if (t->in_nodes[i].first->visited == false)
							// if we haven't calculated value for this node
							s.push(t->in_nodes[i].first);
				}
			}

			for (size_t i=0; i<output_nodes.size(); i++)
				answer[i] = output_nodes[i]->value;

			return answer;
		}

		std::vector<double> 
		evaluate_recurrent(std::vector<double>& input){
			std::vector<double> answer(output_nodes.size(), 0.0);
			for (size_t i=0; i<nodes.size(); i++)
				nodes[i].value = 0.0, nodes[i].visited = false;

			for (size_t i=0; i<input.size(); i++)
				input_nodes[i]->value = input[i],
				input_nodes[i]->visited = true;			

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

			for (size_t i=0; i<output_nodes.size(); i++)
				answer[i] = output_nodes[i]->value;

			return answer;
		}
	

	public:
		neuralnet(){}

		void from_genome(const neat::genome& a){
			unsigned int input_size = a.network_info.input_size;
			unsigned int output_size = a.network_info.output_size;
			unsigned int bias_size = a.network_info.bias_size;
		
			nodes.resize(input_size + bias_size + output_size);

			for (unsigned int i=0; i<input_size; i++){
				this->input_nodes.push_back(&(nodes[i]));
				nodes[i].type = 1;				
			}
			for (unsigned int i=0; i<bias_size; i++){
				this->bias_nodes.push_back(&(nodes[input_size + i]));
				nodes[i].type = 0;
			}
			for (unsigned int i=0; i<output_size; i++){
				this->output_nodes.push_back(&(nodes[input_size + bias_size + i]));
				nodes[i].type = 2;
			}


			std::unordered_map<unsigned int, unsigned int> table;
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
					table[a.genes[i].from_node] = node_counter++;
				}
				if (table.find(a.genes[i].to_node) == table.end()){
					nodes.push_back(n);
					table[a.genes[i].to_node] = node_counter++;
				}				
			}

			for (size_t i=0; i<a.genes.size(); i++)
				nodes[a.genes[i].to_node].in_nodes.push_back(
						std::make_pair(&(nodes[a.genes[i].from_node]), a.genes[i].weight));

			for (size_t i=0; i<nodes.size(); i++)
				nodes[i].number = i;
			
		}

		std::vector<double> evaluate(std::vector<double>& input, bool recurrent){
			if (recurrent)
				return this->evaluate_recurrent(input);
			else
				return this->evaluate_nonrecurrent(input);
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
					unsigned int input_size, output_size, type; // 0 = ordinal, 1 = input, 2 = output
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

			unsigned int neuron_number;
			o << neuron_number << std::endl;					

			for (unsigned int i=0; i<neuron_number; i++){
				unsigned int input_size, output_size, type; // 0 = ordinal, 1 = input, 2 = output
				o << nodes[i].type << " ";

				if (type == 1)
					input_nodes.push_back(&(nodes[i]));
				if (type == 2)
					output_nodes.push_back(&(nodes[i]));

				o << nodes[i].in_nodes.size() << std::endl;
				for (unsigned int j=0; j<input_size; j++)
					o << nodes[i].in_nodes[i].first->number << " " 
						<< nodes[i].in_nodes[i].second << " ";
				o << std::endl;				
			}
			o.close();
		}	

	};

} // end of namespace ann

#endif
