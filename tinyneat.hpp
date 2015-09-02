#ifndef __TINYNEAT_HPP__
#define __TINYNEAT_HPP__

/* custom defines:
 * INCLUDE_ENABLED_GENES_IF_POSSIBLE  - if during experiment you found that too many genes are
 *                                      disabled, you can use this option.
 * ALLOW_RECURRENCY_IN_NETWORK	      - allowing recurrent links 
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <random>
#include <map>

namespace neat {

	typedef struct {	
		unsigned int innovation_num = -1;
		unsigned int from_node = -1;
		unsigned int to_node = -1;
		double weight = 0.0;
		bool enabled = true;
	} gene;
	
	typedef struct {		
		double connection_mutate_chance = 0.25;
		double perturb_chance = 0.90;
		double crossover_chance = 0.75;
		double link_mutation_chance = 2.0;
		double node_mutation_chance = 0.50;
		double bias_mutation_chance = 0.40;
		double step_size = 0.1;
		double disable_mutation_chance = 0.4;
		double enable_mutation_chance = 0.2;			
	} mutation_rate_container;

	/* predefinition */
	class genome;

	/* a specie is group of genomes which differences is smaller than some threshold */
	typedef struct {
		unsigned int top_fitness = 0;
		unsigned int average_fitness = 0;
		std::vector<genome> genomes;
		double staleness = 0.0;		
	} specie;	


	/* a small world, where individuals (genomes) are making babies and evolving,
	 * becoming better and better after each generation :)
	 */
	class pool {
	private:

		/* important part, only accecible for friend */
		unsigned int innovation_number = 0;
		unsigned int new_innovation(){
			return ++innovation_number;
		};


	public:

		/* mutation parameters */
		mutation_rate_container mutation_rates;

		/* species parameters */
		unsigned int population = 300;
		double delta_disjoint = 2.0;
		double delta_weights = 0.4;
		double delta_threshold = 1.0;
		unsigned int stale_species = 15;		

		/* neural network parameters */
		unsigned int input_size;
		unsigned int bias_size;
		unsigned int output_size;
		unsigned int functional_nodes;		
		
		// pool's local random number generator
		std::random_device rd;		
		std::mt19937 generator;

		// constructor
		pool(unsigned int input, unsigned int output, unsigned int bias = 1){
			input_size = input;
			output_size = output;
			bias_size = bias;
			functional_nodes = input_size + output_size + bias_size;

			// seed the mersenne twister with 
			// a random number from our computer
			generator.seed(rd());	

		}

		/* innovation tracking in current generation, should be cleared after each generation */
		std::map<std::pair<unsigned int, unsigned int>, unsigned int> track;

		/* species */
		std::vector<specie> species;
		
	private:
		/* evolutionary methods */
		genome crossover(const genome& g1, const genome& g2);
		void mutate_weight(genome& g);
		void mutate_enable_disable(genome& g, bool enable);
		void mutate_link(genome& g, bool force_bias);
		void mutate_node(genome& g);		
		void mutate(genome& g);		

		double disjoint(const genome& g1, const genome& g2);
		double weights(const genome& g1, const genome& g2);
		bool is_same_species(const genome& g1, const genome& g2);
	};


	class genome {
	public:
		unsigned int fitness = 0;
		unsigned int adjusted_fitness = 0;
		unsigned int global_rank = 0;

		unsigned int max_neuron;

		mutation_rate_container mutation_rates;
		
		std::vector<gene> genes;
		pool* global_pool;

		genome() = delete;
		genome(pool* glob_pool){			
			global_pool = glob_pool;
			max_neuron = global_pool->functional_nodes;
			mutation_rates = global_pool->mutation_rates;
		}
		genome(const genome&) = default;
		
	};
	
	class species {
	public:
		unsigned int top_fitness = 0;
		unsigned int stale_fitness = 0;
		unsigned int average_fitness = 0;
		std::vector<genome> genomes;

		species(){}
	};



	/* now the evolutionary functions itself */


	genome pool::crossover(const genome& g1, const genome& g2){
		// Make sure g1 has the higher fitness, so we will include only disjoint/excess 
		// genes from the first genome.
		// If the fitness is equal then we will include from both genomes.
		bool include_from_both = (g2.fitness == g1.fitness);				
		if (g2.fitness > g1.fitness)
			return crossover(g2, g1);
		
		genome child(this);

		auto it1 = g1.genes.begin();
		auto it2 = g2.genes.begin();

		// coin flip random number distributor
		std::uniform_int_distribution<int> coin_flip(1, 2);	

		for (; it1 != g1.genes.end() && it2 != g2.genes.end(); it1++){			

			// move forward if not match and include genes from the second genome
			// if their fitness are equal
			while ((*it2).innovation_num < (*it1).innovation_num && it2 != g2.genes.end()){
				if (include_from_both)
					child.genes.push_back(*it2);
				it2++;
			}			
		
			// if innovation marks match, do the crossover, else include from the first
			// genome because its fitness is not smaller than the second's 
			
			if ((*it2).innovation_num == (*it1).innovation_num){
				// do the coin flip
				int coin = coin_flip(this->generator);
				

			// now, after flipping the coin, we do the crossover.
			#ifdef INCLUDE_ENABLED_GENES_IF_POSSIBLE
				if (coin == 2 && (*it2).enabled)
					child.genes.push_back(*it2);
				else
					child.genes.push_back(*it1);
			#else
				if (coin == 2)
					child.genes.push_back(*it2);
				else
					child.genes.push_back(*it1);
			#endif

			} else 			

				// as said before, we include from the first (with larger fitness) otherwise
				child.genes.push_back(*it1);

		}

		child.max_neuron = std::max(g1.max_neuron, g2.max_neuron);		
		return child;	
	}


	/* mutations */
	void pool::mutate_weight(genome& g){		
		double step = this->mutation_rates.step_size;
		std::uniform_real_distribution<double> real_distributor(0.0, 1.0);

		for (size_t i=0; i<g.genes.size(); i++){
			if (real_distributor(this->generator) < this->mutation_rates.perturb_chance) 
				g.genes[i].weight += real_distributor(this->generator) * step * 2 - step;
			else
				g.genes[i].weight =  real_distributor(this->generator)*4 - 2; 
			
		}
	}

	void pool::mutate_enable_disable(genome& g, bool enable){
		std::vector<gene*> v;
		for (size_t i=0; i<g.genes.size(); i++)
			if (g.genes[i].enabled != enable)
				v.push_back(&(g.genes[i]));

		if (v.size() == 0)
			return ;

		std::uniform_int_distribution<int> distributor(0, v.size()-1);
		v[distributor(this->generator)]->enabled = enable;
	}

	void pool::mutate_link(genome& g, bool force_bias){
		auto is_input = [&](unsigned int node) -> bool {
				return node < this->input_size; };
		auto is_output = [&](unsigned int node) -> bool {
				return node < (this->input_size + this->output_size) && node >= this->input_size; };
		auto is_bias = [&](unsigned int node) -> bool {
				return node < this->functional_nodes && node >= 
					(this->input_size + this->output_size); };

		std::uniform_int_distribution<unsigned int> distributor1(0, g.max_neuron-1);
		unsigned int neuron1 = distributor1(this->generator);

		std::uniform_int_distribution<unsigned int> distributor2
			(this->functional_nodes, g.max_neuron-1);
		unsigned int neuron2 = distributor2(this->generator);

		if (is_input (neuron1) && is_input (neuron2)) return;
		if (is_output(neuron1) && is_output(neuron2)) return;
		if (is_bias  (neuron1) && is_bias  (neuron2)) return;

	#ifndef ALLOWING_RECURRENCY_IN_NETWORK
		// check for recurrency using BFS
		bool has_recurrence = false;
		if (is_input(neuron1))
			has_recurrence = false;
		else {
			std::queue<size_t> que;
			std::vector<std::vector<unsigned int>> connections(g.max_neuron);
			for (size_t i=0; i<g.genes.size(); i++)
				connections[g.genes[i].from_node].push_back(g.genes[i].to_node);
		
			for (size_t i=0; i<connections[neuron1].size(); i++)
				que.push(connections[neuron1][i]);

			while (!que.empty()){				
				unsigned int tmp = que.front();
				if (tmp == neuron1){
					has_recurrence = true;
					break;
				}
				que.pop();
				for (size_t i=0; i<connections[tmp].size(); i++)
					que.push(connections[tmp][i]);				
			}
		}
		if (has_recurrence)
			return ;
	#endif
		
		// now we can create a link 
		gene new_gene;
		new_gene.from_node = neuron1;
		new_gene.to_node = neuron2;
	
		if (force_bias){
			std::uniform_int_distribution<unsigned int> bias_choose
				(this->input_size+this->output_size, this->functional_nodes-1);
			new_gene.from_node = bias_choose(this->generator);			
		}

		// if genome already has this connection
		for (size_t i=0; i<g.genes.size(); i++)
			if (g.genes[i].from_node == neuron1 && g.genes[i].to_node == neuron2)
				return ;
		
		// add new innovation if needed
		if (this->track.find(std::make_pair(neuron1, neuron2)) == this->track.end())
			new_gene.innovation_num = 
				this->track[std::make_pair(neuron1, neuron2)] = this->new_innovation();

		// mutate new link
		std::uniform_real_distribution<double> weight_generator(0.0, 1.0);
		new_gene.weight = weight_generator(this->generator);
		
		g.genes.push_back(new_gene);
	}	


	void pool::mutate_node(genome& g){
		if (g.genes.size() == 0)
			return ;

		g.max_neuron++;

		// randomly choose a gene to mutate
		std::uniform_int_distribution<unsigned int> distributor(0, g.genes.size());
		unsigned int gene_id = distributor(this->generator);
		if (g.genes[gene_id].enabled == false)
			return ;

		g.genes[gene_id].enabled = false;

		gene new_gene1;
		new_gene1.from_node = g.genes[gene_id].from_node;
		new_gene1.to_node = g.max_neuron-1; // to the last created neuron
		new_gene1.weight = 1.0;
		new_gene1.innovation_num = this->new_innovation();
		new_gene1.enabled = true;

		gene new_gene2;
		new_gene2.from_node = g.max_neuron-1; // from the last created neuron
		new_gene2.weight = g.genes[gene_id].weight;
		new_gene2.innovation_num = this->new_innovation();
		new_gene2.enabled = true;

		g.genes.push_back(new_gene1);
		g.genes.push_back(new_gene2);
	}

	void pool::mutate(genome& g){
		double coefficient[2] = {0.95, 1.05263};

		std::uniform_int_distribution<int> coin_flip(0, 1);		

		g.mutation_rates.enable_mutation_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.disable_mutation_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.connection_mutate_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.node_mutation_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.link_mutation_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.bias_mutation_chance *= coefficient[coin_flip(this->generator)];		
		g.mutation_rates.crossover_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.perturb_chance *= coefficient[coin_flip(this->generator)];

		std::uniform_real_distribution<double> mutate_or_not_mutate(0.0, 1.0);

		if (mutate_or_not_mutate(this->generator) < g.mutation_rates.connection_mutate_chance)
			this->mutate_weight(g);

		double p;

	    p = g.mutation_rates.link_mutation_chance;
		while (p > 0.0) {		
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_link(g, false);
			p = p - 1.0;
		}

		p = g.mutation_rates.bias_mutation_chance;
		while (p > 0.0) {		
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_link(g, true);
			p = p - 1.0;
		}

		p = g.mutation_rates.node_mutation_chance;
		while (p > 0.0) {		
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_node(g);
			p = p - 1.0;
		}

		p = g.mutation_rates.enable_mutation_chance;;
		while (p > 0.0) {		
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_enable_disable(g, true);
			p = p - 1.0;
		}

		p = g.mutation_rates.disable_mutation_chance;
		while (p > 0.0) {		
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_enable_disable(g, false);
			p = p - 1.0;
		}
		
	}
	


	
	double pool::disjoint(const genome& g1, const genome& g2){
		auto it1 = g1.genes.begin();
		auto it2 = g2.genes.begin();

		unsigned int disjoint_count = 0;
		for (; it1 != g1.genes.end() && it2 != g2.genes.end(); it1++){
			for (; (*it2).innovation_num < (*it1).innovation_num && it2 != g2.genes.end(); it2++)
				disjoint_count++;
	
			if ((*it1).innovation_num != (*it2).innovation_num)
				disjoint_count++;			
		}

		return (1. * disjoint_count) / (1. * std::max(g1.genes.size(), g2.genes.size()));
	}
	
	double pool::weights(const genome& g1, const genome& g2){
		auto it1 = g1.genes.begin();
		auto it2 = g2.genes.begin();

		double sum = 0.0;
		unsigned int coincident = 0;
		for (; it1 != g1.genes.end() && it2 != g2.genes.end(); it1++){
			while ((*it2).innovation_num < (*it1).innovation_num && it2 != g2.genes.end()) 
				it2++;

			
			if ((*it1).innovation_num == (*it2).innovation_num){
				coincident++;
				sum += std::abs((*it1).weight - (*it2).weight);
			}
		}

		return sum / (1. * coincident);
	}

	bool pool::is_same_species(const genome& g1, const genome& g2){
		double dd = this->delta_disjoint * disjoint(g1, g2);
		double dw = this->delta_weights * weights(g1, g2);
		return dd + dw < this->delta_threshold;
	}

		

} // end of namespace neat

#endif
