#ifndef __TINYNEAT_HPP__
#define __TINYNEAT_HPP__

/* custom defines:
 * INCLUDE_ENABLED_GENES_IF_POSSIBLE  - if during experiment you found that too many genes are
 *                                      disabled, you can use this option.
 * ALLOW_RECURRENCY_IN_NETWORK	      - allowing recurrent links 
 *
 * GIVING_NAMES_FOR_SPECIES           - giving species unique names (need a dictionary with 
 *                                      names in a file "specie_names.dict"
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <random>
#include <map>
#include <algorithm>
#include <list>
#include <string>

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

	mutation_rate_container read_mutation_rates(std::ifstream& o);
	void write_mutation_rates(std::ofstream& o, mutation_rate_container& rates, std::string prefix);

	class genome {
	private:
		genome(){};

	public:
		unsigned int fitness = 0;
		unsigned int adjusted_fitness = 0;
		unsigned int global_rank = 0;

		unsigned int max_neuron;

		mutation_rate_container mutation_rates;		
		std::vector<gene> genes;

		genome(unsigned int functional_nodes, mutation_rate_container& rates){
			max_neuron = functional_nodes;
			mutation_rates = rates;
		}
		
		genome(const genome&) = default;
	};


	/* a specie is group of genomes which differences is smaller than some threshold */
	typedef struct {
		unsigned int top_fitness = 0;
		unsigned int average_fitness = 0;
		unsigned int staleness = 0;

	#ifdef GIVING_NAMES_FOR_SPECIES
		std::string name;
	#endif
		std::vector<genome> genomes;
	} specie;	


	/* a small world, where individuals (genomes) are making babies and evolving,
	 * becoming better and better after each generation :)
	 */
	class pool {
	private:
		pool(){};

		/* important part, only accecible for friend */
		unsigned int innovation_number = 0;	
		unsigned int new_innovation(){
			return ++innovation_number;
		};		

		unsigned int generation_number = 1;

	public:
		/* pool parameters */
		unsigned int max_fitness = 0;

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

			// create a basic generation with default genomes
			for (unsigned int i = 0; i<this->population; i++){
				genome new_genome(this->functional_nodes, this->mutation_rates);
				this->mutate(new_genome);
				this->add_to_species(new_genome);
			}

		}

		/* innovation tracking in current generation, should be cleared after each generation */
		std::map<std::pair<unsigned int, unsigned int>, unsigned int> track;

		/* species */
		std::list<specie> species;
		
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

		/* specie ranking */
		void rank_globally();
		void calculate_average_fitness(specie& s);
		unsigned int total_average_fitness();

		/* evolution */
		void cull_species(bool cut_to_one);
		genome breed_child(specie& s);
		void remove_stale_species();
		void remove_weak_species();
		void add_to_species(genome& child);

	public:
		/* next generation */
		void new_generation();
		unsigned int generation() { return this->generation_number; }

		/* calculate fitness */
		std::vector<std::pair<specie*, genome*>> get_genomes(){
			std::vector<std::pair<specie*, genome*>> genomes;
			for (auto s = this->species.begin(); s != this->species.end(); s++)
				for (size_t i=0; i<(*s).genomes.size(); i++)
					genomes.push_back(std::make_pair(&(*s), &((*s).genomes[i])));
			return genomes;
		}

		/* import and export */
		void import_fromfile(std::string filename);
		void export_tofile(std::string filename);		
	};
	

	/* now the evolutionary functions itself */


	genome pool::crossover(const genome& g1, const genome& g2){
		// Make sure g1 has the higher fitness, so we will include only disjoint/excess 
		// genes from the first genome.
		// If the fitness is equal then we will include from both genomes.
		bool include_from_both = (g2.fitness == g1.fitness);				
		if (g2.fitness > g1.fitness)
			return crossover(g2, g1);
		
		genome child(this->functional_nodes, this->mutation_rates);

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
				g.genes[i].weight += real_distributor(this->generator) * step * 2.0 - step;
			else
				g.genes[i].weight =  real_distributor(this->generator)*4 - 2.0; 
			
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
		/* network encoding:
		 * | input nodes | bias | output nodes |
		 */		
		auto is_input = [&](unsigned int node) -> bool {
				return node < this->input_size; };
		auto is_output = [&](unsigned int node) -> bool {
				return node < this->functional_nodes && node >= 
					(this->input_size + this->bias_size); };
		auto is_bias = [&](unsigned int node) -> bool {
				return node < (this->input_size + this->bias_size) && node >= this->input_size; };

		std::uniform_int_distribution<unsigned int> distributor1(0, g.max_neuron-1);
		unsigned int neuron1 = distributor1(this->generator);

		std::uniform_int_distribution<unsigned int> distributor2
			(this->input_size + this->bias_size, g.max_neuron-1);
		unsigned int neuron2 = distributor2(this->generator);
			
		if (is_output(neuron1) && is_output(neuron2))
			return;
		
		if (is_output(neuron1))
			std::swap(neuron1, neuron2);		

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
				(this->input_size, this->input_size + this->output_size-1);
			new_gene.from_node = bias_choose(this->generator);			
		}

		// if genome already has this connection
		for (size_t i=0; i<g.genes.size(); i++)
			if (g.genes[i].from_node == neuron1 && g.genes[i].to_node == neuron2)
				return ;
		
		// add new innovation if needed
		auto it = this->track.find(std::make_pair(neuron1, neuron2));
		if (it == this->track.end())
			new_gene.innovation_num = 
				this->track[std::make_pair(neuron1, neuron2)] = this->new_innovation();
		else
			new_gene.innovation_num = (*it).second;

		// mutate new link
		std::uniform_real_distribution<double> weight_generator(0.0, 1.0);
		new_gene.weight = weight_generator(this->generator) * 4.0 - 2.0;		
		
		g.genes.push_back(new_gene);
	}	


	void pool::mutate_node(genome& g){
		if (g.genes.size() == 0)
			return ;

		g.max_neuron++;

		// randomly choose a gene to mutate
		std::uniform_int_distribution<unsigned int> distributor(0, g.genes.size()-1);
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
		new_gene2.to_node = g.genes[gene_id].to_node;
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


	void pool::rank_globally(){
		std::vector<genome*> global;
		for (auto s = this->species.begin(); s != this->species.end(); s++)
			for (size_t i=0; i < (*s).genomes.size(); i++)
				global.push_back(&((*s).genomes[i]));			

		std::sort(global.begin(), global.end(), [](genome*& a, genome*& b) -> bool {
					return a->fitness < b->fitness; 
				});
		for (size_t j=0; j<global.size(); j++)
			global[j]->global_rank = j+1;
	}	

	void pool::calculate_average_fitness(specie& s){
		unsigned int total = 0;
		for (size_t i=0; i<s.genomes.size(); i++)
			total += s.genomes[i].global_rank;
		s.average_fitness = total / s.genomes.size();
	}

	unsigned int pool::total_average_fitness(){
		unsigned int total = 0;
		for (auto s = this->species.begin(); s != this->species.end(); s++)
			total += (*s).average_fitness;
		return total;
	}


	void pool::cull_species(bool cut_to_one) {
		for (auto s = this->species.begin(); s != this->species.end(); s++) {
			std::sort((*s).genomes.begin(), (*s).genomes.end(),
					[](genome& a, genome& b){ return a.fitness < b.fitness; });

			unsigned int remaining = std::ceil((*s).genomes.size() * 1. / 2.);
			if (cut_to_one)
				remaining = 1;
			while ((*s).genomes.size() > remaining)
				(*s).genomes.pop_back();
		}
	}

	genome pool::breed_child(specie &s){
		genome child(this->functional_nodes, this->mutation_rates);
		std::uniform_real_distribution<double> distributor(0.0, 1.0);
		std::uniform_int_distribution<unsigned int> choose_genome(0, s.genomes.size()-1);
		if (distributor(this->generator) < this->mutation_rates.crossover_chance){
			genome& g1 = s.genomes[choose_genome(this->generator)];
			genome& g2 = s.genomes[choose_genome(this->generator)];
			// QUESTION: if g1 == g2, then you can make a baby by fapping?			
			child = this->crossover(g1, g2);
		}
		else 
		{
			genome& g = s.genomes[choose_genome(this->generator)];
			child = g;
		}

		this->mutate(child);
		return child;
	
	}

	void pool::remove_stale_species(){
		auto s = this->species.begin();
		while (s != this->species.end()){
			genome& g = *(std::max_element((*s).genomes.begin(), (*s).genomes.end(), 
					[](genome& a, genome& b) -> bool { return a.fitness < b.fitness; }));
								
			if (g.fitness > (*s).top_fitness){
				(*s).top_fitness = g.fitness;
				(*s).staleness = 0;
			}
			else
				(*s).staleness++;

			if (! ((*s).staleness < this->stale_species || (*s).top_fitness >= this->max_fitness))
				this->species.erase(s++);
			else
				s++;			
		}
	}

	void pool::remove_weak_species(){
		unsigned int sum = this->total_average_fitness();
		auto s = this->species.begin();
		while (s != this->species.end()){
			double breed = std::floor((1. * (*s).average_fitness) / (1. * sum * this->population));
			if (breed >= 1.0)
				s++;
			else
				this->species.erase(s++);
		}
	}

	void pool::add_to_species(genome& child){		
		auto s = this->species.begin();
		while (s != this->species.end()){
			if (this->is_same_species(child, (*s).genomes[0])){
				(*s).genomes.push_back(child);
				break;
			}
			++s;
		}

		if (s == this->species.end()){
			specie new_specie;
			new_specie.genomes.push_back(child);
			this->species.push_back(new_specie);
		}										
	}

	void pool::new_generation(){
		this->cull_species(false);			
		this->rank_globally();
		this->remove_stale_species();

		for (auto s = this->species.begin(); s != this->species.end(); s++)
			this->calculate_average_fitness((*s));
		this->remove_weak_species();
		
		std::vector<genome> children;
		unsigned int sum = this->total_average_fitness();
		for (auto s = this->species.begin(); s != this->species.end(); s++){
			unsigned int breed = 
				std::floor((1. * (*s).average_fitness) / (1. * sum * this->population));
			for (unsigned int i = 0; i < breed; i++)
				children.push_back(this->breed_child(*s));
		}

		this->cull_species(true); // now in each species we have only one genome

		// preparing for MAKING BABIES <3		
		std::uniform_int_distribution<unsigned int> choose_specie(0, this->species.size()-1);
		std::vector<specie*> species_pointer;
		for (auto s = this->species.begin(); s != this->species.end(); s++)
			species_pointer.push_back(&(*s));

		while (children.size() + this->species.size() < this->population)
			children.push_back(this->breed_child(*species_pointer[choose_specie(this->generator)]));
	
		for (size_t i=0; i<children.size(); i++)
			this->add_to_species(children[i]);				

		this->generation_number++;
	}	

	void pool::import_fromfile(std::string filename){
		std::ifstream input;
		input.open(filename);
		if (!input.is_open()){
			std::cerr << "cannot open file '" << filename << "' !";
			return ;
		}
	
		this->species.clear();
		try {
			// current state
			input >> this->innovation_number;
			input >> this->generation_number;
			input >> this->max_fitness;

			// network information
			input >> this->input_size >> this->output_size >> this->bias_size;	
			this->functional_nodes = input_size + output_size + bias_size;			

			// population information
			input >> this->population;
			input >> this->delta_disjoint;
			input >> this->delta_weights;
			input >> this->delta_threshold;
			input >> this->stale_species;

			// mutation parameters
			this->mutation_rates = read_mutation_rates(input);	

			// species information
			unsigned int species_number;
			input >> species_number;
			this->species.clear();			

			for (unsigned int c = 0; c < species_number; c++){
				specie new_specie;
			#ifdef GIVING_NAMES_FOR_SPECIES
				input >> new_specie.name;
			#endif
				input >> new_specie.top_fitness;
				input >> new_specie.average_fitness;
				input >> new_specie.staleness;

				unsigned int specie_population;
				input >> specie_population;

				for (unsigned int i=0; i<specie_population; i++){
					genome new_genome(this->functional_nodes, this->mutation_rates);
					input >> new_genome.fitness;
					input >> new_genome.adjusted_fitness;
					input >> new_genome.global_rank;
					
					new_genome.mutation_rates = read_mutation_rates(input);

					unsigned int gene_number;
					input >> new_genome.max_neuron >> gene_number;

					for (unsigned int j=0; j<gene_number; j++){
						gene new_gene;
						input >> new_gene.innovation_num;
						input >> new_gene.from_node;
						input >> new_gene.to_node;
						input >> new_gene.weight;
						input >> new_gene.enabled;
						new_genome.genes.push_back(new_gene);
					}
					
					new_specie.genomes.push_back(new_genome);
				}

				this->species.push_back(new_specie);
			}
			
		}
		catch (std::string error_message){
			std::cerr << error_message;
		}

		input.close();
	}

	void pool::export_tofile(std::string filename){
		std::ofstream output;
		output.open(filename);
		if (!output.is_open()){
			std::cerr << "cannot open file '" << filename << "' !";
			return ;
		}

		// current state
		output << this->innovation_number << std::endl;
		output << this->generation_number << std::endl;
		output << this->max_fitness << std::endl;

		// network information
		output << this->input_size << " " <<  this->output_size << " " <<
		   	this->bias_size << std::endl;	
		this->functional_nodes = input_size + output_size + bias_size;

		// population information
		output << this->population << std::endl;
		output << this->delta_disjoint << std::endl;
		output << this->delta_weights << std::endl;
		output << this->delta_threshold << std::endl;
		output << this->stale_species << std::endl;

		// mutation parameters
		write_mutation_rates(output, this->mutation_rates, "");

		// species information
		output << this->species.size() << std::endl;
		for (auto s = this->species.begin(); s != this->species.end(); s++){
			output << "   ";
		#ifdef GIVING_NAMES_FOR_SPECIES
			output << (*s).name << " ";
		#endif
			output << (*s).top_fitness << " ";
			output << (*s).average_fitness << " ";
			output << (*s).staleness << std::endl;

			output << "   " << (*s).genomes.size() << std::endl;
			for (size_t i=0; i<(*s).genomes.size(); i++){
				write_mutation_rates(output, (*s).genomes[i].mutation_rates, "      ");
				output << "      ";
			    output << (*s).genomes[i].fitness << " ";
				output << (*s).genomes[i].adjusted_fitness << " ";
				output << (*s).genomes[i].global_rank << std::endl;
				
				output << "      " << (*s).genomes[i].max_neuron << " " <<
				   	(*s).genomes[i].genes.size() << std::endl;
				for (size_t j=0; j<(*s).genomes[i].genes.size(); j++){
					gene& g = (*s).genomes[i].genes[j];
					output << "         ";
					output << g.innovation_num << " " << g.from_node << " " << g.to_node << " "
						<< g.weight << " " << g.enabled << std::endl;
				}
			}

			output << std::endl << std::endl;
		}


	}

	mutation_rate_container read_mutation_rates(std::ifstream& o){
		mutation_rate_container result;
		o >> result.connection_mutate_chance;
		o >> result.perturb_chance;
		o >> result.crossover_chance;
		o >> result.link_mutation_chance;
		o >> result.node_mutation_chance;
		o >> result.bias_mutation_chance;
		o >> result.step_size;
		o >> result.disable_mutation_chance;
		o >> result.enable_mutation_chance;
		return result;
	}

	void write_mutation_rates(std::ofstream& o, mutation_rate_container& rates, std::string prefix){
		o << prefix << rates.connection_mutate_chance << std::endl;
		o << prefix << rates.perturb_chance << std::endl;
		o << prefix << rates.crossover_chance << std::endl;
		o << prefix << rates.link_mutation_chance << std::endl;
		o << prefix << rates.node_mutation_chance << std::endl;
		o << prefix << rates.bias_mutation_chance << std::endl;
		o << prefix << rates.step_size << std::endl;
		o << prefix << rates.disable_mutation_chance << std::endl;
		o << prefix << rates.enable_mutation_chance << std::endl;
	}	

} // end of namespace neat

#endif
