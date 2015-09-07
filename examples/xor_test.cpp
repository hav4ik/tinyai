#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "../src/tinyneat.hpp"
#include "../src/tinyann.hpp"

// returns the fitness.
unsigned int xor_test(ann::neuralnet& n){
	std::vector<double> input(2, 0.0);
	std::vector<double> output(1, 0.0);
	unsigned int fitness = 0;
	double score;
	double answer;

	std::cerr << "     > begin xor test" << std::endl;

	input[0] = 0.0, input[1] = 0.0, answer = 0.0;
	n.evaluate(input, output);
	fitness += std::min(std::abs(1.0 / std::abs(answer - output[0])), 2.0);

	input[0] = 0.0, input[1] = 1.0, answer = 1.0;
	n.evaluate(input, output);
	fitness += std::min(std::abs(1.0 / std::abs(answer - output[0])), 2.0);

	input[0] = 1.0, input[1] = 0.0, answer = 1.0;
	n.evaluate(input, output);
	fitness += std::min(std::abs(1.0 / std::abs(answer - output[0])), 2.0);

	input[0] = 1.0, input[1] = 1.0, answer = 0.0;
	n.evaluate(input, output);
	fitness += std::min(std::abs(1.0 / std::abs(answer - output[0])), 2.0);

	return fitness;
}

int main(){
	neat::pool p(2, 1, 1);
//	ann::neuralnet n;
	srand(time(NULL));
	unsigned int max_fitness = 0;
	while (max_fitness < 8){
		unsigned int current_fitness = 0;		
		for (auto s = p.species.begin(); s != p.species.end(); s++)
			for (size_t i=0; i<(*s).genomes.size(); i++){
				ann::neuralnet n;
				neat::genome& g = (*s).genomes[i];
				std::cerr << " ++ generating neural network from genome..." << std::endl;
				n.from_genome(g);
				std::cerr << " -- generated." << std::endl;
				std::cerr << " ++ doing the xor test..." << std::endl;
				current_fitness = xor_test(n);
				std::cerr << " -- fitness calculated! (" << current_fitness << ")" << std::endl;
				if (current_fitness > max_fitness){
					max_fitness = current_fitness;
					std::string fname = "fit = " + std::to_string(current_fitness);
//					std::cerr << " ++ begin exporting to file" << std::endl;
//					n.export_tofile(fname);
//					std::cerr << " -- exported to file. New record: " << current_fitness << std::endl;
				}
				g.fitness = current_fitness;
				std::cerr << " ** not record, move to next :3" << std::endl;
			}
		std::cerr << " NEW GENERATION ..." << std::endl;
		p.new_generation();
		std::cerr << " ... CREATED! " << std::endl;
	}

	return 0;
}
