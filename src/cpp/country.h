#pragma once
#include "data.h"

#include <map>

class Country{
public:
    Country(const Args& args_, int N, int M): args(args_),
                                              agents(N), hypergraph(N,M, args_.random_seed),
                                              generator(args_.random_seed),
                                              p_inf(0.0,1.0){
        // TODO init SEIR
        set_random_generators();
        read_data();
    }

    void infection(std::array<int, SEIR_SIZE>& stats, bool second_wave, long double I_sum);
    void init_stats(std::array<int, SEIR_SIZE>& stats, bool verbose);
    void simulate();

    // INIT
    void init_agents(std::vector<int>& inf_nodes, int inf_node_num);
    // Read input
    void read_data();
    void read_args(int argc, char* argv[]);

    // Helper functions
    int handle_S(unsigned int agent);
    int handle_E(unsigned int agent);
    int handle_I(unsigned int agent);

    void set_random_generators(){
        E_time = std::geometric_distribution<int>(args.sigma);
        I_time = std::geometric_distribution<int>(args.mu);
    }

    Args args;
private:
    AgentData agents;
    Hypergraph hypergraph;

    std::mt19937 generator;
    std::geometric_distribution<int> E_time;
    std::geometric_distribution<int> I_time;
    std::uniform_real_distribution<long double> p_inf;

};
