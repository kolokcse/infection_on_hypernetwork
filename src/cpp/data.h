#pragma once

#include <vector>
#include <random>
#include <array>

struct Args{
    double beta;
    double mu, sigma;
    double second_launch=1.0;
    int random_seed;
    int max_sim;
    int awareness_treshold;
    bool verbose = false;
};

enum SEIR{S,E,I,R,SEIR_SIZE};
enum COLS {location, seed, states, timers, home, work_location, SIZE};

class AgentData{
public:
    AgentData(unsigned int N_): N(N_){
        rows[location].resize(N);
        rows[home].resize(N);
        rows[work_location].resize(N, -1);

        rows[seed].resize(N, -1);
        rows[states].resize(N, 0);
        rows[timers].resize(N, 0);
    }
    
    std::vector<int>& operator[](COLS col){
        return rows[col];
    }
    unsigned int N;
private:
    std::vector<int> rows[COLS::SIZE];
};

class Hypergraph{
    friend class Country;
public:
    Hypergraph(int N, int M, int random_seed): nodeNum(N), edgeNum(M), SEIR(M), infection_prob(M),
                    generator(random_seed),contactnumber(N), edgeI(M),
                    neigh_edge_indexes(N, std::vector<int>()), neigh_node_indexes(M, std::vector<int>())
                    {
        
    }



    int& get_neigh_edge(int node, int index){
        return neigh_edge_indexes[node][index];
    }
    
    int& get_neigh_node(int edge, int index){
        return neigh_node_indexes[edge][index];
    }

    void resize_node(int node, int degree){
        neigh_edge_indexes[node].resize(degree);
    }
    
    void resize_edge(int edge, int size){
        neigh_node_indexes[edge].resize(size);
    }

    unsigned int nodeNum;
    unsigned int edgeNum;
    std::vector<std::array<int,(int)SEIR_SIZE>> SEIR;
    std::vector<long double> infection_prob;
private:
    std::mt19937 generator;
    std::vector<int> contactnumber;
    std::vector<int> edgeI;
    std::vector<std::vector<int>> neigh_edge_indexes;
    std::vector<std::vector<int>> neigh_node_indexes;
};