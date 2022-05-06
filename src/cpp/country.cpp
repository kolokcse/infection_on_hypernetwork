#include <iostream>
#include <vector>
#include <numeric>
#include <stdlib.h>
#include <assert.h>
#include <math.h> 

#include "country.h"

void Country::init_agents(std::vector<int>& inf_nodes, int inf_node_num){


    for(int i=0;i<inf_node_num;i++){
        agents[states][i] = I;
        agents[timers][i] = std::min(20,I_time(generator)+1);
    }
    //std::cout<<all_agent_num<<" "<<agents.N<<std::endl;
    assert(hypergraph.nodeNum == agents.N);
}


// === Read input ===
void Country::read_data(){
    // 1. Read nodes: intersecting edges
    for(unsigned int node=0;node<hypergraph.nodeNum;node++){

        unsigned int degree;
        std::cin>>degree;
        hypergraph.resize_node(node, degree+1);
        for(unsigned int j=0;j<degree;j++){
            std::cin>>hypergraph.get_neigh_edge(node, j);
        }
    }
    // 2. Read edges: intersecting nodes
    for(unsigned int edge=0; edge<hypergraph.edgeNum;edge++){

        unsigned int edgesize;
        std::cin>>edgesize;
        //unsigned int contactnumber;
        //std::cin>>contactnumber;
        hypergraph.resize_edge(edge, edgesize+1);
        for(unsigned int j=0;j<edgesize;j++){
            std::cin>>hypergraph.get_neigh_node(edge, j);
            hypergraph.SEIR[edge][S] = edgesize;
        }
    }

    // 2. Read infected nodes
    int inf_node_num;
    std::cin>>inf_node_num;
    std::vector<int> inf_nodes;
    inf_nodes.resize(inf_node_num+1);
    for(int i=0;i<inf_node_num;i++){
        int inf_node;
        std::cin>>inf_node;
        inf_nodes[i]=inf_node;
        std::vector<int> neigh_edges;
        neigh_edges = hypergraph.neigh_edge_indexes[inf_node];
        int degree=neigh_edges.size();
        for(int i=0;i<degree;i++){
            hypergraph.SEIR[neigh_edges[i]][S] -= 1;
            hypergraph.SEIR[neigh_edges[i]][I] += 1;
        }
    }

    // === Init data with args ===
    init_agents(inf_nodes, inf_node_num);
    //std::cout<<"Read cities\n";
}

/*
void Country::go_home(){
    for(unsigned int agent=0; agent<agents.N; agent++){
        if(agents[work_location][agent] != -1 &&
           agents[location][agent] == agents[work_location][agent]){
            int old_loc = agents[location][agent];
            int new_loc = agents[home][agent];
            // === Go home ===
            agents[location][agent] = new_loc;
            // Update SEIR in city
            SEIR state = (SEIR) agents[COLS::states][agent];
            graph.SEIR[old_loc][state] -= 1;
            graph.SEIR[new_loc][state] += 1;
        }
    }
}


void Country::move(){
    std::vector<int> new_loc(agents[location]);
    for(unsigned int agent=0; agent<agents.N; agent++){
        int old_loc = agents[location][agent];
        // 1. Worker Agents
        if(agents[work_location][agent] != -1){
            // If at work ==> go home;
            // If at home ==> go work
            if(agents[location][agent] == agents[work_location][agent]){
                new_loc[agent]=agents[home][agent];
            }
            else if(p_moving(generator)){
                new_loc[agent]=agents[work_location][agent];
            }
        }
        // 2. Traveller Agents
        else if(p_moving(generator)){
            int act_city = agents[location][agent];
            // 2. Travellers
            int move_to = graph.get_random_neigh(act_city);
            new_loc[agent] = move_to;
        }
        // Update SEIR in city
        SEIR state = (SEIR) agents[COLS::states][agent];
        graph.SEIR[old_loc][state] -= 1;
        graph.SEIR[new_loc[agent]][state] += 1;
    }

    for(unsigned int i=0;i<agents.N;i++){
        agents[location][i] = new_loc[i];
    }
}
*/

int Country::handle_S(unsigned int agent){
    int change = 0;
    long double p = 1;
    std::vector<int> neigh_edges;
    neigh_edges = hypergraph.neigh_edge_indexes[agent];
    int degree=neigh_edges.size();
    for(int i=0;i<degree;i++){
        p = p*(1-hypergraph.infection_prob[neigh_edges[i]]);
    }
    p = 1 - p;
    if(p_inf(generator) < p){
        agents[states][agent] = E;
        agents[timers][agent] = E_time(generator);
        change = 1;
        for(int i=0;i<degree;i++){
            hypergraph.SEIR[neigh_edges[i]][S] -= 1;
            hypergraph.SEIR[neigh_edges[i]][E] += 1;
        }
    }
    return change;
}

int Country::handle_E(unsigned int agent){
    int change = 0;
    if(agents[timers][agent] == 0){
        change = 1;
        agents[states][agent] = I;
        agents[timers][agent] = std::min(20, I_time(generator)+1);
        //agents[timers][agent] = std::min(10, I_time(generator)+1);

        std::vector<int> neigh_edges;
        neigh_edges = hypergraph.neigh_edge_indexes[agent];
        int degree=neigh_edges.size();
        for(int i=0;i<degree;i++){
            hypergraph.SEIR[neigh_edges[i]][E] -= 1;
            hypergraph.SEIR[neigh_edges[i]][I] += 1;
        }
    }
    else{
        agents[timers][agent]-=1;
    }
    return change;
}

int Country::handle_I(unsigned int agent){
    int change = 0;
    if(agents[timers][agent] == 0){
        change = 1;
        agents[states][agent] = R;
        agents[timers][agent] = 0;

        std::vector<int> neigh_edges;
        neigh_edges = hypergraph.neigh_edge_indexes[agent];
        int degree=neigh_edges.size();
        for(int i=0;i<degree;i++){
            hypergraph.SEIR[neigh_edges[i]][I] -= 1;
            hypergraph.SEIR[neigh_edges[i]][R] += 1;
        }
    }
    else{
        agents[timers][agent]-=1;
    }
    return change;
}

void operator+=(std::array<int,SEIR_SIZE>& base, const std::array<int,SEIR_SIZE>& data){
    for(int i=0;i<SEIR_SIZE;i++){
        base[(SEIR)i] += data[(SEIR) i];
    }
}

void operator<<(std::ostream& os, std::array<int,SEIR_SIZE>& base){
    for(int i=0;i<SEIR_SIZE;i++){
        os<<base[(SEIR)i]<<" ";
    }
    os<<std::endl;
}

void operator<<(std::ostream& os, const std::array<int,SEIR_SIZE>& base){
    for(int i=0;i<SEIR_SIZE;i++){
        os<<base[(SEIR)i]<<" ";
    }
    os<<";";
}

void Country::infection(std::array<int, SEIR_SIZE>& stats, bool second_wave, long double I_sum){
    // === Compute infection probability ===



    for(unsigned int edge=0;edge<hypergraph.edgeNum;edge++){
        long double I_edge = 0;
        int size = hypergraph.neigh_node_indexes[edge].size();
        I_edge = hypergraph.SEIR[edge][I];

        //long double contactNum = hypergraph.contactnumber[edge];
        long double contactNum = 1.0;
        long double p;

        //std::cout<<"I "<<I_city<<std::endl;

        if(I_edge == 0){
            p=0;
        }
        else{
            double beta = args.beta;
            //double beta = args.beta * exp(1-I_sum/args.K);
            p = (contactNum*beta*I_edge/(size-1));
            //p = args.beta*I_city/graph.population[city];
            //std::cout.precision(17);
            //std::cout<<p<<" "<<I_sum<<" "<<std::endl;
        }

        hypergraph.infection_prob[edge] = p;
    }

    // === Infect/heal agents ===
    for(unsigned int agent=0; agent<agents.N; agent++){

        int change = 0;
        if(agents[states][agent] == S){
            change = handle_S(agent);
            stats[S]-=change;
            stats[E]+=change;
        }

        if(agents[states][agent] == E){
            change = handle_E(agent);
            stats[E]-=change;
            stats[I]+=change;
        }

        if(agents[states][agent] == I){
            change = handle_I(agent);
            stats[I]-=change;
            stats[R]+=change;
        }
    }
}

void debug(int iteration_num, AgentData& agents){
    std::cout<<iteration_num<<". iter\n";
    for(unsigned int i=0;i<agents.N;i++){
        std::cout<<agents[states][i];
    }
    std::cout<<std::endl;
}

void log_edges(const Hypergraph& hypergraph){
    for(unsigned int edge=0;edge<hypergraph.edgeNum;edge++){
        std::cout<<hypergraph.SEIR[edge];
    }
    std::cout<<std::endl;
}

void Country::init_stats(std::array<int,SEIR_SIZE>& base, bool verbose){
    for(unsigned int agent=0; agent<agents.N; agent++){
        base[agents[states][agent]]+=1;
    }
    std::cout<<base;
    if(verbose) log_edges(hypergraph);
}

void Country::simulate(){
    // === INIT ===
    std::vector<std::array<int, SEIR_SIZE>> stats(args.max_sim+1, {0,0,0,0});
    init_stats(stats[0], args.verbose);

    for(int iteration_num=1;iteration_num<=args.max_sim; iteration_num++){
        // 1. Move agents
        //move();

        // 2. Infect agents
        bool second_wave = (stats[iteration_num-1][I]+stats[iteration_num-1][R] > agents.N*args.second_launch);
        stats[iteration_num] = stats[iteration_num-1];
        infection(stats[iteration_num], second_wave, stats[iteration_num-1][I]);

        // 3. Go home, and infect there too: TODO: take care for travaller agents
        //go_home();
        //infection(stats[iteration_num], second_wave);

        // === Log ===
        std::cout<<stats[iteration_num];
        if(args.verbose) log_edges(hypergraph);
	    if(stats[iteration_num][I]==0 && stats[iteration_num][E]==0) break;
    }
}
