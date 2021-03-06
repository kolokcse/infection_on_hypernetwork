#include <iostream>
#include <vector>

#include "country.h"

void read_args(int argc, char* argv[], Args& args){
    for(int i=0;i<argc;i++){
        std::string act_param = argv[i];
        if(act_param=="--beta") args.beta = std::stod(argv[++i]);
        else if(act_param=="--second_launch") args.second_launch = std::stod(argv[++i]);
        else if(act_param=="--seed") args.random_seed = std::stoi(argv[++i]);
        else if(act_param=="--sigma") args.sigma = std::stod(argv[++i]);
        else if(act_param=="--mu") args.mu = std::stod(argv[++i]);
        else if(act_param=="--max_sim") args.max_sim = std::stoi(argv[++i]);
        else if(act_param=="--awareness_treshold") args.awareness_treshold = std::stoi(argv[++i]);
        else if(act_param=="--verbose"){ args.verbose=true;++i;}
    }
}
// TODO:
//     * continuous infecting
//     * super_infected
int main(int argc, char* argv[]){
    Args args;
    unsigned int N,M;
    std::cin>>N>>M;
    read_args(argc, argv, args);
    Country country(args, N, M);
    country.simulate();
}