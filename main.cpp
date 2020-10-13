#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "args.hxx"
#include "sdsl/bit_vectors.hpp"
#include "pg-pathcomp.hpp"
#include "xg.hpp"
#include "vg/vg.pb.h"
#include "vg/io/stream.hpp"

using namespace std;
using namespace xg;
using namespace vg;

int main(int argc, char** argv) {
    args::ArgumentParser parser("pg-pathcomp: compuate matrix of pairwise distances of graph");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "analyze the graph in this GFA file", {'g', "gfa-in"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "analyze the graph in this XG file", {'x', "xg-in"});
    args::ValueFlag<std::string> snarls_in(parser, "FILE", "snarls made with vg view -Fv graph.vg | vg snarls", {'r', "snarls-in"});
    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build from GFA", {'b', "base"});
    args::Flag validate(parser, "validate", "validate construction from GFA", {'V', "validate"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }

    XG graph;
    if (!args::get(xg_in).empty()) {
        cerr << "Loading XG" << endl;
        std::ifstream in(args::get(xg_in));
        graph.deserialize(in);
    } else if (!args::get(gfa_in).empty()) {
        cerr << "Importing GFA into XG" << endl;
        graph.from_gfa(args::get(gfa_in), args::get(validate),
                       args::get(base).empty() ? args::get(gfa_in) : args::get(base));
    }

    // load the snarls into a map
    map<handle_t, handle_t> snarl_map;
    if (!args::get(snarls_in).empty()) {
        std::ifstream in(args::get(snarls_in));
        vg::io::for_each<Snarl>(in, [&](Snarl& snarl) {
                handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
                handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
                snarl_map[start_handle] = end_handle;
                snarl_map[graph.flip(end_handle)] = snarl_map[graph.flip(start_handle)];
            });
    }
    
    cerr << "Computing Path Coverage" << endl;
    pair<vector<string>, vector<vector<float>>> path_comp = all_paths_jaccard(graph);

    // print header row (columns are in same order)
    cout << ".";
    for (size_t i = 0; i < path_comp.first.size(); ++i) {
        if (path_comp.first[i].substr(0, 4) != "HG01") continue;
        cout << "\t" << path_comp.first[i];
    }
    cout << endl;

    // print the matrix
    for (size_t i = 0; i < path_comp.second.size(); ++i) {
        if (path_comp.first[i].substr(0, 4) != "HG01") continue;
        cout << path_comp.first[i];
        for (size_t j = 0; j < path_comp.second[i].size(); ++j) {
            if (path_comp.first[j].substr(0, 4) != "HG01") continue;
            cout << "\t" << path_comp.second[i][j];
        }
        cout << endl;
    }    
    
    return 0;
}

