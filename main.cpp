#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "args.hxx"
#include "sdsl/bit_vectors.hpp"
#include "pg-pathcomp.hpp"
#include "xg.hpp"
#include "vg.pb.h"
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
    args::ValueFlag<std::string> ref_path(parser, "NAME", "reference path name for SV pileup", {'p', "ref-path"});
    args::ValueFlag<std::string> mat_out(parser, "FILE", "write path similarity matrix to this file", {'m', "mat-out"});
    args::ValueFlag<std::string> pileup_out(parser, "FILE", "write sv pileup (WIG format) to this file", {'u', "pileup-out}"});
    args::ValueFlag<std::string> coverage_stats(parser, "NAME", "print coverage stats for given prefix", {'s', "coverage-stats"});
    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build from GFA", {'b', "base"});
    args::Flag validate(parser, "validate", "validate construction from GFA", {'V', "validate"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<uint64_t> min_sv_length(parser, "M", "min length for SV to be included in pileup [10]", {'M', "min-sv-len"});
    args::ValueFlag<uint64_t> window_size(parser, "w", "window size for pileup [1]", {'w', "window-size"});
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

    size_t min_sv_len = args::get(min_sv_length);
    if (!min_sv_length) {
        min_sv_len = 10;
    }

    size_t window = args::get(window_size);
    if (!window) {
        window = 5;
    }

    // todo: we don't actually support window yet
    window = 1;

    XG graph;
    if (!args::get(xg_in).empty()) {
        cerr << "Loading XG" << endl;
        std::ifstream in(args::get(xg_in));
        graph.deserialize(in);
    } else if (!args::get(gfa_in).empty()) {
        cerr << "Importing GFA into XG" << endl;
        graph.from_gfa(args::get(gfa_in), args::get(validate),
                       args::get(base).empty() ? args::get(gfa_in) : args::get(base));
    } else {
        cerr << "Error: Input graph in either XG or GFA format is required" << endl;
        return 1;
    }
                                                    
    // load the snarls into a map
    map<handle_t, handle_t> snarl_map;
    vector<Snarl> top_level_snarls;
    if (!args::get(snarls_in).empty()) {
        std::ifstream in(args::get(snarls_in));
        // todo: not actually using the map -- can get rid of 
        vg::io::for_each<Snarl>(in, [&](Snarl& snarl) {
                handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
                handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
                snarl_map[start_handle] = end_handle;
                snarl_map[graph.flip(end_handle)] = snarl_map[graph.flip(start_handle)];
                if (!snarl.has_parent()) {
                    top_level_snarls.push_back(snarl);
                }
            });
        cerr << "Loaded " << top_level_snarls.size() << " top-level snarls" << endl;
    }

    if (!args::get(mat_out).empty()) {
        ofstream mat_file(args::get(mat_out));
        if (!mat_file) {
            cerr << "Error opening " << args::get(mat_out) << endl;
            return 1;
        }
        
        cerr << "Computing Path Coverage" << endl;
        pair<vector<string>, vector<vector<float>>> path_comp = all_paths_jaccard(graph);

        // print header row (columns are in same order)
        mat_file << ".";
        for (size_t i = 0; i < path_comp.first.size(); ++i) {
            mat_file << "\t" << path_comp.first[i];
        }
        mat_file << endl;

        // print the matrix
        for (size_t i = 0; i < path_comp.second.size(); ++i) {
            mat_file << path_comp.first[i];
            for (size_t j = 0; j < path_comp.second[i].size(); ++j) {
                mat_file << "\t" << path_comp.second[i][j];
            }
            mat_file << endl;
        }    
    }

    if (!args::get(pileup_out).empty()) {
        if (args::get(snarls_in).empty() || args::get(ref_path).empty()) {
            cerr << "Error: ref-path and snarls file required for pileup" << endl;
            return 1;
        }
        
        ofstream pu_file(args::get(pileup_out));
        if (!pu_file) {
            cerr << "Error opening " << args::get(pileup_out) << endl;
            return 1;
        }

        if (!graph.has_path(args::get(ref_path))) {
            cerr << "Error: " << args::get(ref_path) << " is not in the graph" << endl;
            return 1;
        }

        cerr << "computing pileup" << endl;
        vector<float> pileup = sv_pileup(graph, top_level_snarls, args::get(ref_path), window, min_sv_len);

        pu_file << "fixedStep" << "\t" << "chrom=" << args::get(ref_path) << "\t" << "step=" << window << endl;
        for (size_t i = 0; i < pileup.size(); ++i) {
            pu_file << pileup[i] << "\n";
        }
        pu_file << endl;
    }

    if (!args::get(coverage_stats).empty()) {
        print_path_coverage(graph, args::get(coverage_stats)); 
    }
    
    return 0;
}

