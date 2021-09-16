//
// Created by Margarita Castro on 2019-10-01.
//

#ifndef DD_CONIC_LIFTING_BDD_CONIC_GENERAL_H
#define DD_CONIC_LIFTING_BDD_CONIC_GENERAL_H



#include "bdd_conic_general.h"
#include "read_data.h"
#include <iostream>
#include <list>
#include <utility>
#include <tuple>
#include <ctime>			//measure time

using namespace std;

// BDD Graph for SOCs with diagonal matrix

class bdd_chance {

public:
    bdd_chance(CQChanceConstraint* _cqc, int _id, int _width):
            cqc(_cqc), bdd_id(_id), num_vars(_cqc->vars.size()), max_width(_width) {
        is_constructed = false;
        bdd_has_change = false;

        debug = true;

        cost = nullptr;
        slack = nullptr;

        epsilon = 0.001;
        int_min = INT32_MIN/2;
        int_max = INT32_MAX/2;

        overflow = false;
    }

    ~bdd_chance();

    //Structures
    struct bdd_node;
    struct bdd_arc;


    // Variables ======================================
private:
    //Debug mode
    bool debug;

    //Precision
    double epsilon;
    int int_min;
    int int_max;

    //Conic Constraint
    const CQChanceConstraint* cqc;
    const int bdd_id;

    //BDD size variables
    const int num_vars;
    const int max_width;

    //Cost & slack variables
    double* cost;
    double* slack;

    vector<double> cost_edge;

    //BDD status variables
    bool is_constructed;
    bool bdd_has_change;

    //Graph structure
    vector< vector <bdd_arc*> > arcs;
    vector< vector <bdd_node*> > nodes;

    // Nods availability
    vector<int> nodes_available;
    vector<int> first_node_available;

    vector< tuple<double, int, int> > candidate_arcs;
    vector<int> candidate_nodes;

    //Active nodes id (i.e., non-free)
    vector< vector<int> > active_nodes;

    //Exact layer
    vector<bool> is_layer_exact;

    //Measure time
    double constr_time;

    //Center point variables
    vector<double> center_point;
    vector< vector<long long> > count_paths_top;
    vector< vector<long long> > count_paths_bottom;

    bool overflow;

    // Functions ======================================
public:
    void construct_bdd();

    bool separate_min_cut(double* _cost, double* pi, double& pi0);

    bool get_bdd_slacks(double* _cost, double* _slacks, double& rhs);

    void compute_longest_path(double* _cost, int* path);

    void export_graph(const char* filename);

    double get_constr_time() {
        return constr_time;
    }

    void set_debug(bool _debug) {
        debug = _debug;
    }

    int get_num_vars() {
        return num_vars;
    }

    const vector< vector <bdd_node*> >& get_nodes() {
        return nodes;
    }

    const vector< vector <bdd_arc*> >& get_arcs() {
        return arcs;
    }

    vector<double>& get_center_point(bool& _overflow) {
        if (center_point.empty() ) compute_center_point();
        _overflow = overflow;

        return  center_point;
    }

private:

    void initialize_bdd_space();
    void width_one_bdd();

    void procedure_top_down(int round);
    void procedure_bottom_up(int round);

    void update_node_top_down(int layer, int node_id);
    void update_node_bottom_up(int layer, int node_id);

    bool is_arc_infeasible(bdd_arc* arc);
    bool is_arc_infeasible_qcc(bdd_arc *arc);

    bool is_arc_feasible(bdd_arc* arc);

    void compute_min_sigma(bdd_arc* arc, int& sigma);
    void compute_max_sigma(bdd_arc* arc, int& sigma);

    void refine(int layer);

    void reduce_bdd();
    void reduce_layer(int layer, bool final_bdd);
    void reduce_nodes(int layer, int node_id, list<int>& similar);
    int get_nodes_to_reduce(int layer);
    int get_nodes_to_reduce_final(int layer);

    void get_candidate_arcs(int layer, int& total);
    bool are_arcs_equal(bdd_arc* arc1, bdd_arc* arc2);
    bool are_in_arcs_min_equal(int layer, int id);

    void check_exact_layer(int layer);

    void create_new_node(int layer, int old_id, list<int>& in_arcs);

    void update_first_node_available(int layer);

    void remove_arc(bdd_arc* arc);
    void remove_node(int layer, int node_id);

    // Min-cut separation functions
    void set_values_max_flow();
    void compute_max_flow(double& max_flow);
    bool find_augmenting_path();
    void reset_visited_nodes();

    void get_min_cut_arcs(double* pi, double& pi0);
    void compute_visited_nodes();

    bool is_constraint_violated(double* pi, double& pi0);

    void get_active_nodes();

    // Slack functions
    void compute_cost_top();
    void compute_cost_bottom();
    void compute_slacks();

    void compute_cost_node_top(int layer, int id);
    void compute_cost_node_bottom(int layer, int id);

    //Longest path functions
    void extract_longest_path(int* path);

    //Center BDD point function (Needed for target cuts)
    void compute_center_point();
    void count_path_top();
    void count_path_bottom();

};

//---------------------------------- BDD Node Structure --------------------------------------//

struct bdd_chance::bdd_node {

    // Top_down information
    int lin_min_top;
    int lin_max_top;

    vector<int> sigma_min_top;
    vector<int> sigma_max_top;

    // Bottom_up information
    int lin_min_bottom;
    int lin_max_bottom;

    vector<int> sigma_min_bottom;
    vector<int> sigma_max_bottom;

    // Node status
    bool is_feasible;       //all path passing throught the node are feasible
    bool is_free;
    bool is_exact;
    bool is_reduce;

    // Incoming arcs;
    list<int> in_arcs;

    // Cost computation variables
    double cost_top;
    double cost_bottom;
    int cost_arc_id;

    // Max-flow parameters
    bool visited;
    pair<int,int> parent; // (layer, node)

    bdd_node(const int& num_vars) {
        is_feasible = false;
        is_free = true;
        is_exact = false;
        is_reduce = false;

        lin_min_top = 0;
        lin_max_top = 0;
        sigma_min_top.resize(num_vars, 0);
        sigma_max_top.resize(num_vars, 0);

        lin_min_bottom = INT32_MIN/2;   //minimum value it can have (can be negative)
        lin_max_bottom = INT32_MAX/2;
        sigma_min_bottom.resize(num_vars, INT32_MIN/2);  //minimum value it can have (can be negative)
        sigma_max_bottom.resize(num_vars, INT32_MAX/2);

        cost_top = 0;
        cost_bottom = 0;
        cost_arc_id = -1;

        visited = false;
        parent = make_pair(-1,-1);
    }

};

//---------------------------------- BDD Arc Structure --------------------------------------//
struct bdd_chance::bdd_arc {

    // Layer
    const int layer;
    // Source node
    const int source;
    // Target node
    int target;
    // Value of arc (label: 0 or 1)
    const int value;

    // Max-flow values
    double res_value;
    double res_rev_value;

    //Constructor
    bdd_arc(int _layer, int _source, int _value)
            : layer(_layer), source(_source), target(-1), value(_value)
    { }
};



#endif //DD_CONIC_LIFTING_BDD_CONIC_GENERAL_H
