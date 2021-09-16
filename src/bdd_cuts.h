//
// Created by Margarita Castro on 2019-08-16.
//

#ifndef DD_CONIC_LIFTING_BDD_CUTS_H
#define DD_CONIC_LIFTING_BDD_CUTS_H

#include <iostream>
#include <ctime>			//measure time

#include <ilcplex/ilocplex.h>
#include "bdd_conic_diag.h"
#include "bdd_conic_general.h"
#include "bdd_lp.h"
#include "read_data.h"
#include "cut_generators.h"

using namespace std;

//Class to generate and lift class w.r.t. BDDs
class bdd_cuts {

public:

    // Constructor for diagonal chance constraints
    bdd_cuts(vector<CQDConstraint*>& _cqd, int _bdd_width, int _bdd_cut_type):
        num_diag(_cqd.size()), num_chance(0), bdd_width(_bdd_width), type_diag(true), bdd_cut_type(_bdd_cut_type) {

        debug = true;
        stop_weak_cuts = false;

        cq_diag.resize(num_diag, nullptr);
        bdds_diag.resize(num_diag, nullptr);
        bdds_lp.resize(num_diag, nullptr);
        slack.resize(num_diag, nullptr);
        pi.resize(num_diag, nullptr);
        x_value.resize(num_diag, nullptr);

        pi_t = nullptr;
        x_t = nullptr;
        x_best = nullptr;
        x_diff = nullptr;

        constraint_margin.resize(num_diag, make_pair(0,0) );

        for (int i = 0; i < num_diag; i++) cq_diag[i] = _cqd[i];

        bdd_time = 0;
        sep_time = 0;
        sep_time_avg = 0;
        sep_it = 0;
        lift_time = 0;
        lift_time_avg = 0;
        lift_it = 0;
    }

    // For general chance constriants
    bdd_cuts(vector<CQChanceConstraint*>& _cqc, int _bdd_width, int _bdd_cut_type):
            num_diag(0), num_chance(_cqc.size()), bdd_width(_bdd_width), type_diag(false), bdd_cut_type(_bdd_cut_type) {

        debug = true;
        stop_weak_cuts = false;

        cq_chance.resize(num_chance, nullptr);
        bdds_chance.resize(num_chance, nullptr);
        bdds_lp.resize(num_chance, nullptr);
        slack.resize(num_chance, nullptr);
        pi.resize(num_chance, nullptr);
        x_value.resize(num_chance, nullptr);

        pi_t = nullptr;
        x_t = nullptr;
        x_best = nullptr;
        x_diff = nullptr;

        constraint_margin.resize(num_chance, make_pair(0,0) );

        for (int i = 0; i < num_chance; i++) cq_chance[i] = _cqc[i];

        bdd_time = 0;
        sep_time = 0;
        sep_time_avg = 0;
        sep_it = 0;
        lift_time = 0;
        lift_time_avg = 0;
        lift_it = 0;
    }

    ~bdd_cuts(){
        for (int i =0; i < num_diag; i++) {
            if (bdds_diag[i] == nullptr) continue;

            delete bdds_diag[i];
            delete slack[i];
            delete pi[i];
            delete x_value[i];
            delete bdds_lp[i];
        }

        for (int i =0; i < num_chance; i++) {
            if (bdds_chance[i] == nullptr) continue;

            delete bdds_chance[i];
            delete slack[i];
            delete pi[i];
            delete x_value[i];
            delete bdds_lp[i];
        }

        if (bdd_cut_type == 3) {
            delete  pi_t;
            delete  x_t;
            delete  x_best;
            delete  x_diff;
        }

    }


    // Variables ======================================
private:
    //Debug mode
    bool debug;

    //Parameters
    const int num_diag;
    const int num_chance;
    const int bdd_width;
    const bool type_diag;
    int bdd_cut_type;

    bool stop_weak_cuts;

    // BDDs
    vector<bdd_diag*> bdds_diag;
    vector<bdd_chance*> bdds_chance;

    // BDD LP models
    vector<bdd_lp*> bdds_lp;

    //Slack and x_values vectors
    vector<double*> slack;
    vector<double*> x_value;

    // Constructed BDDs
    list<int> bdd_constructed;

    // To sort margin of constraints
    vector< pair<double, int> > constraint_margin;

    //BDD time
    double bdd_time;

    // BDD separation time
    double sep_time;        // Accummulated time for a single BDD to separate a point
    double sep_time_avg;    // Average time
    int sep_it;             // Total # of separation iterations

    // BDD lifting time
    double lift_time;
    double lift_time_avg;
    int lift_it;

    clock_t start_time;
    clock_t end_time;

    //For projected sub-gradient cuts
    double*  pi_t;
    int*  x_t;
    int*  x_best;
    double* x_diff;


public:
    vector<CQDConstraint*> cq_diag;
    vector<CQChanceConstraint*> cq_chance;
    vector<double*> pi;


    // Functions ======================================
public:

    int get_bdd_cut(IloNumArray& x_vals, double& pi0, bool multi_cuts, bool& finish, bool& equality);

    bool lift_bdd_constraint(int id, double& pi0);
    bool lift_cover_qcd(list<int>& cover, double& pi0, int id);

    void set_debug(bool _debug) {
        debug = _debug;
    }

    void set_weak_bdd_cuts() {
        bdd_cut_type = 1;
    };

    double get_bdd_time() {
        return bdd_time;
    }

    double get_separation_time() {
        return sep_time_avg;
    }

    double get_lifting_time() {
        return lift_time_avg;
    }

    int get_num_constraints() {
        if (type_diag) return num_diag;

        return num_chance;
    }


private:

    void create_bdd_diag(int id);
    void create_bdd_chance(int id);

    bool compute_weak_bdd_cut(int id, IloNumArray& x_vals, double& pi0);
    bool compute_cglp_bdd_cut(int id, IloNumArray& x_vals, double& pi0);
    bool compute_projected_bdd_cut(int id, IloNumArray& x_vals, double& pi0);
    bool compute_target_bdd_cut(int id, IloNumArray& x_vals, double& pi0, bool& equality);

    bool projected_gradient_routine(int id, int num_vars, double& pi0);

    int choose_slack(int id);

    void get_sorted_bdd_list(list<int>& bdd_list, IloNumArray& x_vals);
    double get_margin_diagonal(int id, IloNumArray& x_vals);
    double get_margin_chance(int id, IloNumArray& x_vals);

    void print_inequality(int id, double&pi0);

};

// ============================================
// Call-back class BDD cover cuts
// ============================================

class BDDCoverLiftI : public IloCplex::UserCutCallbackI {
    IloIntVarArray x;
    bdd_cuts* bdd_lift;
    cover_cuts* cover;
    IloNum& lift_count;

public:
    IloCplex::CallbackI* duplicateCallback() const {
        return (new(getEnv()) BDDCoverLiftI(*this));
    }

    BDDCoverLiftI(IloEnv env, IloIntVarArray _x, bdd_cuts* _bdd_lift, cover_cuts* _cover,  IloNum& _lift_count) :
            IloCplex::UserCutCallbackI(env), x(_x), bdd_lift(_bdd_lift), cover(_cover),  lift_count(_lift_count) {

        rhs = 0;
        num_vars = x.getSize();
    }
    void main();

private:
    // Auxiliary variables for computing cuts
    IloInt num_vars;
    double rhs;
    int id_con;
};

IloCplex::Callback
BDDCoverLift(IloEnv env, IloIntVarArray _x, bdd_cuts* _bdd_lift,  cover_cuts* _cover, IloNum& _lift_count);

// ============================================
// Call-back class BDD flow cut and lift
// ============================================

class BDDCutLiftI : public IloCplex::UserCutCallbackI {
    IloIntVarArray x;
    bdd_cuts* bdd;
    IloNum& lift_count;


public:
    IloCplex::CallbackI* duplicateCallback() const {
        return (new(getEnv()) BDDCutLiftI(*this));
    }

    BDDCutLiftI(IloEnv env, IloIntVarArray _x, bdd_cuts* _bdd,  IloNum& _lift_count, bool _type_diag, bool _lift, cplex_param* param) :
            IloCplex::UserCutCallbackI(env), x(_x), bdd(_bdd), lift_count(_lift_count), type_diag(_type_diag), lift(_lift),
            multi_cut(param->multi_cut), root_cuts(param->root_cuts), weak_cut_tree(param->weak_tree_cuts)  {

        rhs = 0;
        num_vars = x.getSize();
        cuts_count =0;
    }
    void main();
private:
    // Auxiliary variables for computing cuts
    IloInt num_vars;
    double rhs;
    int id_con;
    const bool type_diag;
    const bool lift;
    const bool multi_cut;
    const bool root_cuts;
    const bool weak_cut_tree;
    int cuts_count;
};

IloCplex::Callback
BDDCutLift(IloEnv env, IloIntVarArray _x, bdd_cuts* _bdd_lift, IloNum& _lift_count, bool _type_diag, bool _lift, cplex_param* param);

#endif //DD_CONIC_LIFTING_BDD_CUTS_H