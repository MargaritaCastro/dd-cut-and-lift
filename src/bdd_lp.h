//
// Created by Margarita Castro on 2019-12-18.
//

#ifndef DD_CONIC_LIFTING_BDD_LP_H
#define DD_CONIC_LIFTING_BDD_LP_H

#include <iostream>
#include <ilcplex/ilocplex.h>
#include "bdd_conic_diag.h"
#include "bdd_conic_general.h"


using namespace std;
typedef IloArray<IloNumVarArray> IloNumVarArray2;

class bdd_lp {

public:

    // Constructor for any BDD structure
    bdd_lp(int _id): bdd_id(_id) {

        IloEnv env;

        model = IloModel(env);
        cplex = IloCplex(env);
        obj = IloObjective(env);

        cplex.extract(model);

        constructred = false;
        debug = false;
        overflow = false;
    }

    ~bdd_lp(){
        model.getEnv().end();
    }

    // Variables
private:

    const int bdd_id;

    //CPLEX
    IloModel model;
    IloCplex cplex;

    //Model elements
    IloObjective obj;
    IloNumVarArray mu;
    IloNumVarArray lambda;
    IloNumVarArray2 w;
    IloNum num_vars;

    //Status variables
    bool constructred;
    bool debug;
    bool overflow;  // for target cuts

    //Center point
    vector<double> center;

    // Variable ID (for rays - target cuts)
    vector<int> mu_id;

    //Functions
public:

    // LP for our general max-flow cuts
    void create_lp_model(bdd_chance* bdd);
    void create_lp_model(bdd_diag* bdd);
    bool solve_lp(double* xval, double* pi, double& pi0);

    // LP for target cuts
    void create_lp_model_target_cuts(bdd_chance* bdd, vector<double>& _center);
    void create_lp_model_target_cuts(bdd_diag* bdd, vector<double>& _center);

    bool solve_lp_target_cuts(double* xval, double* pi, double& pi0, bool& equality);

    void set_debug(bool _debug) {
        debug = _debug;
    }

    void set_overflow() {
        overflow = true;
    }

    bool get_overflow() {
        return  overflow;
    }

private:
    void update_objective(double* xval);
    void update_objective_target_cuts(double* xval);

};
#endif //DD_CONIC_LIFTING_BDD_LP_H
