//
// Created by Margarita Castro on 2019-08-05.
//

#ifndef DD_CONIC_LIFTING_CUT_GENERATORS_H
#define DD_CONIC_LIFTING_CUT_GENERATORS_H


#include "read_data.h"
#include <list>
#include <math.h>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <ctime>			//measure time

using namespace std;

class cover_cuts {

public:
    cover_cuts(vector<CQDConstraint*>& _cqd):
    cqd(_cqd), num_con(_cqd.size()) {

        debug = false;
        is_lift= false;

        weights.resize(num_con, vector<pair <double,int> >() );
        coeff.resize(num_con, vector<double>());
        is_cqp_created.resize(num_con, false);

        for (int i = 0; i < num_con; i++) {
            weights[i].resize(cqd[i]->vars.size(), make_pair(-1,-1));
            coeff[i].resize(cqd[i]->vars.size(), 0);
        }

        //CPLEX objects
        env.resize(num_con, nullptr);
        cplex.resize(num_con, nullptr);
        obj.resize(num_con, nullptr);
        x.resize(num_con, nullptr);

        //Set time variables
        sep_time = 0;
        sep_time_avg = 0;
        sep_it = 0;

        lift_time = 0;
        lift_time_avg = 0;
        lift_it = 0;

    }

    ~cover_cuts(){
        for (int i = 0; i < num_con; i++) {
            if(is_cqp_created[i]) env[i]->end();

            delete env[i];
            delete cplex[i];
            delete obj[i];
            delete x[i];
        }

    }

    // Variables ======================================
public:
    vector<CQDConstraint*> cqd;
    vector< vector<double> > coeff;

    bool debug;
    bool is_lift;

private:
    const int num_con;

    vector<bool> is_cqp_created;
    vector< vector< pair <double,int> > > weights;

    // Separation time
    double sep_time;        // Accummulated time for a single BDD to separate a point
    double sep_time_avg;    // Average time
    int sep_it;             // Total # of separation iterations

    // Lifting time
    double lift_time;
    double lift_time_avg;
    int lift_it;

    //Clocks
    clock_t start_time;
    clock_t end_time;

    // CQP solver
    vector<IloEnv*> env;
    vector<IloCplex*> cplex;
    vector<IloObjective*> obj;
    vector<IloNumVarArray*> x;

    // Functions ======================================
public:

    int get_cover_cut(IloNumArray& x,  list<int>& cover_set);
    void lift_cover_inequality(int id, list<int>& cover_set);

    double get_separation_time() {
        return sep_time_avg;
    }

    double get_lifting_time() {
        return lift_time_avg;
    }

private:

    // Methods to compute cover cut as in Atamturk (2009)
    bool find_cover_set(IloNumArray& x, int id, list<int>& cover_set);
    void compute_sorted_weights(IloNumArray& x, int id,  int i, int j);
    bool is_cover(int id,  list<int>& cover_set);
    bool cuts_solution(IloNumArray& x,int id,  list<int>& cover_set);

    // Methods to compute lifted cover cut as in Atamturk (2009)
    void fix_coefficients(int id, list<int>& cover_set, list<int>& remaining_vars);
    void get_mu_nu(int id, list<int>& S, double& mu, double& nu);

    //Methods to solve and create CQP for lifting
    void create_cqp(int id);
    void update_run_cqp(int id, int var, int cover_size);

    //Auxiliary method to evaluate function
    void evaluate_constraint(int id,  list<int>& S, double& value);
    void diff_constraint(int id,  list<int>& S, int i, double& diff);


};



// ===========================
// Call-back class implementation for cover cuts
// ===========================
class CoverCutsI : public IloCplex::UserCutCallbackI {
    IloIntVarArray x;
    cover_cuts*  cover;

public:
    IloCplex::CallbackI* duplicateCallback() const {
        return (new(getEnv()) CoverCutsI(*this));
    }

    CoverCutsI(IloEnv env, IloIntVarArray _x, cover_cuts* _cover) :
            IloCplex::UserCutCallbackI(env), x(_x), cover(_cover) {
    }
    void main();

};

IloCplex::Callback
CoverCuts(IloEnv env, IloIntVarArray _x, cover_cuts* _cover);


// ===========================
// Call-back class implementation for lifted cover cuts
// ===========================

class LiftCoverCutsI : public IloCplex::UserCutCallbackI {
    IloIntVarArray x;
    cover_cuts*  cover;
    IloNum& lift_count;

public:
    IloCplex::CallbackI* duplicateCallback() const {
        return (new(getEnv()) LiftCoverCutsI(*this));
    }

    LiftCoverCutsI(IloEnv env, IloIntVarArray _x, cover_cuts* _cover, IloNum&  _lift_count) :
            IloCplex::UserCutCallbackI(env), x(_x), cover(_cover),lift_count(_lift_count)  {
    }
    void main();

};

IloCplex::Callback
LiftCoverCuts(IloEnv env, IloIntVarArray _x, cover_cuts* _cover, IloNum&  _lift_count);

#endif //DD_CONIC_LIFTING_CUT_GENERATORS_H
