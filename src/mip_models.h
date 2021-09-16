//
// Created by Margarita Castro on 2019-08-05.
//

#ifndef DD_CONIC_LIFTING_MIP_MODELS_H
#define DD_CONIC_LIFTING_MIP_MODELS_H


//#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>
#include "read_data.h"


//MIP models
void mip_qc_diagonal(char* filename, string outputname, int cover_cut, bdd_param* bp, cplex_param* cp);
void mip_qc_chance(char* filename, string outputname, bdd_param* bp, cplex_param* cp);


//Check feasibility of solutions
void check_solution_qc_diag(IloCplex& cplex, IloIntVarArray& x, vector<CQDConstraint*>& CQDConstraints, int& num_vars);
void check_solution_qc_chance(IloCplex& cplex, IloIntVarArray& x, IloNumVarArray& z, vector<CQChanceConstraint*>& cqc_constrs, int& num_vars);
void check_objective_function(IloCplex& cplex, IloIntVarArray& x,vector<int>& cost);

//CPLEX cuts
void desactivate_cplex_cuts(IloCplex& cplex);
void get_cut_count(IloCplex& cplex, int& n_cuts);

//Call-back class to get information on the root node/first node
class RootNodeInfoI : public IloCplex::MIPInfoCallbackI {
    IloNum& time;
    IloNum& gap;
    IloNum& bound;

public:
    IloCplex::CallbackI* duplicateCallback() const {
        return (new(getEnv()) RootNodeInfoI(*this));
    }

    RootNodeInfoI(IloEnv env, IloNum& _time, IloNum& _gap, IloNum& _bound) :
            IloCplex::MIPInfoCallbackI(env), time(_time), gap(_gap), bound(_bound) {}
    void main();

};

IloCplex::Callback
RootNodeInfo(IloEnv env, IloNum &_time, IloNum &_gap, IloNum& _bound);

#endif //DD_CONIC_LIFTING_MIP_MODELS_H
