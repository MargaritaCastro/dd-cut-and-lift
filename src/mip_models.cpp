//
// Created by Margarita Castro on 2019-08-05.
//

#include <list>
#include "mip_models.h"
#include "bdd_cuts.h"
#include "cut_generators.h"

#define VARNAME_SIZE 256

using namespace std;

// ========================================================
// MIP MODEL - CONIC QUADRATIC
// ========================================================

void mip_qc_diagonal(char* filename, string outputname, int cover_cut, bdd_param* bp, cplex_param* cp) {

    //Read instance
    int num_vars = 0;
    int num_cons = 0;
    vector<int> cost;
    vector<CQDConstraint*> CQDConstraints;
    bdd_cuts* bdd_lift = nullptr;
    cover_cuts * covers = nullptr;

    read_conic_diagonal_input(filename, CQDConstraints, num_vars, cost);

    //Create MIP model
    IloEnv env;
    try {
        char varname[VARNAME_SIZE];

        IloModel model(env);
        IloCplex cplex(env);

        IloNum start_time = cplex.getCplexTime();

        // Create Variables
        IloIntVarArray x(env,num_vars,0,1);

        for (int i =0; i <num_vars; i++) {
            sprintf(varname, "x_%d", i);
            x[i].setName(varname);
            model.add(x[i] >= 0);
        }

        //Create Conic-Quadratic constraints
        num_cons = CQDConstraints.size();
        int num_qvars = 0;

        for (int i = 0; i < num_cons; i ++) {
            num_qvars = CQDConstraints[i]->vars.size();

            //Create auxiliary variables
            IloNumVarArray y(env,num_qvars, 0, IloInfinity);

            for (int j = 0; j <num_qvars; j++) {
                sprintf(varname, "y_%d_%d", i, CQDConstraints[i]->vars[j] );
                y[j].setName(varname);
            }

            sprintf(varname, "z_%d", i);
            IloNumVar z(env, 0, IloInfinity, ILOFLOAT, varname);

            //Add quadratic term
            model.add(IloScalProd(y,y) <= z*z);

            //Relate auxiliary variables to normal variables
            for (int j =0; j < num_qvars; j++) {
                model.add(y[j] == CQDConstraints[i]->diag[j] * x[CQDConstraints[i]->vars[j]]);
            }

            // Add linear term
            IloExpr qConstraint(env);
            for (int j =0; j < num_qvars; j++){
                qConstraint += CQDConstraints[i]->linear[j] * x[CQDConstraints[i]->vars[j]];
            }

            qConstraint += CQDConstraints[i]->omega * z;
            model.add(qConstraint <= CQDConstraints[i]->rhs );
        }

        //Objective function
        IloExpr obj(env);
        for (IloInt i = 0; i < num_vars; ++i) {
            obj += x[i]*cost[i];
        }
        model.add(IloMaximize(env, obj));
        obj.end();

        //Solver parameters
        cplex.extract(model);

        cplex.setParam(IloCplex::Param::MIP::Strategy::MIQCPStrat, cp->strategy);
        cplex.setParam(IloCplex::Param::Threads, 1);
        //cplex.setParam(IloCplex::Param::TimeLimit, 600); // 10 minutes time limit
        cplex.setParam(IloCplex::Param::TimeLimit, 3600); // 1 hour time limit
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
        cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);

        //cplex.exportModel("test.lp");

        IloNum root_time = 0;
        IloNum root_gap = 0;
        IloNum root_bound = 0;
        IloNum lift_count = 0;

        cplex.use(RootNodeInfo(env, root_time, root_gap, root_bound));

        if (bp->width > 0) {
            bdd_lift = new bdd_cuts(CQDConstraints, bp->width, bp->cut_type);
            bdd_lift->set_debug(false);
        }
        if (cover_cut > 0) {
            covers = new cover_cuts(CQDConstraints);
            covers->debug = false;
        }

        // De-activate CPLEX cuts
        if (!cp->cplex_cuts) desactivate_cplex_cuts(cplex);

        // Cover cuts
        if (bp->width == 0 && cover_cut > 0) {
            if (cover_cut == 1)
                cplex.use(CoverCuts(env, x, covers));
            else
                cplex.use(LiftCoverCuts(env, x, covers, lift_count));
        }
        if (bp->width > 0) {
            if (cover_cut == 1)
                cplex.use(BDDCoverLift(env, x, bdd_lift, covers, lift_count));
            else
                cplex.use(BDDCutLift(env, x, bdd_lift, lift_count, true, bp->lift, cp));
        }


        //Solve problem
        cplex.solve();

        IloNum end_time = cplex.getCplexTime();

        // Get total cuts
        int n_cuts = 0;
        int n_cuts_user = cplex.getNcuts(IloCplex::CutType::CutUser) ;
        get_cut_count(cplex, n_cuts);

        // Get BDD time
        double bdd_time = 0;
        double sep_time = 0;
        double lift_time = 0;

        if (bdd_lift != nullptr) {
            bdd_time = bdd_lift->get_bdd_time();
            sep_time = bdd_lift->get_separation_time();
            lift_time = bdd_lift->get_lifting_time();

            if (cover_cut > 0) sep_time = covers->get_separation_time();

        }
        else if (cover_cut > 0) {
            sep_time = covers->get_separation_time();
            lift_time = covers->get_lifting_time();
        }

        //Print statistics
        cout << "Strategy:\t" << cp->strategy << endl;
        cout << "Status:\t" << cplex.getStatus() << endl;
        cout << "Objective:\t" << cplex.getObjValue(-1) << endl;
        cout << "Gap:\t" << cplex.getMIPRelativeGap() << endl;
        cout << "Gap root:\t" << root_gap << endl;
        cout << "UB root:\t" << root_bound << endl;
        cout << "Time:\t" << end_time - start_time << endl;
        cout << "Time root:\t" << root_time - start_time << endl;
        cout << "Time BDD:\t" << bdd_time << endl;
        cout << "Nodes:\t" << cplex.getNnodes() << endl;
        cout << "CPLEX Cuts:\t" << n_cuts <<endl;
        cout << "User Cuts:\t" << n_cuts_user << endl;
        cout << "Lifted Cuts:\t" << lift_count << endl;

        // Write results
        if (!outputname.empty()) {
            fstream output;

            output.open(outputname.c_str(), fstream::in | fstream::out | fstream::app);
            output << filename << "\t"
                    << num_vars << "\t"
                    << num_cons << "\t"
                    << CQDConstraints[0]->omega << "\t"
                    << cp->strategy << "\t"
                    << cover_cut << "\t"
                    << bp->lift << "\t"
                    << bp->cut_type                 << "\t"
                    << cp->cplex_cuts               << "\t"          // Signal for using CPLEX cuts or not
                    << cp->multi_cut                << "\t"
                    << bp->width                   << "\t"
                    << cplex.getNnodes()            << "\t"
                    << cplex.getStatus()		    << "\t"
                    << cplex.getObjValue()          << "\t"
                    << cplex.getBestObjValue()      << "\t"
                    << root_time - start_time       << "\t"         // Root node time
                    << bdd_time                     << "\t"         // total BDD time
                    << end_time - start_time        << "\t"         // Total time
                    << root_gap                     << "\t"         // Root gap
                    << root_bound                   << "\t"         // UB in the root node
                    << cplex.getMIPRelativeGap()    << "\t"         // Finish gap
                    << n_cuts                       << "\t"
                    << n_cuts_user                  << "\t"
                    << lift_count                   << "\t"
                    << sep_time                     << "\t"
                    << lift_time                    << endl;

            output.close();
        }

        //Check validity of the solution
        check_solution_qc_diag(cplex, x, CQDConstraints, num_vars);
        check_objective_function(cplex, x, cost);

    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }

    //Delete objects
    delete bdd_lift;
    delete covers;

    for (int i =0; i < num_cons; i++){
        delete CQDConstraints[i];
    }

    env.end();
}

// ========================================================
// MIP MODEL - CHANCE CONSTRAINT
// ========================================================

void mip_qc_chance(char* filename, string outputname, bdd_param* bp, cplex_param* cp){

    //Read instance
    int num_vars = 0;
    int num_cons = 0;
    vector<int> cost;
    vector<CQChanceConstraint*> cqc_constrs;
    bdd_cuts* bdd_lift = nullptr;

    read_conic_chance_input(filename, cqc_constrs, num_vars, cost);

    //Create MIP model
    IloEnv env;

    try {
        char varname[VARNAME_SIZE];

        IloModel model(env);
        IloCplex cplex(env);

        IloNum start_time = cplex.getCplexTime();

        // Create Variables
        IloIntVarArray x(env,num_vars,0,1);

        for (int i =0; i <num_vars; i++) {
            sprintf(varname, "x_%d", i);
            x[i].setName(varname);
            model.add(x[i] >= 0);
        }

        //Create Conic-Quadratic constraints
        num_cons = cqc_constrs.size();
        int num_qvars = 0;

        IloNumVarArray z(env, num_cons, 0, IloInfinity, ILOFLOAT);

        for (int i = 0; i < num_cons; i ++) {
            num_qvars = cqc_constrs[i]->vars.size();

            //Create auxiliary variables
            IloNumVarArray y(env,num_qvars, 0, IloInfinity);

            for (int j = 0; j <num_qvars; j++) {
                sprintf(varname, "y_%d_%d", i, cqc_constrs[i]->vars[j] );
                y[j].setName(varname);
            }

            sprintf(varname, "z_%d", i);
            z[i].setName(varname);
            model.add(z[i] >= 0);

            //Add quadratic term
            model.add(IloScalProd(y,y) <= z[i]*z[i]);

            //Relate auxiliary variables to normal variables

            for (int j =0; j < num_qvars; j++) {
                IloExpr std(env);

                for (int k =0; k < num_qvars; k++) {
                    std += cqc_constrs[i]->sigma[j][k] * x[cqc_constrs[i]->vars[k]];
                }

                model.add( std <= y[j]);
                model.add( -std <= y[j]);
            }

            // Add linear term
            IloExpr qConstraint(env);
            for (int j =0; j < num_qvars; j++){
                qConstraint += cqc_constrs[i]->linear[j] * x[cqc_constrs[i]->vars[j]];
            }

            qConstraint += cqc_constrs[i]->omega * z[i];
            model.add(qConstraint <= cqc_constrs[i]->rhs );

        }

        //Objective function
        IloExpr obj(env);
        for (IloInt i = 0; i < num_vars; ++i) {
            obj += x[i]*cost[i];
        }
        model.add(IloMaximize(env, obj));
        obj.end();

        //Solver parameters
        cplex.extract(model);
        //cplex.exportModel("chanceConstraints.lp");

        //Relaxation: quadratic or linearization
        cplex.setParam(IloCplex::Param::MIP::Strategy::MIQCPStrat, cp->strategy);

        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::TimeLimit, 3600); // 1 hour time limit
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
        cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);

        IloNum root_time = 0;
        IloNum root_gap = 0;
        IloNum root_bound = 0;
        IloNum lift_count = 0;

        cplex.use(RootNodeInfo(env, root_time, root_gap, root_bound));

        // De-activate CPLEX cuts
        if (!cp->cplex_cuts) desactivate_cplex_cuts(cplex);

        if (bp->width > 0 && bp->cut_type >= 1 ) {
            bdd_lift = new bdd_cuts(cqc_constrs, bp->width, bp->cut_type);
            bdd_lift->set_debug(false);

            cplex.use(BDDCutLift(env, x, bdd_lift, lift_count, false, bp->lift, cp));
        }

        //Solve problem
        cplex.solve();

        IloNum end_time = cplex.getCplexTime();

        // Get total cuts
        int n_cuts = 0;
        int n_cuts_user = cplex.getNcuts(IloCplex::CutType::CutUser) ;
        get_cut_count(cplex, n_cuts);

        // Get BDD time
        double bdd_time = 0;
        double sep_time = 0;
        double lift_time = 0;

        if (bdd_lift != nullptr) {
            bdd_time = bdd_lift->get_bdd_time();
            sep_time = bdd_lift->get_separation_time();
            lift_time = bdd_lift->get_lifting_time();
        }

        //Print statistics
        cout << "Strategy:\t" << cp->strategy << endl;
        cout << "Status:\t" << cplex.getStatus() << endl;
        cout << "Objective:\t" << cplex.getObjValue(-1) << endl;
        cout << "Gap:\t" << cplex.getMIPRelativeGap() << endl;
        cout << "Gap root:\t" << root_gap << endl;
        cout << "UB root:\t" << root_bound << endl;
        cout << "Time:\t" << end_time - start_time << endl;
        cout << "Time root:\t" << root_time - start_time << endl;
        cout << "Time BDD:\t" << bdd_time << endl;
        cout << "BDD width:\t" << bp->width << endl;
        cout << "BDD lift:\t" << bp->lift << endl;
        cout << "BDD cut type:\t" << bp->cut_type << endl;
        cout << "Nodes:\t" << cplex.getNnodes() << endl;
        cout << "CPLEX Cuts:\t" << n_cuts <<endl;
        cout << "User Cuts:\t" << n_cuts_user << endl;
        cout << "Lifted Cuts:\t" << lift_count << endl;


        // Write results
        if (!outputname.empty()) {
            fstream output;

            output.open(outputname.c_str(), fstream::in | fstream::out | fstream::app);
            output << filename << "\t"
                    << num_vars << "\t"
                    << num_cons << "\t"
                    << cqc_constrs[0]->omega << "\t"
                    << cp->strategy << "\t"
                    << bp->lift << "\t"
                    << bp->width                   << "\t"
                    << cplex.getNnodes()            << "\t"
                    << cplex.getStatus()		    << "\t"
                    << cplex.getObjValue()          << "\t"
                    << cplex.getBestObjValue()      << "\t"
                    << root_time - start_time       << "\t"         // Root node time
                    << bdd_time                     << "\t"         // total BDD time
                    << end_time - start_time        << "\t"         // Total time
                    << root_gap                     << "\t"         // Root gap
                    << root_bound                   << "\t"         // UB in the root node
                    << cplex.getMIPRelativeGap()    << "\t"         // Finish gap
                    << n_cuts                       << "\t"
                    << n_cuts_user                  << "\t"
                    << lift_count                   << "\t"
                    << cp->cplex_cuts               << "\t"          // Signal for using CPLEX cuts or not
                    << cp->multi_cut                << "\t"
                    << bp->cut_type                 << "\t"
                    << sep_time                     << "\t"
                    << lift_time                    << endl;

            output.close();
        }

        //Check validity of the solution
        check_objective_function(cplex, x, cost);
        check_solution_qc_chance(cplex, x, z, cqc_constrs, num_vars);

    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }

    //Delete objects
    delete bdd_lift;

    for (int i =0; i < num_cons; i++){
        delete cqc_constrs[i];
    }

    env.end();
}

// ========================================================
// CALLBACK ROOT NODE INFO
// ========================================================

IloCplex::Callback
RootNodeInfo(IloEnv env, IloNum &_time, IloNum &_gap, IloNum& _bound) {
    return (IloCplex::Callback(new(env) RootNodeInfoI(env, _time, _gap, _bound)));
}

void RootNodeInfoI::main() {
    // Don't call this after the first node
    if (getNnodes() > 1)
        return;

    time = getCplexTime();
    gap = getMIPRelativeGap();
    bound = getBestObjValue();
}

// ========================================================
// CPLEX CUT PARAMETERS
// ========================================================

void desactivate_cplex_cuts(IloCplex& cplex) {

    cplex.setParam(IloCplex::Param::MIP::Cuts::BQP, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Disjunctive, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::GUBCovers	, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::Implied	, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj	, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::MCFCut	, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut	, -1);
    cplex.setParam(IloCplex::Param::MIP::Cuts::PathCut	, -1 );
    cplex.setParam(IloCplex::Param::MIP::Cuts::RLT	, -1 );
    cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
}

void get_cut_count(IloCplex& cplex, int& n_cuts) {
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutBenders);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutBQP);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutClique);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutCover);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutDisj);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutFlowCover);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutFlowPath);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutFrac);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutGubCover);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutImplBd);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutLiftProj);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutLocalCover);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutLocalImplBd);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutMCF);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutMir);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutObjDisj);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutRLT);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutSolnPool);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutTable);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutTighten);
    n_cuts += cplex.getNcuts(IloCplex::CutType::CutZeroHalf);
}

// ========================================================
// CHECK VALIDITY OF SOLUTIONS
// ========================================================


void check_solution_qc_diag(IloCplex& cplex, IloIntVarArray& x, vector<CQDConstraint*>& CQDConstraints, int& num_vars) {

    // Check Feasibility of QC_Diagonal Constraints
    int num_cons = CQDConstraints.size();
    double rhs_linear = 0.0;
    double rhs_diag = 0.0;
    int var_id = 0;
    int num_qvars = 0;

    //Get variables values
    vector<double> x_vals(num_vars);

    for (IloInt i = 0; i < num_vars; i++) {
        x_vals[i] = cplex.getValue(x[i]);
    }

    //Print Solution
    cout << "Selected x values :" << endl;
    for (IloInt i = 0; i < num_vars; i++) {
        if (x_vals[i] == 1)
            cout << "\tx[" << i << "]";
    }

//    cout << " z values :" << endl;
//    for (IloInt i = 0; i < num_cons; i++) {
//        cout << "\tz[" << i << "] =" << cplex.getValue(z[i]) << endl;
//    }
    cout << endl;

    // Check SOC constraints
    for (int i = 0; i < num_cons; i ++) {
        num_qvars = CQDConstraints[i]->vars.size();
        rhs_linear = 0;
        rhs_diag = 0;

        for (int j = 0; j <num_qvars; j++) {
            var_id = CQDConstraints[i]->vars[j];
            rhs_linear += x_vals[var_id] *  CQDConstraints[i]->linear[j];
            rhs_diag += x_vals[var_id] *  CQDConstraints[i]->diag2[j];
        }
        
        if (rhs_linear + sqrt(rhs_diag) <= CQDConstraints[i]->rhs) {
            cout << " - CQ Diagonal Constraint " << i << " is correct" << endl;
        }
        else {
            cout << " - Error CQ Diagonal Constraint " << i << endl;
            cout << "\t "<< rhs_linear << " + sqrt( " << rhs_diag << " ) = " << rhs_linear + sqrt(rhs_diag);
            cout << " > " << CQDConstraints[i]->rhs << endl;
            assert(false);
        }
    }

}

void check_solution_qc_chance(IloCplex& cplex, IloIntVarArray& x, IloNumVarArray& z, vector<CQChanceConstraint*>& cqc_constrs, int& num_vars){

    //Check Feasibility of QC_Diagonal Constraints
    int num_cons = cqc_constrs.size();
    double rhs_linear = 0.0;
    vector<double> rhs_std;
    double rhs_quad = 0;
    int var_id = 0;
    int num_qvars = 0;


    //Get variables values
    vector<double> x_vals(num_vars);

    for (IloInt i = 0; i < num_vars; i++) {
        x_vals[i] = cplex.getValue(x[i]);
    }

    //Print Solution
    cout << "Selected x values :" << endl;
    for (IloInt i = 0; i < num_vars; i++) {
        if (x_vals[i] >= 0.95)
            cout << "\tx[" << i << "]";
    }
    cout << endl;

    cout << " z values :" << endl;
    for (IloInt i = 0; i < num_cons; i++) {
        cout << "\tz[" << i << "] =" << cplex.getValue(z[i]) << endl;
    }


    for (int i = 0; i < num_cons; i ++) {
        num_qvars = cqc_constrs[i]->vars.size();
        rhs_linear = 0;
        rhs_quad = 0;
        rhs_std.resize(num_qvars, 0);
        rhs_std.assign(num_qvars, 0);


        for (int j = 0; j <num_qvars; j++) {
            var_id = cqc_constrs[i]->vars[j];
            rhs_linear += x_vals[var_id] *  cqc_constrs[i]->linear[j];

            for(int k = 0; k <num_qvars; k++){
                var_id = cqc_constrs[i]->vars[k];
                rhs_std[j] += x_vals[var_id] *  cqc_constrs[i]->sigma[j][k];
            }
            rhs_quad += rhs_std[j]*rhs_std[j];
        }

        if (rhs_linear + cqc_constrs[i]->omega*sqrt(rhs_quad) <= cqc_constrs[i]->rhs) {
            cout << " - CQ Diagonal Constraint " << i << " is correct" << endl;
            cout << "\t "<< rhs_linear << " + " << cqc_constrs[i]->omega << " * sqrt( " << rhs_quad << " ) = " << rhs_linear + sqrt(rhs_quad);
            cout << " <= " << cqc_constrs[i]->rhs << endl;
        }
        else {
            cout << " - Error CQ Diagonal Constraint " << i << endl;
            cout << "\t "<< rhs_linear << " + " << cqc_constrs[i]->omega << " * sqrt( " << rhs_quad << " ) = " << rhs_linear + sqrt(rhs_quad);
            cout << " > " << cqc_constrs[i]->rhs << endl;
            assert(false);
        }


    }
}

void check_objective_function(IloCplex& cplex, IloIntVarArray& x,vector<int>& cost){

    int num_vars = cost.size();
    double obj_value = cplex.getObjValue();
    double obj_mine = 0;

    for (IloInt i = 0; i < num_vars; i++) {
        obj_mine += cplex.getValue(x[i])*cost[i];
    }

    if (abs(obj_value - obj_mine) < 0.001) {
        cout << " - Objective Value is consistent!" << endl;
    }
    else {
        cout << " - Error in objective value " << endl;
        cout << "\t solver = " << obj_value << endl;
        cout << "\t mine = " << obj_mine << endl;
        assert(false);
    }

}


