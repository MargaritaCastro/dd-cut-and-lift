//
// Created by Margarita Castro on 2019-12-18.
//

#include "bdd_lp.h"


void bdd_lp::create_lp_model(bdd_chance* bdd){

    char varname[256];
    IloEnv env = model.getEnv();

    // Get BDD num_vars and nodes
    num_vars = bdd->get_num_vars();
    const vector< vector <bdd_chance::bdd_node*> >& nodes = bdd->get_nodes();
    const vector< vector <bdd_chance::bdd_arc*> >& arcs = bdd->get_arcs();

    // Construct variables
    mu      = IloNumVarArray(env, num_vars, 0, IloInfinity);
    lambda  = IloNumVarArray(env, num_vars, 0, IloInfinity);
    w = IloNumVarArray2(env);

    for (int i = 0; i < num_vars; i++) {
        sprintf(varname, "mu_%d", i);
        mu[i].setName(varname);

        sprintf(varname, "lambda_%d", i);
        lambda[i].setName(varname);

        model.add(lambda[i] >= 0);
        model.add(mu[i] >= 0);
    }

    //Map nodes to variables
    vector<vector<int> > nodes_to_vars(num_vars + 1, vector<int>() );
    int count = 0;
    int total = 0;

    for (int l = 0; l < num_vars + 1; l++) {
        total = nodes[l].size();
        nodes_to_vars[l].resize(total, -1);
        count = 0;

        //map between variable index and node index.
        for (int i =0; i < total; i++) {
            if(nodes[l][i]->is_free) continue;

            nodes_to_vars[l][i] = count;
            count++;
        }
        w.add(IloNumVarArray(env, count, -IloInfinity, IloInfinity));

        for (int i = 0; i < count; i++) {
            sprintf(varname, "w_%d_%d", l, i);
            w[l][i].setName(varname);
        }
    }

    //Constraints: first layer and last layer
    model.add( w[0][0] == 1);
    model.add( w[num_vars][0] == 0);

    //Constraints: all layers (in_arcs to terminal)
    bdd_chance::bdd_arc* arc = nullptr;
    for (int l = 0; l < num_vars; l++) {
        total = arcs[l].size();
        for (int i =0; i < total; i++) {
            arc = arcs[l][i];
            if(arc->target == -1) continue;

            if(arc->value == 0)
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]] + mu[l] >= 0);
            else
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]] + lambda[l] >= 0);
        }
    }

    //Empty Objective
    obj.setSense(IloObjective::Minimize);
    model.add(obj);

    constructred = true;
}

void bdd_lp::create_lp_model(bdd_diag* bdd){

    char varname[256];
    IloEnv env = model.getEnv();

    // Get BDD num_vars and nodes
    num_vars = bdd->get_num_vars();
    const vector< vector <bdd_diag::bdd_node*> >& nodes = bdd->get_nodes();
    const vector< vector <bdd_diag::bdd_arc*> >& arcs = bdd->get_arcs();

    // Construct variables
    mu      = IloNumVarArray(env, num_vars, 0, IloInfinity);
    lambda  = IloNumVarArray(env, num_vars, 0, IloInfinity);
    w = IloNumVarArray2(env);

    for (int i = 0; i < num_vars; i++) {
        sprintf(varname, "mu_%d", i);
        mu[i].setName(varname);

        sprintf(varname, "lambda_%d", i);
        lambda[i].setName(varname);

        model.add(lambda[i] >= 0);
        model.add(mu[i] >= 0);
    }

    //Map nodes to variables
    vector<vector<int> > nodes_to_vars(num_vars + 1, vector<int>() );
    int count = 0;
    int total = 0;

    for (int l = 0; l < num_vars + 1; l++) {
        total = nodes[l].size();
        nodes_to_vars[l].resize(total, -1);
        count = 0;

        //map between variable index and node index.
        for (int i =0; i < total; i++) {
            if(nodes[l][i]->is_free) continue;

            nodes_to_vars[l][i] = count;
            count++;
        }
        w.add(IloNumVarArray(env, count, -IloInfinity, IloInfinity));

        for (int i = 0; i < count; i++) {
            sprintf(varname, "w_%d_%d", l, i);
            w[l][i].setName(varname);
        }
    }

    //Constraints: first layer and last layer
    model.add( w[0][0] == 1);
    model.add( w[num_vars][0] == 0);

    //Constraints: all layers (in_arcs to terminal)
    bdd_diag::bdd_arc* arc = nullptr;
    for (int l = 0; l < num_vars; l++) {
        total = arcs[l].size();
        for (int i =0; i < total; i++) {
            arc = arcs[l][i];
            if(arc->target == -1) continue;

            if(arc->value == 0)
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]] + mu[l] >= 0);
            else
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]] + lambda[l] >= 0);
        }
    }

    //Empty Objective
    obj.setSense(IloObjective::Minimize);
    model.add(obj);

    constructred = true;
}


bool bdd_lp::solve_lp(double* xval, double* pi, double& pi0) {

    update_objective(xval);

    if(!debug) cplex.setOut(cplex.getEnv().getNullStream());

    cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.solve();

    //Get Pi values
    pi0 = -1;
    for(int i =0; i < num_vars; i++){
        pi[i] = cplex.getValue(mu[i]) - cplex.getValue(lambda[i]);
        pi0 += cplex.getValue(mu[i]);
    }

    if(debug){
        cout << "Status:\t" << cplex.getStatus() << endl;
        cout << "Objective:\t" << cplex.getObjValue() << endl;

        cout << "Lambda = ";
        for (int i =0; i < num_vars; i++) {
            cout << cplex.getValue(lambda[i]) << " ";
        }
        cout << endl;

        cout << "Mu = ";
        for (int i =0; i < num_vars; i++) {
            cout << cplex.getValue(mu[i]) << " ";
        }
        cout << endl;
    }

    return abs(cplex.getObjValue() - 1) >= 0.001;
}

void bdd_lp::update_objective(double* xval){

    model.remove(obj);
    IloExpr objExpr = obj.getExpr();
    objExpr.clear();

    for(int i =0; i < num_vars; i++) {
        objExpr += xval[i]*lambda[i] + (1 - xval[i])*mu[i];
    }

    obj.setExpr(objExpr);
    model.add(obj);
    objExpr.end();

    if(debug){
        cplex.exportModel("bdd.lp");
        cout << "Xval = ";
        for (int i =0; i < num_vars; i++) {
            cout << xval[i] << " ";
        }
        cout << endl;
    }
}

//------------- FUNCTIONS TARGET CUTS -------------------//

void bdd_lp::create_lp_model_target_cuts(bdd_chance* bdd, vector<double>& _center) {

    char varname[256];
    IloEnv env = model.getEnv();

    //Setting center point
    center = _center;

    // Get BDD num_vars and nodes
    num_vars = bdd->get_num_vars();
    const vector< vector <bdd_chance::bdd_node*> >& nodes = bdd->get_nodes();
    const vector< vector <bdd_chance::bdd_arc*> >& arcs   = bdd->get_arcs();

    // Construct variables
    mu = IloNumVarArray(env, num_vars, -IloInfinity, IloInfinity);
    w  = IloNumVarArray2(env);

    mu_id.resize(num_vars, 0);

    for (int i = 0; i < num_vars; i++) {
        sprintf(varname, "mu_%d", i);
        mu[i].setName(varname);
        mu_id[i] = mu[i].getId();

        model.add( mu[i] >= -IloInfinity);
    }

    //Map nodes to variables
    vector<vector<int> > nodes_to_vars(num_vars + 1, vector<int>() );
    int count = 0;
    int total = 0;

    for (int l = 0; l < num_vars + 1; l++) {
        total = nodes[l].size();
        nodes_to_vars[l].resize(total, -1);
        count = 0;

        //map between variable index and node index.
        for (int i =0; i < total; i++) {
            if(nodes[l][i]->is_free) continue;

            nodes_to_vars[l][i] = count;
            count++;
        }
        w.add(IloNumVarArray(env, count, -IloInfinity, IloInfinity));

        for (int i = 0; i < count; i++) {
            sprintf(varname, "w_%d_%d", l, i);
            w[l][i].setName(varname);
        }
    }

    //Constraints: first layer
    IloExpr aux(env);
    aux += 1 - w[0][0];
    for (int i = 0; i < num_vars; i++) aux += mu[i]*center[i];

    model.add( aux <= 0);
    model.add( aux >= 0);

    //Constraints: last layer
    model.add( w[num_vars][0] == 0);

    //Constraints: all layers (in_arcs to terminal)
    bdd_chance::bdd_arc* arc = nullptr;
    for (int l = 0; l < num_vars; l++) {
        total = arcs[l].size();
        for (int i =0; i < total; i++) {
            arc = arcs[l][i];
            if(arc->target == -1) continue;

            if(arc->value == 0)
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]]  <= 0);
            else
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]] + mu[l] <= 0);
        }
    }

    //Empty Objective
    obj.setSense(IloObjective::Maximize);
    model.add(obj);

    constructred = true;
}

void bdd_lp::create_lp_model_target_cuts(bdd_diag* bdd, vector<double>& _center) {

    char varname[256];
    IloEnv env = model.getEnv();

    //Setting center point
    center = _center;

    // Get BDD num_vars and nodes
    num_vars = bdd->get_num_vars();
    const vector< vector <bdd_diag::bdd_node*> >& nodes = bdd->get_nodes();
    const vector< vector <bdd_diag::bdd_arc*> >& arcs   = bdd->get_arcs();

    // Construct variables
    mu = IloNumVarArray(env, num_vars, -IloInfinity, IloInfinity);
    w  = IloNumVarArray2(env);

    mu_id.resize(num_vars, 0);

    for (int i = 0; i < num_vars; i++) {
        sprintf(varname, "mu_%d", i);
        mu[i].setName(varname);
        mu_id[i] = mu[i].getId();

        model.add( mu[i] >= -IloInfinity);
    }

    //Map nodes to variables
    vector<vector<int> > nodes_to_vars(num_vars + 1, vector<int>() );
    int count = 0;
    int total = 0;

    for (int l = 0; l < num_vars + 1; l++) {
        total = nodes[l].size();
        nodes_to_vars[l].resize(total, -1);
        count = 0;

        //map between variable index and node index.
        for (int i =0; i < total; i++) {
            if(nodes[l][i]->is_free) continue;

            nodes_to_vars[l][i] = count;
            count++;
        }
        w.add(IloNumVarArray(env, count, -IloInfinity, IloInfinity));

        for (int i = 0; i < count; i++) {
            sprintf(varname, "w_%d_%d", l, i);
            w[l][i].setName(varname);
        }
    }

    //Constraints: first layer
    IloExpr aux(env);
    aux += 1 - w[0][0];
    for (int i = 0; i < num_vars; i++) aux += mu[i]*center[i];

    model.add( aux <= 0);
    model.add( aux >= 0);

    //Constraints: last layer
    model.add( w[num_vars][0] == 0);

    //Constraints: all layers (in_arcs to terminal)
    bdd_diag::bdd_arc* arc = nullptr;
    for (int l = 0; l < num_vars; l++) {
        total = arcs[l].size();
        for (int i =0; i < total; i++) {
            arc = arcs[l][i];
            if(arc->target == -1) continue;

            if(arc->value == 0)
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]]  <= 0);
            else
                model.add( w[l+1][nodes_to_vars[l+1][arc->target] ] - w[l][nodes_to_vars[l][arc->source]] + mu[l] <= 0);
        }
    }

    //Empty Objective
    obj.setSense(IloObjective::Maximize);
    model.add(obj);

    constructred = true;

}

bool bdd_lp::solve_lp_target_cuts(double* xval, double* pi, double& pi0, bool& equality) {

    update_objective_target_cuts(xval);

    //if(!debug)
        cplex.setOut(cplex.getEnv().getNullStream());

    cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.solve();

    if (cplex.getStatus() ==  IloAlgorithm::Optimal) {
        //Get Pi values of extreme point constraint
        pi0 = 1.00001;
        for(int i =0; i < num_vars; i++){
            pi[i] = cplex.getValue(mu[i]);
            pi0 += cplex.getValue(mu[i])*center[i];
        }
    }
    else if(cplex.getStatus() ==  IloAlgorithm::Unbounded) {
        equality = true;

        IloEnv env = model.getEnv();
        IloNumArray vals(env,num_vars);
        IloNumVarArray vars(env,num_vars, 0, IloInfinity);
        cplex.getRay(vals, vars);

        //Update pi values
        pi0 = 0;
        int var_id = 0;
        for (int i =0; i < num_vars; i++) pi[i] = 0;

        for (int i =0; i < vars.getSize(); i++) {
            if (vars[i].getId() >= mu_id[0] && vars[i].getId() <= mu_id[num_vars - 1]) {
                var_id = vars[i].getId() - mu_id[0];
                pi[var_id] = vals[i];
                pi0 += vals[i] * center[var_id];
            }
        }

        if(debug) cout << "Creating Unbounded cut" << endl;

        if(debug) {
            cout << "Ray Vals= ";
            for (int i = 0; i < vals.getSize(); i++) {
                cout << vals[i] << " ";
            }
            cout << "\nRay Vars = ";
            for (int i = 0; i < vars.getSize(); i++) {
                cout << vars[i].getName() << "," << vars[i].getId() << " ";
            }
            cout << "\nPi = ";
            for (int i =0; i < num_vars; i++) {
                cout << pi[i] << " ";
            }
            cout << "\nPi_0 = " << pi0 << endl;
        }

        //Free memory
        vals.end();
        vars.end();
    }

    if(debug){
        cout << "Status:\t" << cplex.getStatus() << endl;
        cout << "Objective:\t" << cplex.getObjValue() << endl;

        cout << "Mu = ";
        for (int i =0; i < num_vars; i++) {
            cout << cplex.getValue(mu[i]) << " ";
        }
        cout << endl;
    }

    if (cplex.getStatus() ==  IloAlgorithm::Unbounded) return true;

    return cplex.getObjValue() >= 1.001;
}

void bdd_lp::update_objective_target_cuts(double* xval) {

    model.remove(obj);
    IloExpr objExpr = obj.getExpr();
    objExpr.clear();

    for (int i =0; i < num_vars; i++) {
        objExpr += mu[i]*(xval[i] - center[i]);
    }

    obj.setExpr(objExpr);
    model.add(obj);
    objExpr.end();

    if(debug){
        cplex.exportModel("bdd.lp");
        cout << "Xval = ";
        for (int i =0; i < num_vars; i++) {
            cout << xval[i] << " ";
        }
        cout << endl;
    }


}