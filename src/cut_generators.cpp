//
// Created by Margarita Castro on 2019-08-05.
//

#include "cut_generators.h"

using namespace std;

// ======================================
// Callback functions - cover cuts
//=======================================

IloCplex::Callback
CoverCuts(IloEnv env, IloIntVarArray _x, cover_cuts* _cover) {
    return (IloCplex::Callback(new(env) CoverCutsI(env, _x, _cover)));
}

//Main loop for cover cut generation
void CoverCutsI::main() {

    if (!isAfterCutLoop()) return;
    if (getNnodes() > 0) return; // Add cuts only in the root node

    IloEnv env = getEnv();
    IloNumArray x_vals(env, x.getSize());
    getValues(x_vals, x);

    list<int> cover_set;

    int id_con = cover->get_cover_cut(x_vals, cover_set);

    if (id_con == -1) {
        abortCutLoop();
        return;
    }

    //Create cover inequality
    IloExpr lhs(env);
    IloNum rhs = cover_set.size() -1;

    for (int i: cover_set) {
        lhs += x[cover->cqd[id_con]->vars[i] ];
    }

    add(lhs <= rhs).end();

    //Free memory
    lhs.end();
    x_vals.end();

//    cout << "Cover set: ";
//    for (int i: cover_set) {
//        cout << i << " ";
//    }
//    cout << endl;
//
//    cout << "x: ";
//    for (int i = 0; i < num_vars; i++){
//        cout << x_vals[i] << " ";
//    }
//    cout << endl;

    return;
}


// ======================================
// Callback functions - lift cover cuts
//=======================================

IloCplex::Callback
LiftCoverCuts(IloEnv env, IloIntVarArray _x, cover_cuts* _cover, IloNum&  _lift_count) {
    return (IloCplex::Callback(new(env) LiftCoverCutsI(env, _x, _cover,_lift_count)));
}

//Main loop for cover cut generation
void LiftCoverCutsI::main() {

    if (!isAfterCutLoop()) return;
    if (getNnodes() > 0) return; // Add cuts only in the root node

    IloEnv env = getEnv();
    IloNumArray x_vals(env, x.getSize());
    getValues(x_vals, x);

    // Get minimal cover
    list<int> cover_set;
    int id_con = cover->get_cover_cut(x_vals, cover_set);

    if (id_con == -1) {
        x_vals.end();
        abortCutLoop();
        return;
    }

    //Create lifted constraint
    IloExpr lhs(env);
    IloNum rhs = cover_set.size() -1;

    if (cover->debug) {
        // Print cover constraints
        int count = 0;
        cout << "Cover constraint:\n";
        for (int i: cover_set) {
            cout << "x[" << cover->cqd[id_con]->vars[i] << "]";

            if (count < rhs) cout << " + ";
            count++;
        }
        cout << " <= " << rhs << endl;
    }

    // Lift cover constraint
    int num_vars = cover->cqd[id_con]->vars.size();
    cover->lift_cover_inequality(id_con, cover_set);

    if (cover->is_lift) lift_count++;

    for (int i = 0; i <num_vars; i++) {
        lhs += x[ cover->cqd[id_con]->vars[i] ]*cover->coeff[id_con][i];
    }

    add(lhs <= rhs).end();

    //Free memory
    lhs.end();
    x_vals.end();


    if (cover->debug){
        // Print lifted constraints
        int count = 0;
        cout << "Lifted cover constraint:\n";
        count = 0;
        for (int i = 0; i < num_vars; i++) {
            cout << cover->coeff[id_con][i] << "*x[" << cover->cqd[id_con]->vars[i] <<"]";
            if (count < num_vars - 1) cout << " + ";

            count++;
        }
        cout << " <= " << rhs << endl;
    }

}

// ======================================
// Get cover cut functions
//=======================================

int cover_cuts::get_cover_cut(IloNumArray& x, list<int>& cover_set) {

    bool found= false;

    for (int i = 0; i < num_con; i++) {
        cover_set.clear();
        found = find_cover_set(x, i, cover_set);

        if (found) return i;
    }

    return -1;
}


bool cover_cuts::find_cover_set(IloNumArray& x, int id, list<int>& cover_set) {

    start_time = clock();

    int num_vars = cqd[id]->vars.size();

    weights[id].assign(num_vars, make_pair(-1,-1));

    bool is_cover_set = false;
    bool to_return = false;

    // Iterate over all the combinations
    for (int i = 0; i < num_vars -1; i ++) {
        for (int j = i+1; j < num_vars; j ++) {
            cover_set.clear();
            is_cover_set = false;

            compute_sorted_weights( x, id, i, j);

            //Add variables to the cover inequality
            for (int k =0; k < num_vars; k++) {
                cover_set.push_back(weights[id][k].second);

                if (!is_cover_set){
                    is_cover_set = is_cover(id, cover_set);
                }

                if (is_cover_set) {
                    if (cuts_solution(x, id, cover_set)) {
                        // We found a cover cut that cuts off the solution!
                        to_return = true;
                        break;
                    }
                }
            }
            if (to_return) break;
        }
        if (to_return) break;
    }

    end_time = clock();

    sep_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
    sep_it++;
    sep_time_avg = sep_time/sep_it;

    if(!to_return) cover_set.clear();

    return to_return;
}

void cover_cuts::compute_sorted_weights(IloNumArray& x, int id, int i, int j){

    int num_vars = cqd[id]->vars.size();

    double lambda = 0.0;
    double mu = 0.0;

    //Compute lambda and mu
    if (cqd[id]->diag[i] == 0 && cqd[id]->diag[j] == 0)
        return;

    double xbar_i = 1 - x[cqd[id]->vars[i]];
    double xbar_j = 1 - x[cqd[id]->vars[j]];

    mu = ( cqd[id]->linear[i] * xbar_j - cqd[id]->linear[j]* xbar_i ) /
            ( cqd[id]->linear[i]*cqd[id]->diag2[j] - cqd[id]->linear[j]*cqd[id]->diag2[i]  ) ;
    lambda = ( (1 - xbar_i) - cqd[id]->diag2[i]*mu ) / cqd[id]->linear[i];


    //Compute weights
    for (int k = 0; k < num_vars; k ++) {
        weights[id][k].first = (1 - x[cqd[id]->vars[k]])/( cqd[id]->linear[k]*lambda + cqd[id]->diag2[k]*mu);
        weights[id][k].second = k;

    }

    sort(weights[id].begin(), weights[id].end());
}

bool cover_cuts::is_cover(int id,  list<int>& cover_set) {

    double linear = 0;
    double diagonal = 0;

    for (auto i: cover_set){
        linear += cqd[id]->linear[i];
        diagonal += cqd[id]->diag2[i];
    }

    return (linear + cqd[id]->omega * sqrt(diagonal) > cqd[id]->rhs);
}

bool cover_cuts::cuts_solution(IloNumArray& x, int id,  list<int>& cover_set) {

    double rhs = 0;
    double lhs = cover_set.size() - 1;

    for (auto i: cover_set){
        rhs += x[cqd[id]->vars[i]];
    }

    return (rhs > lhs + 0.001);

}

// ======================================
// Get  lifted cover cut functions
//=======================================

void cover_cuts::lift_cover_inequality(int id, list<int>& cover_set) {

    start_time = clock();

    list<int> remaining_vars;

    int num_vars = cqd[id]->vars.size();
    int cover_size = cover_set.size();

    coeff[id].assign(num_vars, 0);
    is_lift = false;

    //Set initial coefficients
    for (int i : cover_set) coeff[id][i] = 1;

    //Get remaining coefficients
    for (int i =0; i <num_vars; i++) {
        if (coeff[id][i] == 0) remaining_vars.push_back(i);
    }

    // Pre-fix variables
    fix_coefficients(id, cover_set, remaining_vars);

    // Lift remaining variables
    for (int i : remaining_vars){
        if (!is_cqp_created[id]) create_cqp(id);

        update_run_cqp(id, i, cover_size);
    }

    end_time = clock();

    lift_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
    lift_it++;
    lift_time_avg = lift_time/lift_it;

}


//Fix some lifting coefficients as in Atamturk(2009)
void cover_cuts::fix_coefficients(int id, list<int>& cover_set, list<int>& remaining_vars) {

    int num_vars = cqd[id]->vars.size();

    double mu = 0;
    double nu = 0;
    double diff = 0;
    int item = 0;

    list<int> empty_list;
    list<int> full_list;

    for (int i =0; i < num_vars; i++) full_list.push_back(i);

    get_mu_nu(id, cover_set, mu, nu);

    list<int>::iterator it = remaining_vars.begin();

    for ( ; it != remaining_vars.end(); ) {
        item = *it;
        diff_constraint(id, empty_list, *it, diff);

        //Check if first condition applies
        if (diff <= cqd[id]->rhs - nu) {
            //eliminate i from the remaining values
            it = remaining_vars.erase(it);
            continue;
        }

        //Check if second condition applies
        full_list.remove(item);
        diff_constraint(id, full_list, item, diff);
        full_list.push_back(item);

        if (diff >= mu) {
            // eliminate i from the remaining values
            it = remaining_vars.erase(it);
            coeff[id][item] = cover_set.size() - 1;
            is_lift = true;
            continue;
        }

        it++;
    }
}

//Computes mu_{|S| - 1} and nu_{|S| - 1}
void cover_cuts::get_mu_nu(int id, list<int>& S, double& mu, double& nu) {

    int j = 0;
    list <int> aux = S;
    double value = 0;

    mu = INT_MIN;   // max_value
    nu = INT_MAX;   // min_value

    //Compute values for the function
    for (int i: S) {
        aux.remove(i);

        evaluate_constraint(id, aux, value);

        //add element back
        aux.push_back(i);
        j++;

        //update mu and nu
        if (value > mu)
            mu = value;

        if (value < nu)
            nu = value;
    }

    value = 0;
}


void cover_cuts::create_cqp(int id) {

    if (is_cqp_created[id]) return;

    if (debug) cout << "=== Creating CQP for constraint " << id << " ====" << endl;

    char varname[256];

    env[id] = new IloEnv();
    cplex[id] = new IloCplex(*env[id]);

    IloModel model(*env[id]);
    cplex[id]->extract(model);

    if(!debug) cplex[id]->setOut(env[id]->getNullStream());
    cplex[id]->setParam(IloCplex::Param::Preprocessing::Reduce, 0);
    cplex[id]->setParam(IloCplex::Param::Threads, 1);

    // Create x variables
    IloInt num_vars = cqd[id]->vars.size();

    x[id] = new IloNumVarArray(*env[id],num_vars,0,1);

    for (int i =0; i <num_vars; i++) {
        sprintf(varname, "x_%d", i);
        (*x[id])[i].setName(varname);
    }
    model.add( *x[id]);

    // Create auxiliary variables
    IloNumVarArray y( *env[id], num_vars, 0, IloInfinity);

    for (int j = 0; j <num_vars; j++) {
        sprintf(varname, "y_%d", j );
        y[j].setName(varname);
    }
    model.add(y);

    sprintf(varname, "z");
    IloNumVar z( *env[id], 0, IloInfinity, ILOFLOAT, varname);
    model.add(z);

    //Add quadratic term
    model.add(IloScalProd(y,y) <= z*z);

    //Relate auxiliary variables to normal variables
    for (int j =0; j < num_vars; j++) {
        model.add(y[j] == cqd[id]->diag[j] * (*x[id])[j]);
    }

    // Add linear term
    IloExpr qConstraint(*env[id]);
    for (int j =0; j < num_vars; j++){
        qConstraint += cqd[id]->linear[j] * (*x[id])[j];
    }

    qConstraint += cqd[id]->omega * z;
    model.add(qConstraint <= cqd[id]->rhs );

    //Add dummy objective function
    obj[id] = new IloObjective(*env[id]);

    obj[id]->setSense(IloObjective::Maximize);
    model.add(*obj[id]);

    is_cqp_created[id] = true;

    //Export model
    if (debug) {
        sprintf(varname, "qcp_%d.lp",id);
        cplex[id]->exportModel(varname);
    }

}

// Update CQP model and solves it
void cover_cuts::update_run_cqp(int id,  int var, int cover_size) {

    if (debug)
        cout << "=== Solving CQP for constraint " << id << " and variable " << cqd[id]->vars[var] << " ====" << endl;

    IloModel model = cplex[id]->getModel();
    model.remove(*obj[id]);
    IloExpr objExpr = obj[id]->getExpr();
    objExpr.clear();

    // Update objective value
    IloInt num_vars = cqd[id]->vars.size();

    for (int i =0; i < num_vars; i++) {
        objExpr += (*x[id])[i] * coeff[id][i];
    }

    obj[id]->setExpr(objExpr);
    model.add(*obj[id]);
    objExpr.end();

    // Update bound of variables
    for (int i =0; i < num_vars; i++) {
        if (coeff[id][i] == 0) {
            //Set variable to zero
            (*x[id])[i].setBounds(0, 0);
        }
        else {
            //Set variable to normal bounds
            (*x[id])[i].setBounds(0, 1);
        }
    }

    //Set bounds of target variable
    (*x[id])[var].setBounds(1, 1);

    // Solve the worker LP
    cplex[id]->solve();

    //Update coefficient
    if (cplex[id]->getStatus() == IloAlgorithm::Optimal ) {
        double opt = cplex[id]->getObjValue();
        coeff[id][var] = cover_size -1 - floor(opt);

        if (coeff[id][var] < 0) {
            coeff[id][var] = 0;
        }
        else{
            is_lift = true;
        }

        if (debug) cout << "optimal = " << opt << endl;
    }


}

void cover_cuts::diff_constraint(int id,  list<int>& S, int i, double& diff) {

    double value1 = 0;
    double value2 = 0;

    //Evaluate  S
    evaluate_constraint(id, S, value2);

    //Evaluate S \cup \{ i \}
    S.push_back(i);
    evaluate_constraint(id, S, value1);
    S.pop_back();

    diff = value1 - value2;
}

void cover_cuts::evaluate_constraint(int id,  list<int>& S, double& value) {

    double linear = 0;
    double diag = 0;

    for (int i: S) {
        linear += cqd[id]->linear[i];
        diag += cqd[id]->diag2[i];
    }

    value = linear + cqd[id]->omega*sqrt(diag);
}
