//
// Created by Margarita Castro on 2019-08-16.
//

#include "bdd_cuts.h"
#include <cmath>

using namespace std;


// Create a specific BDD
void bdd_cuts::create_bdd_diag(int id) {

    if (bdds_diag[id] != nullptr) return;

    if (debug) cout << "\n=== Creating BDD (Diagonal)" << id << " =====" << endl;

    //Create BDD
    bdds_diag[id] = new bdd_diag(cq_diag[id], id, bdd_width);

    bdds_diag[id]->set_debug(false);
    bdds_diag[id]->construct_bdd();
    //exit(1);

    //Initialize slack
    int num_vars = cq_diag[id]->vars.size();
    slack[id] = new double[num_vars];
    pi[id] = new double[num_vars];
    x_value[id] = new double[num_vars];

    // Add BDD construction time
    bdd_time += bdds_diag[id]->get_constr_time();

    if(debug) cout << "Time BDD " << id <<"  = " << bdd_time << endl;

}

void bdd_cuts::create_bdd_chance(int id) {

    if (bdds_chance[id] != nullptr) return;

    if (debug) cout << "\n=== Creating BDD (Chance Constraint) " << id << " =====" << endl;

    //Create BDD
    bdds_chance[id] = new bdd_chance(cq_chance[id], id, bdd_width);
    bdds_chance[id]->set_debug(false);
    bdds_chance[id]->construct_bdd();
    //exit(1);

    //Initialize pointer
    int num_vars = cq_chance[id]->vars.size();
    slack[id] = new double[num_vars];
    pi[id] = new double[num_vars];
    x_value[id] = new double[num_vars];

    // Add BDD construction time
    bdd_time += bdds_chance[id]->get_constr_time();

    if(debug) cout << "Time BDD " << id <<"  = " << bdd_time << endl;

}

int bdd_cuts::get_bdd_cut(IloNumArray& x_vals, double& pi0, bool multi_cuts, bool& finish, bool& equality) {

    bool bdd_separate = false;
    int bdd_id = -1;
    int num_cons = 0;

    // Construct all BDDs
    if (type_diag)  for (int i = 0; i < num_diag; i++) create_bdd_diag(i);
    else            for (int i = 0; i < num_chance; i++) create_bdd_chance(i);

    // Fill BDD list if already used
    if (bdd_constructed.empty()) {
        if (type_diag)  num_cons = num_diag;
        else            num_cons = num_chance;

        for (int i = 0; i < num_cons; i++) bdd_constructed.push_back(i);
    }

    // Sort BDD list
    if (!multi_cuts) get_sorted_bdd_list(bdd_constructed, x_vals);


    // Weak BDD cuts
    if (bdd_cut_type == 1 || bdd_cut_type == 2 || bdd_cut_type == 5) {
        for (int id : bdd_constructed) {
            bdd_separate = compute_weak_bdd_cut(id, x_vals, pi0);

            if (bdd_separate) {
                bdd_id = id;
                break;
            }
        }
    }

    // Projected sub-gradient BDD cuts
    if (bdd_cut_type == 3){
        for (int id : bdd_constructed) {
            bdd_separate = compute_projected_bdd_cut(id, x_vals, pi0);

            if (bdd_separate) {
                bdd_id = id;
                break;
            }
        }
    }

    // Target BDD cuts
    if (bdd_id == -1 && (bdd_cut_type == 4 || bdd_cut_type == 5) ){
        for (int id : bdd_constructed) {
            bdd_separate = compute_target_bdd_cut(id, x_vals, pi0, equality);

            if (bdd_separate) {
                bdd_id = id;
                break;
            }
        }
    }

    // General BDD cuts (if all the weak bdd cuts fail)
    if (bdd_id == -1 && (bdd_cut_type == 2 || bdd_cut_type == 6) ) {
        stop_weak_cuts = true;
        // Iterate over constructed BDDs
        for (int id : bdd_constructed) {
            bdd_separate = compute_cglp_bdd_cut(id, x_vals, pi0);

            if (bdd_separate) {
                bdd_id = id;
                break;
            }
        }
    }

    // Remove BDD id from list to avoid choosing them again - only for multi_cut setting
    if (multi_cuts) {
        if(bdd_id == -1){
            bdd_constructed.clear();
        }
        else {
            list<int>::iterator it = find(bdd_constructed.begin(), bdd_constructed.end(), bdd_id);
            it++;
            bdd_constructed.erase(bdd_constructed.begin(), it);
        }

        if (bdd_constructed.empty()) finish = true;
    }

    return bdd_id;
}

bool bdd_cuts::compute_weak_bdd_cut(int id, IloNumArray& x_vals, double& pi0) {

    start_time = clock();

    int num_vars = 0;
    bool bdd_separate;

    if (type_diag)  num_vars = cq_diag[id]->vars.size();
    else            num_vars = cq_chance[id]->vars.size();

    for (int i = 0; i < num_vars; i++) pi[id][i] = 0;

    if (type_diag) {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_diag[id]->vars[i]];
    }
    else {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_chance[id]->vars[i]];
    }

    if (debug) {
        if (type_diag) bdds_diag[id]->set_debug(true);
        else           bdds_chance[id]->set_debug(true);
    }

    if (type_diag) bdd_separate = bdds_diag[id]->separate_min_cut(x_value[id], pi[id], pi0);
    else           bdd_separate = bdds_chance[id]->separate_min_cut(x_value[id], pi[id], pi0);

    if(bdd_cut_type == 1) {
        // Only count time when we are using poorly weak cuts
        end_time = clock();

        sep_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
        sep_it++;
        sep_time_avg = sep_time/sep_it;

//        cout << sep_time << " " << sep_it << " " << sep_time_avg << endl;
    }

    return bdd_separate;
}

bool bdd_cuts::compute_cglp_bdd_cut(int id, IloNumArray& x_vals, double& pi0){

    start_time = clock();

    int num_vars = 0;
    bool bdd_separate = false;

    if (type_diag)  num_vars = cq_diag[id]->vars.size();
    else            num_vars = cq_chance[id]->vars.size();

    for (int i = 0; i < num_vars; i++) pi[id][i] = 0;

    if (type_diag) {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_diag[id]->vars[i]];
    }
    else {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_chance[id]->vars[i]];
    }

    if (bdds_lp[id] == nullptr) {
        bdds_lp[id] = new bdd_lp(id);
        if (type_diag)  bdds_lp[id]->create_lp_model(bdds_diag[id]);
        else            bdds_lp[id]->create_lp_model(bdds_chance[id]);
    }

    if (debug) bdds_lp[id]->set_debug(true);

    bdd_separate = bdds_lp[id]->solve_lp(x_value[id], pi[id], pi0);

    end_time = clock();

    sep_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
    sep_it++;
    sep_time_avg = sep_time/sep_it;

    //cout << sep_time << " " << sep_it << " " << sep_time_avg << endl;

    return bdd_separate;
}

bool bdd_cuts::compute_target_bdd_cut(int id, IloNumArray& x_vals, double& pi0, bool& equality) {

    start_time = clock();

    int num_vars = 0;
    bool bdd_separate = false;

    if (type_diag) num_vars = cq_diag[id]->vars.size();
    else num_vars = cq_chance[id]->vars.size();

    for (int i = 0; i < num_vars; i++) pi[id][i] = 0;

    if (type_diag) {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_diag[id]->vars[i]];
    } else {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_chance[id]->vars[i]];
    }


    if (debug) {
        if (type_diag)  bdds_diag[id]->set_debug(true);
        else            bdds_chance[id]->set_debug(true);
    }

    if (bdds_lp[id] == nullptr) {
        bdds_lp[id] = new bdd_lp(id);

        vector<double> center;
        bool overflow = false;
        if (type_diag)  {
            center = bdds_diag[id]->get_center_point(overflow);
            if (!overflow)  bdds_lp[id]->create_lp_model_target_cuts(bdds_diag[id], center);
            else            bdds_lp[id]->set_overflow();
        }
        else {
            center = bdds_chance[id]->get_center_point(overflow); //TODO - implement overflow here too
            if (!overflow)  bdds_lp[id]->create_lp_model_target_cuts(bdds_chance[id], center);
            else            bdds_lp[id]->set_overflow();
        }
    }

    if (debug) bdds_lp[id]->set_debug(true);

    //Skip if we didn't create the BDD due to overflow of an int
    if ( !bdds_lp[id]->get_overflow()) {
        bdd_separate = bdds_lp[id]->solve_lp_target_cuts(x_value[id], pi[id], pi0, equality);
    }

    end_time = clock();

    sep_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
    sep_it++;
    sep_time_avg = sep_time/sep_it;

    //cout << sep_time << " " << sep_it << " " << sep_time_avg << endl;

    return bdd_separate;
}

bool bdd_cuts::compute_projected_bdd_cut(int id, IloNumArray& x_vals, double& pi0){

    start_time = clock();

    int num_vars = 0;
    bool bdd_separate = false;

    if (type_diag)  num_vars = cq_diag[id]->vars.size();
    else            num_vars = cq_chance[id]->vars.size();

    for (int i = 0; i < num_vars; i++) pi[id][i] = 0;

    if (type_diag) {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_diag[id]->vars[i]];
    }
    else {
        for (int i = 0; i < num_vars; i++) x_value[id][i] = x_vals[cq_chance[id]->vars[i]];
    }

    if (debug) {
        if (type_diag) bdds_diag[id]->set_debug(true);
        else           bdds_chance[id]->set_debug(true);
    }

    bdd_separate = projected_gradient_routine(id, num_vars, pi0);

    end_time = clock();

    sep_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
    sep_it++;
    sep_time_avg = sep_time/sep_it;

//    cout << sep_time << " " << sep_it << " " << sep_time_avg << endl;

    return bdd_separate;
}



// Projected sub-gradient cuts code
bool bdd_cuts::projected_gradient_routine(int id, int num_vars, double& pi0){

    double delta = 0;
    double rho = 1;
    int max_it = 50;  // Change to something else afterwards
    double aux = 0;

    //Initialize pointers
    if (pi_t == nullptr) {
        int max_vars = 0;
        if (type_diag)  {
            for (int i = 0; i < num_diag; i++) {
                if (cq_diag[i]->vars.size() > max_vars) max_vars = cq_diag[i]->vars.size();
            }
        }
        else {
            for (int i = 0; i < num_chance; i++) {
                if (cq_chance[i]->vars.size() > max_vars) max_vars = cq_chance[i]->vars.size();
            }
        }

        pi_t = new double[max_vars];
        x_t = new int[max_vars];
        x_best = new int[max_vars];
        x_diff = new double[max_vars];
    }

    for (int i=0; i < num_vars; i++) pi_t[i] = 0;

    for (int t = 0; t < max_it; t++) {
        //Step 1: Get longest path from BDD
        if (type_diag) bdds_diag[id]->compute_longest_path(pi_t, x_t);
        else           bdds_chance[id]->compute_longest_path(pi_t, x_t);

        //Step 2: Save best inequality (pi) and longest path for now
        for (int i = 0; i < num_vars; i++) x_diff[i] = x_value[id][i] -  x_t[i];

        aux = 0;
        for (int i = 0; i < num_vars; i++) aux +=  pi_t[i]* x_diff[i];

        if (aux > delta) {
            delta = aux;

            for (int i = 0; i < num_vars; i++) {
                x_best[i] = x_t[i];
                pi[id][i] = pi_t[i];
            }
        }

        //Step 3: Update the gradient
        for (int i = 0; i < num_vars; i++) pi_t[i] = pi_t[i] + rho*x_diff[i];

        //Step 4: normalize gradient
        aux = 0;
        for (int i = 0; i < num_vars; i++) aux += pi_t[i]*pi_t[i];
        aux = 1/sqrt(aux);

        for (int i = 0; i < num_vars; i++) pi_t[i] = aux*pi_t[i];

    }

    if (delta >  0.00001) { // To avoid precision errors
        //Compute pi0
        pi0 = 0;
        for (int i = 0; i < num_vars; i++) pi0 += pi[id][i]*x_best[i];

        return true;
    }

    return false;
}


bool  bdd_cuts::lift_bdd_constraint(int id, double& pi0) {

    if (pi[id] == nullptr) {
        cout << "Error - You need to construct a BDD before lifting" << endl;
    }

    bool can_lift = true;
    int lift_id = 0;
    bool lift = false;

    if (debug) {
        cout << "\nStarting constraint:" << endl;
        print_inequality(id, pi0);

        if (type_diag) bdds_diag[id]->set_debug(true);
        else            bdds_chance[id]->set_debug(true);
    }

    while (true) {

        start_time = clock();

        if (type_diag)  can_lift = bdds_diag[id]->get_bdd_slacks(pi[id], slack[id], pi0);
        else            can_lift = bdds_chance[id]->get_bdd_slacks(pi[id], slack[id], pi0);

        end_time = clock();

        lift_time += double(end_time - start_time)/ CLOCKS_PER_SEC;
        lift_it++;
        lift_time_avg = lift_time/lift_it;

//        cout << lift_time << " " << lift_it << " " << lift_time_avg << endl;

        if (!can_lift) {
            if (debug) cout << "Can't lift constraint" << endl;
            break;
        }

        lift_id = choose_slack(id);

        if (lift_id == -1) {
            if (debug) cout << "All BDD slacks equal to 0!" << endl;
            break;
        }

        //Update coefficients and rhs
        pi[id][lift_id] += slack[id][lift_id];

        if (slack[id][lift_id] < 0) pi0 += slack[id][lift_id];

        lift = true;

        if (debug) {
            cout << "Lifted constraint:" << endl;
            print_inequality(id, pi0);
        }
    }

    return lift;
}


// Lifting procedure for cover inequalities
bool bdd_cuts::lift_cover_qcd(list<int>& cover, double& pi0, int id) {

    // Create BDD
    create_bdd_diag(id);

    int num_vars = cq_diag[id]->vars.size();

    if (pi[id] == nullptr) {
        pi[id] = new double[cq_diag[id]->vars.size()];
    }

    // Reset pi values
    for (int i = 0; i < num_vars; i++) pi[id][i] = 0;

    // Set cover values
    for (int i: cover) pi[id][i] = 1;

    return lift_bdd_constraint(id, pi0);
}

int bdd_cuts::choose_slack(int id) {

    //Choose minimum absolute value slack

    int num_vars = 0;
    int candidate = -1;
    double value = INT_MAX;

    if (type_diag)  num_vars = cq_diag[id]->vars.size();
    else            num_vars = cq_chance[id]->vars.size();

    for (int i = 0; i < num_vars; i++) {
        if (abs(slack[id][i]) <  value && abs(slack[id][i]) > 0.001) {
            value = abs(slack[id][i]);
            candidate = i;
        }
    }

    return candidate;
}

void bdd_cuts::get_sorted_bdd_list(list<int>& bdd_list, IloNumArray& x_vals) {

    int count = 0;

    for (int i : bdd_list){
        if (type_diag)  constraint_margin[count].first = get_margin_diagonal(i, x_vals);
        else            constraint_margin[count].first = get_margin_chance(i, x_vals);

        constraint_margin[count].second = i;
        count++;
    }

    //Sort vector
    sort(constraint_margin.begin(), constraint_margin.begin() + count, less<>() );

    //Change the elements according to new order
    list<int>::iterator it = bdd_list.begin();
    count = 0;

    while (it!= bdd_list.end()) {
        *it = constraint_margin[count].second;
        count++;
        it++;
    }

}

double bdd_cuts::get_margin_diagonal(int id, IloNumArray& x_vals) {

    double lhs_linear = 0;
    double lhs_diag = 0;
    int num_var = cq_diag[id]->vars.size();

    for (int i =0; i < num_var; i++) {
        lhs_linear += x_vals[cq_diag[id]->vars[i]] * cq_diag[id]->linear[i];
        lhs_diag += x_vals[cq_diag[id]->vars[i]] * cq_diag[id]->diag2[i];
    }

    return cq_diag[id]->rhs - ( lhs_linear + cq_diag[id]->omega*sqrt(lhs_diag) );
}

double bdd_cuts::get_margin_chance(int id, IloNumArray& x_vals){

    double lhs_linear = 0;
    double lhs_sigma = 0;
    double aux = 0;
    int num_var = cq_chance[id]->vars.size();

    for (int i =0; i < num_var; i++) {
        lhs_linear += x_vals[cq_chance[id]->vars[i]] * cq_chance[id]->linear[i];
        aux = 0;
        //Compute quadratic term
        for (int j = 0; j <num_var; j++) {
            aux += x_vals[cq_chance[id]->vars[i]] * cq_chance[id]->sigma[i][j];
        }
        lhs_sigma += aux*aux;
    }

    return cq_chance[id]->rhs - ( lhs_linear + cq_chance[id]->omega*sqrt(lhs_sigma) );
}


// Print inequality
void bdd_cuts::print_inequality(int id, double&pi0) {

    int num_vars = 0;

    if (type_diag)   num_vars = cq_diag[id]->vars.size();
    else            num_vars = cq_chance[id]->vars.size();

    for (int i = 0; i < num_vars; i++) {
        if (type_diag)  cout << pi[id][i] << "*x[" << cq_diag[id]->vars[i] << "] ";
        else            cout << pi[id][i] << "*x[" << cq_chance[id]->vars[i] << "] ";

        if (i < num_vars - 1) cout << "+ ";
    }

    cout << "<= " << pi0 << endl;

}

// ======================================
// Callback functions bdd-cover
//=======================================

IloCplex::Callback
BDDCoverLift(IloEnv env, IloIntVarArray _x, bdd_cuts* _bdd_lift, cover_cuts* _cover,  IloNum& _lift_count ) {
    return (IloCplex::Callback(new(env) BDDCoverLiftI(env, _x, _bdd_lift, _cover, _lift_count)));
}


//Main loop for cut generation algorithm
void BDDCoverLiftI::main() {

    if (!isAfterCutLoop()) return;
    if (getNnodes() > 0) return; // Add cuts only in the root node

    //Get value of variables
    IloEnv env = getEnv();
    IloNumArray x_vals(env, num_vars);
    getValues(x_vals, x);

    list<int> cover_set;

    id_con = cover->get_cover_cut(x_vals, cover_set);
    rhs = cover_set.size() -1;

    if (id_con == -1) {
        x_vals.end();
        abortCutLoop();
        return;
    }

    // Create BDD and lift inequality
    bool lifted = bdd_lift->lift_cover_qcd(cover_set,rhs,id_con);

    lift_count +=lifted;

    // Add Cut
    IloExpr lhs(env);
    int n = bdd_lift->cq_diag[id_con]->vars.size();

    for (int i = 0; i < n; i++) {
        lhs += x[bdd_lift->cq_diag[id_con]->vars[i]]*bdd_lift->pi[id_con][i];
    }

    add(lhs <= rhs).end();

    //Free memory
    lhs.end();
    x_vals.end();

    return;
}

// ======================================
// Callback functions bdd-cut-and-lift
//=======================================

IloCplex::Callback
BDDCutLift(IloEnv env, IloIntVarArray _x, bdd_cuts* _bdd, IloNum& _lift_count,  bool _type_diag, bool _lift, cplex_param* param) {
    return (IloCplex::Callback(new(env) BDDCutLiftI(env, _x, _bdd, _lift_count,  _type_diag, _lift, param)));
}


//Main loop for cut generation algorithm
void BDDCutLiftI::main() {

    if (!isAfterCutLoop()) return;
    if (getNnodes() > 0){
        if (root_cuts) return; // Add cuts only in the root node

        // Change to weak BDD cuts during search
        if (weak_cut_tree) bdd->set_weak_bdd_cuts();
    }

    //Get value of variables
    IloEnv env = getEnv();
    IloNumArray x_vals(env, num_vars);
    getValues(x_vals, x);

    int n = 0;
    int num_cuts = 1;
    bool finish = false;
    bool equality = false;

    if (multi_cut) num_cuts = bdd->get_num_constraints();

    //Print X
//    cout << "x = ";
//    for (int i = 0; i < num_vars; i++) cout << x_vals[i] << " ";
//    cout << endl;

    // Find appropriate BDD to create
    for(int j =0; j < num_cuts; j++){
        equality = false;
        finish = false;
        id_con = bdd->get_bdd_cut(x_vals, rhs, multi_cut, finish, equality);

        if (id_con == -1) break;

        cuts_count++;

        // Lift Inequality
        if (lift && !equality) {
            bool lifted = bdd->lift_bdd_constraint(id_con, rhs);
            lift_count +=lifted;
        }

        // Add Cut
        IloExpr lhs(env);

        if (type_diag) {
            n = bdd->cq_diag[id_con]->vars.size();
            for (int i = 0; i < n; i++) {
                lhs += x[bdd->cq_diag[id_con]->vars[i]] * bdd->pi[id_con][i];
            }
            //Print Cut
//            cout << "== Cut from BDD " << id_con << endl;
//            for (int i = 0; i < n; i++) {
//                cout << bdd->pi[id_con][i] <<"*x["<< bdd->cq_diag[id_con]->vars[i] << "] + ";
//            }
//            if(equality)    cout << "== " << rhs << endl;
//            else            cout << "<= " << rhs << endl;

            //Check violated constraint
//            int myints[] = {1,2,3,5,6,9,11,12,13,16,18,19,20,21,22,24,25,26,27,
//            28,29,33,36,38,40,41,42,45,46,50,51,52,53,56,57,58,60,61,64,70,71,72,74,
//            76,77,78,80,82,83,85,86,88,89,91,92,93,94,95,97,101,102,103,104,105,106,
//            109,111,112,113,115,118,119,121,123,126,128,129,135,137,139,142,143,146,
//            147,149};
//
//            std::vector<int> index (myints, myints + sizeof(myints) / sizeof(int) );
//
//            double value = 0;
//            for (int i = 0; i < index.size(); i++) {
//                for (int j =0; j < n; j++) {
//                    if(index[i] == bdd->cq_diag[id_con]->vars[j]){
//                        value += bdd->pi[id_con][j];
//                        break;
//                    }
//                }
//            }
//            if (equality) {
//                if(value != rhs){
//                    cout << "Error here - BDD id = " << id_con << endl;
//                    cout << value << "  " << rhs << endl;
//                    exit(1);
//                }
//            }
//            else {
//                if(value > rhs){
//                    cout << "Error here - BDD id = " << id_con << endl;
//                    cout << value << "  " << rhs << endl;
//                    exit(1);
//                }
//            }

        }
        else {
            n = bdd->cq_chance[id_con]->vars.size();
            for (int i = 0; i < n; i++) {
                lhs += x[bdd->cq_chance[id_con]->vars[i]] * bdd->pi[id_con][i];
            }

            //Print Cut
//            cout << "== Cut from BDD " << id_con << endl;
//            for (int i = 0; i < n; i++) {
//                cout << bdd->pi[id_con][i] <<"*x["<< bdd->cq_chance[id_con]->vars[i] << "] + ";
//            }
//
//            if(equality)    cout << "== " << rhs << endl;
//            else            cout << "<= " << rhs << endl;
        }

        if(equality)    add(lhs == rhs).end();
        else            add(lhs <= rhs).end();

        //Free memory
        lhs.end();

        if(finish) break;
    }

    //Free memory
    x_vals.end();

    return;
}