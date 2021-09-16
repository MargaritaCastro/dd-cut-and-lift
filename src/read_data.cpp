//
// Created by Margarita Castro on 2019-08-02.
//

#include <iostream>
#include <fstream>
#include <tuple>
#include <algorithm>
#include "read_data.h"

// ========================================================
// READ - CONIC DIAGONAL (+ GUB)
// ========================================================

void read_conic_diagonal_input(char* filename, vector<CQDConstraint*>& CQDConstraints, int& num_vars, vector<int>& cost) {

    ifstream input(filename);
    if( !input.is_open() ) {
        cerr << "Error: could not open file " << filename << endl;
        exit(1);
    }

    int num_cons = 0;
    int omega = 0;

    input >> num_vars;
    input >> num_cons;
    input >> omega;

    vector<int> b (num_cons,0);
    vector< vector<int> > a (num_cons, vector<int>(num_vars, 0) );
    vector< vector<int> > d (num_cons, vector<int>(num_vars, 0) );

    cost.resize(num_vars);

    // Read cost vector
    for (int i = 0; i < num_vars; i++) input >> cost[i];

    // Read lhs vector
    for (int i = 0; i < num_cons; i++) input >> b[i];

    // Read linear matrix
    for (int i = 0; i < num_cons; i++) {
        for (int j = 0; j < num_vars; j++) input >> a[i][j];
    }

    // Read diagonal matrix
    for (int i = 0; i < num_cons; i++) {
        for (int j = 0; j < num_vars; j++) input >> d[i][j];
    }

    //Create Constraint structures
    CQDConstraints.resize(num_cons, nullptr);

    int nvars = 0;
    int index = 0;
    vector<int> vars;
    vector<int> linear;
    vector<int> diagonal;

    for (int i =0; i < num_cons; i ++){

        //Count variables
        nvars = 0;
        for (int j = 0; j < num_vars; j++) {
            if (a[i][j] > 0) nvars++;
        }

        vars.resize(nvars, 0);
        linear.resize(nvars, 0);
        diagonal.resize(nvars,0);

        //Fill values of vectors
        index = 0;
        for (int j = 0; j < num_vars; j++) {
            if (a[i][j] > 0){
                vars[index] = j;
                linear[index] = a[i][j];
                diagonal[index] = d[i][j];
                index++;
            }
        }

        CQDConstraints[i] = new CQDConstraint(vars, linear, diagonal, b[i], omega);
        sort_vars_cqd(CQDConstraints[i]);
    }

    input.close();
}

// ========================================================
// READ - CHANCE CONSTRAINT
// ========================================================


void read_conic_chance_input(char* filename, vector<CQChanceConstraint*>& cqc_constr, int& num_vars, vector<int>& cost){

    ifstream input(filename);
    if( !input.is_open() ) {
        cerr << "Error: could not open file " << filename << endl;
        exit(1);
    }

    int num_cons = 0;
    int omega = 0;

    input >> num_vars;
    input >> num_cons;
    input >> omega;

    vector<int> b (num_cons,0);
    vector< vector<int> > a (num_cons, vector<int>(num_vars, 0) );
    vector< vector< vector<int> > > s (num_cons, vector< vector<int> >(num_vars, vector<int>(num_vars, 0)) );

    cost.resize(num_vars);


    // Read cost vector
    for (int i = 0; i < num_vars; i++) {
        input >> cost[i];
    }

    // Read lhs vector
    for (int i = 0; i < num_cons; i++) {
        input >> b[i];
    }

    // Read linear matrix
    for (int i = 0; i < num_cons; i++) {
        for (int j = 0; j < num_vars; j++) {
            input >> a[i][j];
        }
    }

    // Read sigma matrix
    for (int i = 0; i < num_cons; i++) {
        for (int j = 0; j < num_vars; j++) {
            for (int k = 0; k < num_vars; k++) {
                input >> s[i][j][k];
            }
        }
    }

    //Create Constraint structures
    cqc_constr.resize(num_cons, nullptr);

    int nvars = 0;
    int index = 0;
    vector<int> vars;
    vector<int> linear;
    vector< vector<int> > std;

    for (int i =0; i < num_cons; i ++) {

        //Count variables
        nvars = 0;
        for (int j = 0; j < num_vars; j++) {
            if (a[i][j] != 0) nvars++;
            else {
                for (int k = 0; k < num_vars; k++) {
                    if (s[i][j][k] != 0) {
                        nvars++;
                        break;
                    }
                }
            }
        }

        // Resize vectors
        vars.resize(nvars, 0);
        linear.resize(nvars, 0);
        std.resize(nvars, vector<int>(nvars,0) );

        for (int j = 0; j < nvars; j++) std[j].resize(nvars, 0);

        //Fill values of vectors
        index = 0;
        for (int j = 0; j < num_vars; j++) {
            if (a[i][j] != 0) {
                vars[index] = j;
                linear[index] = a[i][j];
                index++;
            }
            else{
                for (int k = 0; k < num_vars; k++) {
                    if (s[i][j][k] != 0) {
                        vars[index] = j;
                        linear[index] = a[i][j];
                        index++;
                        break;
                    }
                }
            }
        }

        //Fill values of matrix
        for (int j = 0; j < nvars; j++) {
            for (int k = j; k < nvars; k++) {
                std[j][k] = s[i][ vars[j] ][ vars[k] ];
                std[k][j] = s[i][ vars[k] ][ vars[j] ];
            }
        }

        cqc_constr[i] = new CQChanceConstraint(vars, linear, std, b[i], omega);

    }

}

// ========================================================
// SORTING FUNCTIONS
// ========================================================

void sort_vars_cqd(CQDConstraint* cqd){

    //Check dominate part : linear or diagonal
    int nvars = cqd->vars.size();
    int max_linear = 0;
    double max_diagonal = 0;

    for(int i =0 ; i < nvars; i++){
        max_linear += cqd->linear[i];
        max_diagonal += cqd->diag2[i];
    }

    max_diagonal = cqd->omega*sqrt(max_diagonal);

    vector< tuple<int, int, int> > sorting(nvars, make_tuple(-1,-1,-1));

    if (max_linear >= max_diagonal) {
        cqd->set_linear_dominant();

        //Linear first
        for (int i =0 ; i < nvars; i++) {
            get<0>(sorting[i]) = cqd->linear[i];
            get<1>(sorting[i]) = cqd->diag[i];
            get<2>(sorting[i]) = cqd->vars[i];
        }

        sort(sorting.begin(), sorting.end(), greater<tuple<int,int,int> >() );

        //Set values back
        for (int i =0 ; i < nvars; i++) {
            cqd->linear[i] = get<0>(sorting[i]);
            cqd->diag[i] = get<1>(sorting[i]);
            cqd->vars[i] = get<2>(sorting[i]);
            cqd->diag2[i] = pow(cqd->diag[i], 2);
        }
    }
    else {
        cqd->set_diag_dominant();

        //Diagonal first
        for (int i =0 ; i < nvars; i++) {
            get<0>(sorting[i]) = cqd->diag[i];
            get<1>(sorting[i]) = cqd->linear[i];
            get<2>(sorting[i]) = cqd->vars[i];
        }

        sort(sorting.begin(), sorting.end(), greater<tuple<int,int,int> >() );

        //Set values back
        for (int i =0 ; i < nvars; i++) {
            cqd->diag[i] = get<0>(sorting[i]);
            cqd->linear[i] = get<1>(sorting[i]);
            cqd->vars[i] = get<2>(sorting[i]);
            cqd->diag2[i] = pow(cqd->diag[i], 2);
        }
    }

}

void sort_vars_gub(vector< vector<int> >& gub, CQDConstraint* cqd){

    int num_gub = gub.size();
    int nvars = 0;
    int index = -1;
    int count = 0;
    int count_group = 0;

    vector<int> group_sizes;

    std::vector<int>::iterator it;

    cqd->set_linear_dominant();

    // === Order variables for each group
    // 0 : group id, 1: linear, 2: diag, 3: index
    vector< tuple<int, int, int, int> > sorting(cqd->vars.size(), make_tuple(-1,-1,-1,-1));

    //Sort variables for each group in the DD
    for (int k = 0; k < num_gub; k++) {
        nvars = gub[k].size();
        count_group = 0;

        for (int i = 0; i < nvars; i++) {
            //Check is variable is part of the constraint
            it = find(cqd->vars.begin(), cqd->vars.end(), gub[k][i]);

            if (it == cqd->vars.end()) {
                continue;
            }

            index = std::distance(cqd->vars.begin(), it);

            get<0>(sorting[count]) = -k;
            get<1>(sorting[count]) = cqd->linear[index];
            get<2>(sorting[count]) = cqd->diag[index];
            get<3>(sorting[count]) = cqd->vars[index];

            count++;
            count_group++;
        }
        if (count_group > 0) group_sizes.push_back(count_group);
    }

    sort(sorting.begin(), sorting.end(), greater<tuple<int,int,int, int> >() );

    //Set values back
    for (int i =0 ; i < cqd->vars.size(); i++) {
        cqd->linear[i] = get<1>(sorting[i]);
        cqd->diag[i] = get<2>(sorting[i]);
        cqd->vars[i] = get<3>(sorting[i]);
        cqd->diag2[i] = pow(cqd->diag[i], 2);
    }

    //Set groups in constraint
    cqd->set_gub_groups(group_sizes);
}