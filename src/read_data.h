//
// Created by Margarita Castro on 2019-08-02.
//

#ifndef DD_CONIC_LIFTING_READ_DATA_H
#define DD_CONIC_LIFTING_READ_DATA_H


#include <cassert>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>

using namespace std;

//Parameter structure
struct bdd_param {
    int cut_type;
    bool lift;
    int width;

    bdd_param() {
        cut_type = 0;
        lift = 0;
        width = 0;
    }
};

struct cplex_param {
    bool cplex_cuts;
    bool multi_cut;
    bool root_cuts;
    bool weak_tree_cuts;
    int strategy;

    cplex_param() {
        cplex_cuts = true;
        root_cuts = true;
        multi_cut = false;
        weak_tree_cuts = false;
        strategy = 0;
    }
};

/*
 * Structure for Conic Quadratic Diagonal Constraints
 */
struct CQDConstraint {

    //Index of variables involve in constraints
    vector<int> vars;
    //Linear coefficients
    vector<int> linear;
    //Diagonal coefficients
    vector<int> diag;
    //Diagonal coefficients square
    vector<int> diag2;
    //LHS constant
    int rhs;
    //Omega constant
    int omega;
    //linear dominates
    bool dom_linear;

    // GUB extra parameters
    // Number of variables in each GUB group (follows the order in vars)
    vector<int> gub_groups;
    // Using GUB or nor
    bool gub;

    // Constructor
    CQDConstraint(vector<int> _vars, vector<int> _linear, vector<int> _diag, int _rhs, int _omega)
            : vars(_vars), linear(_linear), diag(_diag), rhs(_rhs), omega(_omega)
    {
        int num_vars = diag.size();
        dom_linear = true;

        diag2.resize(num_vars);

        for (int i =0; i < num_vars; i++){
            diag2[i] = pow(diag[i], 2);
        }

        gub = false;
    }

    ~CQDConstraint()= default;


    void set_linear_dominant(){
        dom_linear = true;
    }
    void set_diag_dominant(){
        dom_linear = false;
    }

    void set_gub_groups(vector<int> _gub_groups) {
        gub_groups = _gub_groups;
        gub = true;
    }
};


/*
 * Structure for Conic Quadratic Chance Constraints
 */
struct CQChanceConstraint {

    //Index of variables involve in constraints
    vector<int> vars;
    //Linear coefficients
    vector<int> linear;
    //Diagonal coefficients
    vector< vector<int> > sigma;
    //LHS constant
    int rhs;
    //Omega constant
    int omega;
    //linear domainat
    bool dom_linear;

    // Constructor
    CQChanceConstraint(vector<int> _vars, vector<int> _linear, vector< vector<int> > _std, int _rhs, int _omega)
            : vars(_vars), linear(_linear), sigma(_std), rhs(_rhs), omega(_omega)
    {
        dom_linear = true;
    }

    ~CQChanceConstraint()= default;

    void set_linear_dominant(){
        dom_linear = true;
    }
    void set_std_dominant(){
        dom_linear = false;
    }
};

void read_conic_diagonal_input(char* filename, vector<CQDConstraint*>& CQDConstraints, int& num_vars, vector<int>& cost);
void read_conic_chance_input(char* filename, vector<CQChanceConstraint*>& cqc_constr, int& num_vars, vector<int>& cost);

//Sorting functions
void sort_vars_cqd(CQDConstraint* cqd);
void sort_vars_gub(vector< vector<int> >& gub, CQDConstraint* CQDConstraints);

#endif //DD_CONIC_LIFTING_READ_DATA_H
