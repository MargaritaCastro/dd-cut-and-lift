// ------------------------------------------------------------
// BDD lifting and cut generation for SOC
// ------------------------------------------------------------

#include <iostream>
#include "mip_models.h"

using namespace std;

//
// Main
//
int main(int argc, char* argv[]) {

    if (true) {
        cout << "\nUsage: bdd_conic [filename (soc)] [flags] \n";
        cout << "\tflags:\n";
        cout << "\t\t-d = SOC with diagonal matrix \n";
        cout << "\t\t-g = SOC with general matrix \n";
        cout << "\t\t-c = Turn off CPLEX cuts\n";
        cout << "\t\t-s# = CPLEX startegy; 0 (default), 1(QCP) or 2(LP) \n";
        cout << "\t\t-bw = Activate Weak BDD cuts \n";
        cout << "\t\t-bg = Activate General BDD cuts \n";
        cout << "\t\t-bgw = Activate Weak + General BDD cuts \n";
        cout << "\t\t-bp = Activate Projected Sub-gradient BDD cuts \n";
        cout << "\t\t-bt = Activate Target BDD cuts \n";
        cout << "\t\t-btw = Activate Weak + Target BDD cuts \n";
        cout << "\t\t-bl = Activate BDD lifting cuts (only works in conjuction with -u1, -bw, or -bg flags)\n";
        cout << "\t\t-u# = Cover cuts (only for SOC knapsack). 1: Cover cuts, 2: Cover cuts and continuous lifting \n";
        cout << "\t\t-w# = BDD max_width. 0: no BDD, -1: let the code decide \n";
        cout << "\t\t-m# = Activate adding multiple BDD cuts in each iteration (default = single cut) \n";
        cout << "\t\t-t  = Add user cuts at the tree search \n";
        cout << "\t\t-tw  = Add Weak BDD cuts at the tree search \n";
        cout << "\t\t-o[filename] = Name of output file for results";
        cout << "\n\n";
    }

    cplex_param* cplex_p = new cplex_param();
    bdd_param* bdd_p = new bdd_param();

    int cover_cuts = 0;
    char* filename = argv[1];
    string outputname = "";

    int argcount = 2;

    int type = 0;

    while(argcount < argc) {
        if(argv[argcount][0] == '-') {
            switch(argv[argcount][1]) {
                case 'b':
                {
                    if (argv[argcount][2] == 'l') {
                        // Activate bdd lifting procedure
                        bdd_p->lift = true;
                    }
                    else if (argv[argcount][2] == 'w') {
                        // Activate weak bdd cuts
                        bdd_p->cut_type = 1;
                    }
                    else if (argv[argcount][2] == 'g') {
                        // Activate general bdd cuts
                        bdd_p->cut_type = 6;

                        // Activate weak + target bdd cuts
                        if(argv[argcount][3] == 'w') {
                            bdd_p->cut_type = 2;
                        }
                    }
                    else if (argv[argcount][2] == 'p') {
                        // Activate projected bdd cuts
                        bdd_p->cut_type = 3;
                    }
                    else if (argv[argcount][2] == 't') {
                        // Activate target bdd cuts
                        bdd_p->cut_type = 4;

                        // Activate weak +  target bdd cuts
                        if(argv[argcount][3] == 'w') {
                            bdd_p->cut_type = 5;
                        }
                    }
                    else{
                        cout << " no flag named " << argv[argcount] << endl;
                        exit(1);
                    }

                    break;
                }
                case 'c':
                {
                    //Turn off CPLEX cuts
                    cplex_p->cplex_cuts = false;
                    break;
                }
                case 'd':
                {
                    //SOC with diagonal matrix
                    type = 0;
                    break;
                }
                case 'g':
                {
                    //SOC with general matrix
                    type = 1;
                    break;
                }
                case 'm':
                {
                    //activate multiple BDD cuts
                    cplex_p->multi_cut = true;
                    break;
                }
                case 't':
                {
                    if (argv[argcount][2] == 'w') {
                        // Only weak cuts in the tree search
                        cplex_p->weak_tree_cuts = true;
                    }
                    //add cuts also at the tree search
                    cplex_p->root_cuts = false;
                    break;
                }
                case 'o':
                {
                    // output file name
                    outputname = &(argv[argcount][2]);
                    break;
                }
                case 's':
                {
                    //CPLEX strategy
                    cplex_p->strategy = atoi(&(argv[argcount][2]));
                    break;
                }
                case 'u':
                {
                    // Activate user Cover cuts for knapsack SOC
                    cover_cuts = atoi(&(argv[argcount][2]));
                    break;
                }
                case 'w':
                {
                    // set BDD maximum Width
                    bdd_p->width = atoi(&(argv[argcount][2]));
                    break;
                }
                default:
                {
                    cout << "No new parameters" << endl;
                    break;
                }
            }
        }
        ++argcount;
    }


    // === Check correct use of flags
    //If using any BDD flag, you need to specify a width
    if ( (bdd_p->lift || bdd_p->cut_type > 0) && bdd_p->width == 0)  {
        cout << "Error: Need to specify a BDD width to use the BDD cutting and/or lifting procedures " << endl;
        exit(1);
    }
    // bdd_lifting restrictions
    if (bdd_p->lift && cover_cuts == 2)  {
        cout << "Error: Can't use BDD lifting on top of the continuous lifting procedure" << endl;
        exit(1);
    }
    // bdd_lifting can be used only with other type of  user cuts
    if (bdd_p->lift && cover_cuts != 1 && bdd_p->cut_type == 0)  {
        cout << "Error: You need to specify a user cut (bdd or cover) to use the bdd lifting on" << endl;
        exit(1);
    }
    // Check correct entry of strategy
    if (!(cplex_p->strategy == 0 || cplex_p->strategy == 1 || cplex_p->strategy == 2) ){
        cout << "Error: use either -s1 or -s2 to set the cplex solving strategy" << endl;
        exit(1);
    }
    // Check correct entry of cover cuts
    if (!(cover_cuts == 0 || cover_cuts == 1 || cover_cuts == 2) ){
        cout << "Error: use either -u1 or -u2 for cover cuts" << endl;
        exit(1);
    }
    // Cover cuts only for diagonal case
    if (cover_cuts > 0 && type == 1 ){
        cout << "Error: can't use cover cuts for general SOC constraints" << endl;
        exit(1);
    }


    // === Run MIP  models
    if (type == 0)  mip_qc_diagonal(filename, outputname, cover_cuts, bdd_p, cplex_p);
    else            mip_qc_chance(filename, outputname, bdd_p, cplex_p);

    return 0;
}


