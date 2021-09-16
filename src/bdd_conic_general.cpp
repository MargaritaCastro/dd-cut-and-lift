//
// Created by Margarita Castro on 2019-10-01.
//


#include "bdd_conic_general.h"

#include <fstream>
#include <algorithm>
#include <queue>


// Construct BDD from scratch
void bdd_chance::construct_bdd() {

    if (is_constructed){
        cout << "BDD is already constructed" << endl;
        return;
    }

    clock_t start_constr = clock();

    //Initialize space and all necessary variable
    initialize_bdd_space();

    // Create initial bdd
    width_one_bdd();

    procedure_bottom_up(0);

    int i = 1;
    bdd_has_change = true;

    // Iterative procedures
    while (bdd_has_change && i <= 10) {
        bdd_has_change = false;

        procedure_top_down(i);

        if (!bdd_has_change) break;

        bdd_has_change = false;

        procedure_bottom_up(i);

        i++;
    }

    //Reduce final BDD
    reduce_bdd();

    is_constructed = true;

    clock_t end_constr = clock();
    constr_time = double(end_constr - start_constr) / CLOCKS_PER_SEC;

}


void bdd_chance::initialize_bdd_space() {

    int size = 0;

    // Create and initialize nodes
    nodes.resize(num_vars + 1, vector<bdd_node*>());

    for (int l=0 ; l < num_vars + 1; l++) {
        if (l == 0 || l == num_vars) {
            nodes[l].resize(1, nullptr);
            nodes[l][0] = new bdd_node(num_vars);
        }
        else {
            size = 2*nodes[l-1].size();
            size = min(size, max_width);
            nodes[l].resize(size, nullptr);

            for (int i = 0; i < size ; i ++){
                nodes[l][i] = new bdd_node(num_vars);
            }
        }
    }

    // Create and initialize arcs
    arcs.resize(num_vars, vector<bdd_arc*>());

    for (int l=0 ; l < num_vars; l++) {
        if (l == 0){
            arcs[l].resize(2, nullptr);
            arcs[l][0] = new bdd_arc(l, 0, 0);
            arcs[l][1] = new bdd_arc(l, 0, 1);
        }
        else{
            size = 2*nodes[l].size();
            arcs[l].resize(size, nullptr);

            for (int i = 0; i < size ; i ++) {
                arcs[l][i] = new bdd_arc(l, (int)floor(i/2), i%2);
            }
        }
    }

    // Initialize nodes available
    nodes_available.resize(num_vars + 1,max_width);

    for (int l=0 ; l < num_vars + 1; l++) {
        nodes_available[l] = nodes[l].size();
    }

    first_node_available.resize(num_vars + 1, 0);

    // Initialize arc & node candidate
    candidate_arcs.resize(max_width*2, make_tuple(-1,-1,-1));
    candidate_nodes.resize(max_width, -1);

    // Initialize exact layer
    is_layer_exact.resize(num_vars + 1, false);

    // Initialize edge cost vector (for slack computation)
    cost_edge.resize(2);

}

void bdd_chance::width_one_bdd() {
    // Update root node
    nodes[0][0]->is_exact = true;
    nodes[0][0]->is_free = false;

    nodes_available[0] = 0;
    first_node_available[0] = -1;

    // Update every layer - one node and two arcs per layer
    for (int l = 0; l < num_vars ; l++) {

        //Direct arcs of first node
        arcs[l][0]->target = 0;
        arcs[l][1]->target = 0;

        // Activate first node in next layer
        nodes[l+1][0]->is_free = false;

        nodes[l+1][0]->in_arcs.push_back(0);
        nodes[l+1][0]->in_arcs.push_back(1);

        update_node_top_down(l+1, 0);

        //Update nodes available
        nodes_available[l+1]--;
        first_node_available[l+1] = 1;
    }

    //Update terminal node bottom-up information
    nodes[num_vars][0]->lin_min_bottom = 0;
    nodes[num_vars][0]->lin_max_bottom = 0;


    nodes[num_vars][0]->sigma_min_bottom.assign(num_vars, 0);
    nodes[num_vars][0]->sigma_max_bottom.assign(num_vars, 0);

    first_node_available[num_vars] = -1;
    is_layer_exact[0] = true;

    if (debug) {
        char name[256];
        sprintf(name, "graphs/bddc%d_width0.gml", bdd_id);
        export_graph(name);
    }

}

void bdd_chance::procedure_top_down(int round) {

    int size = 0;

    if (debug) {
        cout << "\n== Starting top-down procedure "<< round <<" ==\n";
    }

    for (int l = 0; l < num_vars; l++) {
        size = arcs[l].size();

        //Filter arcs in current layer
        for (int i = 0; i < size; i++) {
            if (arcs[l][i]->target == -1) continue;

            if (is_arc_infeasible(arcs[l][i])) {
                remove_arc(arcs[l][i]);
                if (debug) cout << "Filter arc " << i << " in layer " << l << endl;
            }
        }

        // Check nodes and see if there is any available
        size = nodes[l + 1].size();
        for (int i = 0; i < size; i++) {
            if (nodes[l + 1][i]->is_free) continue;
            if (nodes[l + 1][i]->in_arcs.empty() ) {
                remove_node(l + 1, i);
                continue;
            }
            if (l + 1 < num_vars && arcs[l + 1][2*i]->target == -1 && arcs[l + 1][2*i+1]->target == -1){
                remove_node(l + 1, i);
                continue;
            }
        }

        // Refine nodes in next layer
        if (nodes_available[l + 1] >= 1 && !is_layer_exact[l + 1]) refine(l + 1);

        // Update nodes in next layer
        for (int i = 0; i < size; i++) {
            if (nodes[l+1][i]->is_free || nodes[l+1][i]->is_exact) continue;

            if (nodes[l+1][i]->in_arcs.empty()) {
                remove_node(l + 1, i);
                continue;
            }

            update_node_top_down(l+1, i);
        }

        check_exact_layer(l+1);
    }

    if (debug){
        cout << "\n== End top-down procedure "<< round <<" ==\n";

        char name[256];
        sprintf(name, "graphs/bdd%d_top%d.gml", bdd_id,round);
        export_graph(name);
    }

}


void bdd_chance::procedure_bottom_up(int round) {

    if (debug){
        cout << "\n== Starting bottom-up procedure "<< round <<" ==\n";
    }

    int size = 0;

    for (int l = num_vars-1; l >= 0; l--) {
        size = arcs[l].size();

        //Filter arcs is the current layer
        for (int i = 0; i < size; i++) {
            if (arcs[l][i]->target == -1) continue;

            if (is_arc_infeasible(arcs[l][i])) {
                remove_arc(arcs[l][i]);
                if(debug) cout << "Filter arc " << i << " in layer " << l << endl;
            }
        }

        //Reduce nodes in the layer
        reduce_layer(l, false);

        // Update nodes in the current layer
        size = nodes[l].size();
        for (int i = 0; i < size; i++) {
            if (nodes[l][i]->is_free) continue;
            if (arcs[l][2*i]->target == -1 && arcs[l][2*i+1]->target == -1){
                remove_node(l, i);
                continue;
            }
            if (l > 0  && nodes[l][i]->in_arcs.empty()) {
                remove_node(l, i);
                continue;
            }

            update_node_bottom_up(l, i);
        }
    }

    if (debug){
        cout << "\n== End bottom-up procedure "<< round <<" ==\n";

        char name[256];
        sprintf(name, "graphs/bdd%d_bottom%d.gml", bdd_id, round);
        export_graph(name);
    }

}

void bdd_chance::update_node_top_down(int layer, int node_id) {

    bdd_node* node = nodes[layer][node_id];
    bdd_node* source = nullptr;
    bdd_arc* arc = nullptr;

    int lin_value = 0;
    int sigma_value = 0;

    node->lin_min_top = int_max;
    node->lin_max_top = int_min;
    node->sigma_min_top.assign(num_vars, int_max);
    node->sigma_max_top.assign(num_vars, int_min);

    assert(!node->in_arcs.empty());

    for (int i : node->in_arcs) {
        arc = arcs[layer - 1][i];
        source = nodes[layer - 1][arc->source];

        assert(!source->is_free);

        lin_value = arc->value*cqc->linear[layer - 1];

        // Update linear values
        if (lin_value + source->lin_min_top < node->lin_min_top)
            node->lin_min_top = lin_value  + source->lin_min_top ;

        if (lin_value + source->lin_max_top  > node->lin_max_top)
            node->lin_max_top = lin_value + source->lin_max_top ;

        // Update sigma matrix values
        for (int j = 0; j< num_vars; j++) {
            sigma_value = arc->value*cqc->sigma[j][layer - 1];

            if (sigma_value + source->sigma_min_top[j] < node->sigma_min_top[j])
                node->sigma_min_top[j] = sigma_value  + source->sigma_min_top[j] ;

            if (sigma_value + source->sigma_max_top[j]  > node->sigma_max_top[j])
                node->sigma_max_top[j] = sigma_value + source->sigma_max_top[j] ;
        }

    }

    // Check if node is exact or not
    node->is_exact = true;
    if (node->lin_min_top == node->lin_max_top) {
        for (int j = 0; j< num_vars; j++) {
            if (node->sigma_min_top[j] != node->sigma_max_top[j]) {
                node->is_exact = false;
                break;
            }
        }
    }
    else {
        node->is_exact = false;
    }

    //TODO: Check if all the solutions are feasible or not
    node->is_feasible = true;
    if (!node->is_reduce) {
        //check incoming arcs independently for feasibility
        for (int i : node->in_arcs) {
            if (!is_arc_feasible(arcs[layer-1][i])) {
                node->is_feasible = false;
                break;
            }
        }
    }


    if (debug){
        //Print node top-down information
        cout << "\nNode " << node_id << " in layer " << layer << " -- top-down" << endl;
        cout << "\tin_arcs  = ";
        for (int i : node->in_arcs) {
            cout << i << " ";
        }
        cout << "\n\tlinear   = [" << node->lin_min_top << " , " << node->lin_max_top << "]"<< endl;
        for (int j = 0; j< num_vars; j++) {
            cout << "\tsigma_" << j << " = [" << node->sigma_min_top[j] << " , " << node->sigma_max_top[j] << "]" << endl;
        }

        if (node->is_exact)
            cout <<"\tNode is exact" << endl;

        if (node->is_feasible)
            cout <<"\tNode has only feasible paths" << endl;
    }


}

void bdd_chance::update_node_bottom_up(int layer, int node_id) {

    bdd_node* node = nodes[layer][node_id];

    //Outgoing arcs
    bdd_arc* arc0 = arcs[layer][node_id*2];
    bdd_arc* arc1 = arcs[layer][node_id*2+1];

    if (arc0->target != -1) {
        node->lin_min_bottom = nodes[layer + 1][arc0->target]->lin_min_bottom;
        node->lin_max_bottom = nodes[layer + 1][arc0->target]->lin_max_bottom;

        for (int j = 0; j < num_vars; j++) {
            node->sigma_min_bottom[j] = nodes[layer + 1][arc0->target]->sigma_min_bottom[j];
            node->sigma_max_bottom[j] = nodes[layer + 1][arc0->target]->sigma_max_bottom[j];
        }

        if (arc1->target != -1) {
            node->lin_min_bottom = min(node->lin_min_bottom,
                                       nodes[layer + 1][arc1->target]->lin_min_bottom + cqc->linear[layer]);
            node->lin_max_bottom = max(node->lin_max_bottom,
                                       nodes[layer + 1][arc1->target]->lin_max_bottom + cqc->linear[layer]);

            for (int j = 0; j < num_vars; j++) {
                node->sigma_min_bottom[j] = min(node->sigma_min_bottom[j],
                                                nodes[layer + 1][arc0->target]->sigma_min_bottom[j] +
                                                cqc->sigma[j][layer]);
                node->sigma_max_bottom[j] = max(node->sigma_max_bottom[j],
                                                nodes[layer + 1][arc0->target]->sigma_max_bottom[j] +
                                                cqc->sigma[j][layer]);
            }
        }
    }
    else {
        node->lin_min_bottom = nodes[layer + 1][arc1->target]->lin_min_bottom + cqc->linear[layer];
        node->lin_max_bottom = nodes[layer + 1][arc1->target]->lin_max_bottom + cqc->linear[layer];

        for( int j = 0; j < num_vars; j++){
            node->sigma_min_bottom[j] = nodes[layer + 1][arc1->target]->sigma_min_bottom[j] + cqc->sigma[j][layer];
            node->sigma_max_bottom[j] = nodes[layer + 1][arc1->target]->sigma_max_bottom[j] + cqc->sigma[j][layer];
        }
    }

    if (debug){
        //Print node top-down information
        cout << "\nNode " << node_id << " in layer " << layer << " -- bottom-up" << endl;
        cout << "\tout_arcs  = ";
        if (arc0->target != -1) cout << node_id*2 << " ";
        if (arc1->target != -1) cout << node_id*2 + 1 << " ";

        cout << "\n\tlinear   = [" << node->lin_min_bottom << " , " << node->lin_max_bottom  << "]"<< endl;
        for (int j = 0; j< num_vars; j++) {
            cout << "\tsigma_" << j << " = [" << node->sigma_min_bottom[j] << " , " << node->sigma_max_bottom[j] << "]" << endl;
        }
    }
}


// Filter arcs
bool bdd_chance::is_arc_infeasible(bdd_arc* arc) {

    if (is_arc_infeasible_qcc(arc)){
        return true;
    }

    return false;
}

//Filter arc w.r.t to QCC constraint
bool bdd_chance::is_arc_infeasible_qcc(bdd_arc *arc){

    const int layer = arc->layer;

    int linear = nodes[layer][arc->source]->lin_min_top + arc->value*cqc->linear[layer]
                 + nodes[layer + 1][arc->target]->lin_min_bottom;

    int sigma = 0;

    // Compute minimum value of quadratic term
    compute_min_sigma(arc, sigma);

    if (linear + cqc->omega*sqrt(sigma) >= epsilon + cqc->rhs) {
        return true;
    }

    return false;
}

// Check if all paths over an arc are feasible or not
bool bdd_chance::is_arc_feasible(bdd_arc* arc){

    const int layer = arc->layer;

    int linear = nodes[layer][arc->source]->lin_max_top + arc->value*cqc->linear[layer]
                 + nodes[layer + 1][arc->target]->lin_max_bottom;

    int sigma = 0;

    // Compute maximum value of quadratic term
    compute_max_sigma(arc, sigma);

    if(linear + cqc->omega*sqrt(sigma) <= cqc->rhs) {
        return true;
    }

    return false;
}

void bdd_chance::compute_min_sigma(bdd_arc* arc, int& sigma) {

    const int layer = arc->layer;
    int aux_min = 0;
    int aux_max = 0;

    sigma = 0;

    // Compute minimum value of quadratic term
    for(int j = 0; j < num_vars; j++){
        aux_min = nodes[layer][arc->source]->sigma_min_top[j] + arc->value*cqc->sigma[j][layer]
                  + nodes[layer + 1][arc->target]->sigma_min_bottom[j];
        aux_max = nodes[layer][arc->source]->sigma_max_top[j] + arc->value*cqc->sigma[j][layer]
                  + nodes[layer + 1][arc->target]->sigma_max_bottom[j];

        if (aux_min >= 0)       sigma += aux_min*aux_min;
        else if (aux_max <= 0)  sigma += aux_max*aux_max;

    }

}

void bdd_chance::compute_max_sigma(bdd_arc* arc, int& sigma) {

    const int layer = arc->layer;
    int aux_min = 0;
    int aux_max = 0;

    sigma = 0;

    // Compute maximum value of quadratic term
    for(int j = 0; j < num_vars; j++){
        aux_min = nodes[layer][arc->source]->sigma_min_top[j] + arc->value*cqc->sigma[j][layer]
                  + nodes[layer + 1][arc->target]->sigma_min_bottom[j];
        aux_max = nodes[layer][arc->source]->sigma_max_top[j] + arc->value*cqc->sigma[j][layer]
                  + nodes[layer + 1][arc->target]->sigma_max_bottom[j];

        sigma += max(aux_min*aux_min, aux_max*aux_max);
    }

}


void bdd_chance::refine(int layer) {


    if (debug) cout << "\n= Refining layer " << layer << " =" << endl;

    // Get candidate arcs (lhs, linear, arc_id)
    int total = 0;
    int arc_id = 0;
    get_candidate_arcs(layer, total);

    int old_target_id = 0;
    list<int> in_arcs;

    // Iterate over arcs and create new nodes (util we reach max_width)
    for (int i = 0; i < total; i++) {
        if (get<0>(candidate_arcs[i]) == -1) continue;

        arc_id = get<2>(candidate_arcs[i]);
        old_target_id = arcs[layer - 1][arc_id]->target;

        //Skip arc if it is the last in the node
        if(nodes[layer][old_target_id]->in_arcs.size() == 1) continue;

        in_arcs.clear();
        in_arcs.push_back(get<2>(candidate_arcs[i]));

        remove_arc(arcs[layer - 1][arc_id]);

        //Get other arcs with similar value
        for (int j= i+1; j< total; j++) {

            // Check if arcs are equal
            if (get<0>(candidate_arcs[i]) == get<0>(candidate_arcs[j])) {
                if (are_arcs_equal(arcs[layer - 1][arc_id],
                                    arcs[layer - 1][get<2>(candidate_arcs[j])]) ) {

                    // Add equal arc to the same new node
                    in_arcs.push_back(get<2>(candidate_arcs[j]));
                    get<0>(candidate_arcs[j]) = -1;
                }
            }
            else {
                break;
            }
        }

        // Remove edges from old_target nodes
        for (int j : in_arcs) remove_arc(arcs[layer - 1][j]);

        // Create new node for the arcs
        create_new_node(layer,old_target_id, in_arcs);

        // Check availability
        if (nodes_available[layer] == 0) break;
    }

}

void bdd_chance::get_candidate_arcs(int layer, int& total) {
    int size = nodes[layer].size();
    double value = 0;
    int linear_value = 0;
    int sigma_value = 0;

    bdd_arc* arc = nullptr;

    candidate_arcs.assign(max_width*2, make_tuple(-1,-1,-1));
    total = 0;

    // Iterate over nodes in the layer
    for (int i = 0; i < size; i++) {
        if (nodes[layer][i]->is_free || nodes[layer][i]->is_exact || nodes[layer][i]->is_feasible) continue;

        // No refinement for nodes with one incoming arc
        if (nodes[layer][i]->in_arcs.size() == 1) continue;

        // Skip nodes that have all arcs equal (w.r.t. min value)
        if (are_in_arcs_min_equal(layer, i)) continue;

        // Iterate over their incoming arcs
        for (auto j: nodes[layer][i]->in_arcs) {
            arc = arcs[layer - 1][j];

            //Ignore arcs that are always feasible
            if (is_arc_feasible(arc)) continue;

            linear_value = arc->value*cqc->linear[layer - 1] + nodes[layer - 1][ arc->source]->lin_min_top;
            compute_min_sigma(arc, sigma_value);
            value = linear_value + cqc->omega*sqrt(sigma_value);

            //Add arc to the candidate list
            get<0>(candidate_arcs[total]) = value;
            get<1>(candidate_arcs[total]) = linear_value;
            get<2>(candidate_arcs[total]) = j;
            total++;
        }

    }

    //Sort the vector in decreasing order
    sort(candidate_arcs.begin(), candidate_arcs.begin() + total, greater<>() );

}

bool bdd_chance::are_in_arcs_min_equal(int layer, int id) {

    bdd_arc* first_arc = nullptr;

    for (int j: nodes[layer][id]->in_arcs) {
        if (first_arc == nullptr) {
            first_arc = arcs[layer - 1 ][j];
            continue;
        }

        if (!are_arcs_equal(first_arc, arcs[layer - 1 ][j] ) ) return false;
    }

    return true;
}

// Function to compare the state ( linear, sigma_1,sigma_2, ...) of two arcs  in the same layer
bool bdd_chance::are_arcs_equal(bdd_arc* arc1, bdd_arc* arc2){

    int aux1 = 0;
    int aux2 = 0;
    int layer = arc1->layer;

    assert(layer == arc2->layer);

    // Compare linear component
    aux1 = arc1->value*cqc->linear[layer] + nodes[layer][ arc1->source ]->lin_min_top;
    aux2 = arc2->value*cqc->linear[layer] + nodes[layer][ arc2->source ]->lin_min_top;

    if (aux1 != aux2) return false;

    // Compare sigma components
    for (int i = 0; i < num_vars; i++) {
        aux1 = arc1->value*cqc->sigma[i][layer] + nodes[layer][ arc1->source ]->sigma_min_top[i];
        aux2 = arc2->value*cqc->sigma[i][layer] + nodes[layer][ arc2->source ]->sigma_min_top[i];

        if (aux1 != aux2) return  false;
    }

    return true;
}


void bdd_chance::check_exact_layer(int layer){

    int size = nodes[layer].size();

    //Check if layer is exact or not
    if (!is_layer_exact[layer]) {

        is_layer_exact[layer] = true;

        for (int i = 0; i < size; i++) {
            if (nodes[layer][i]->is_free) continue;

            if (!nodes[layer][i]->is_exact) {
                is_layer_exact[layer] = false;
                break;
            }
        }
    }
}

void bdd_chance::create_new_node(int layer, int old_id, list<int>& in_arcs){
    assert(first_node_available[layer] >= 0);

    bdd_node* new_node = nodes[layer][first_node_available[layer]];

    //Activate new node
    new_node->is_free = false;
    new_node->in_arcs.assign( in_arcs.begin(), in_arcs.end());

    for (int i: in_arcs){
        arcs[layer - 1][i]->target = first_node_available[layer];
    }

    //Set up outgoing arcs
    int new_target = 0;

    for(int j =0; j < 2; j++){
        new_target = arcs[layer][2*old_id + j]->target;

        if (new_target >= 0) {
            arcs[layer][2*first_node_available[layer] + j]->target = new_target;
            nodes[layer + 1][new_target]->in_arcs.push_back(2*first_node_available[layer] + j);
        }
    }

    //Update available nodes
    update_first_node_available(layer);
    nodes_available[layer]--;

    bdd_has_change = true;
}

void bdd_chance::update_first_node_available(int layer){

    //Find new first node available
    for (int i = first_node_available[layer]; i < nodes[layer].size(); i++){
        if (nodes[layer][i]->is_free){
            first_node_available[layer] = i;
            return;
        }
    }

    first_node_available[layer] = -1;

    assert(nodes_available[layer] <= 1);
}

void bdd_chance::reduce_bdd() {
    //Procedure to reduce the final BDD

    for (int l = num_vars-1; l >= 0; l--) {
        //Reduce nodes in the layer
        reduce_layer(l, true);
    }

    if (debug){
        cout << "\n== End final BDD reduction  ==\n";

        char name[256];
        sprintf(name, "graphs/bdd%d_reduced.gml", bdd_id);
        export_graph(name);
    }

}

void bdd_chance::reduce_layer(int layer, bool final_bdd) {

    int total = 0;

    if(final_bdd)   total = get_nodes_to_reduce_final(layer);
    else            total = get_nodes_to_reduce(layer);


    if (total == 0) return;

    list<int> similar;

    //Iterate over the list add
    for (int i =0; i < total; i++) {
        if(candidate_nodes[i] == -1) continue;

        similar.clear();

        //Find similar nodes
        for (int j = i+1; j < total; j++) {
            if (candidate_nodes[j] == -1) continue;

            if (arcs[layer][2*candidate_nodes[i] ]->target == arcs[layer][2*candidate_nodes[j] ]->target &&
                arcs[layer][2*candidate_nodes[i] +1 ]->target == arcs[layer][2*candidate_nodes[j] +1 ]->target) {
                similar.push_back(candidate_nodes[j]);
                candidate_nodes[j] = - 1;
            }
        }

        //Reduce nodes
        if (!similar.empty()) reduce_nodes(layer,candidate_nodes[i], similar);
    }

}

int bdd_chance::get_nodes_to_reduce(int layer) {
    //Get list of possible nodes to reduce
    //i.e., exact/feasible nodes that have exact/feasible target nodes (or terminal node)

    int size = nodes[layer].size();
    int total = 0;

    candidate_nodes.assign(max_width, -1);

    for (int i = 0; i < size; i++) {
        if (nodes[layer][i]->is_free) continue;

        //Remove nodes that are disconnect
        if (arcs[layer][2*i]->target == -1 && arcs[layer][2*i+1]->target == -1){
            remove_node(layer, i);
            continue;
        }

        if (layer > 0  && nodes[layer][i]->in_arcs.empty()) {
            remove_node(layer, i);
            continue;
        }

        // Check that node is either exact or feasible
        if (!nodes[layer][i]->is_exact && !nodes[layer][i]->is_feasible) continue;

        // Check that node have exact or feasible targets
        if (layer == num_vars - 1) {
            candidate_nodes[total] = i;
            total++;
        }
        else if (arcs[layer][2*i]->target == -1 ||
                 (nodes[layer+1][arcs[layer][2*i]->target]->is_exact || nodes[layer+1][arcs[layer][2*i]->target]->is_feasible) ){
            if(arcs[layer][2*i+1]->target == -1 ||
               (nodes[layer+1][arcs[layer][2*i+1]->target]->is_exact || nodes[layer+1][arcs[layer][2*i+1]->target]->is_feasible)){
                candidate_nodes[total] = i;
                total++;
            }
        }
    }

    return total;
}

int bdd_chance::get_nodes_to_reduce_final(int layer){
    //Get list of possible nodes to reduce
    //In this case we don't care about the status of the target nodes


    int size = nodes[layer].size();
    int total = 0;

    candidate_nodes.assign(max_width, -1);

    for (int i = 0; i < size; i++) {
        if (nodes[layer][i]->is_free) continue;

        //Remove nodes that are disconnect
        if (arcs[layer][2*i]->target == -1 && arcs[layer][2*i+1]->target == -1){
            remove_node(layer, i);
            continue;
        }

        if (layer > 0  && nodes[layer][i]->in_arcs.empty()) {
            remove_node(layer, i);
            continue;
        }

        candidate_nodes[total] = i;
        total++;

    }

    return total;

}

void bdd_chance::reduce_nodes(int layer, int node_id, list<int>& similar) {

    bdd_node* node = nodes[layer][node_id];

    for (auto i: similar) {
        // Redirect arcs
        for (auto arc_id : nodes[layer][i]->in_arcs) {
            node->in_arcs.push_back(arc_id);
            arcs[layer-1][arc_id]->target = node_id;
        }

        // Remove node
        nodes[layer][i]->in_arcs.clear();
        remove_node(layer, i);
    }

    bdd_has_change = true;

    node->is_reduce = true;
    node->is_feasible = true;
    node->is_exact = false;

    if (debug){
        cout << "Merge nodes: ";
        for (auto i: similar) cout << i << " ";
        cout << "in layer "<< layer <<" into node " << node_id << endl;
    }
}



void bdd_chance::remove_arc(bdd_arc* arc) {

    if (arc->target == -1) return;

    //Remove arc from target node
    list<int>::iterator it = nodes[arc->layer + 1][arc->target]->in_arcs.begin();
    list<int>::iterator end = nodes[arc->layer + 1][arc->target]->in_arcs.end();
    int id = 2*arc->source + arc->value;

    while( it != end){
        if( *it == id ){
            it = nodes[arc->layer + 1][arc->target]->in_arcs.erase(it);
            break;
        }
        it++;
    }

    // Desactivate arc
    arc->target = -1;

    bdd_has_change = true;
}

void bdd_chance::remove_node(int layer, int node_id) {

    if (debug) {
        cout << "Removing node "<< node_id << " in layer " << layer << endl;
    }

    bdd_node* node = nodes[layer][node_id];

    //Reset node values
    node->is_feasible = false;
    node->is_free = true;
    node->is_exact = false;
    node->is_reduce = false;

    node->lin_min_top = 0;
    node->lin_max_top = 0;
    node->sigma_min_top.assign(num_vars, 0);
    node->sigma_max_top.assign(num_vars, 0);

    node->lin_min_bottom = int_min;
    node->lin_max_bottom = int_max;
    node->sigma_min_bottom.assign(num_vars, int_min);
    node->sigma_max_bottom.assign(num_vars, int_max);

    // Reset incoming arcs (if any)
    for (int i : node->in_arcs) {
        arcs[layer - 1][i]->target = -1;
    }

    node->in_arcs.clear();

    // Reset outgoing arcs
    remove_arc(arcs[layer][node_id*2]);
    remove_arc(arcs[layer][node_id*2 +1 ]);

    //Update available nodes
    nodes_available[layer]++;
    if (node_id < first_node_available[layer] || first_node_available[layer] == -1) {
        first_node_available[layer] = node_id;
    }

    bdd_has_change = true;
}

// =============== MIN-CUT SEPARATION FUNCTIONS - BEGIN

bool bdd_chance::separate_min_cut(double* _cost, double* pi, double& pi0) {

    //Get active nodes
    get_active_nodes();

    // Set values
    cost = _cost;
    set_values_max_flow();

    if(debug){
        cout << "Capacity BDD "<< bdd_id << ": ";
        for(int i = 0; i < cqc->vars.size(); i++){
            cout << cost[i] << " ";
        }
        cout << endl;
    }

    // Compute max-flow
    double max_flow = 0;

    compute_max_flow(max_flow);

    if (max_flow >= 1 - epsilon) {

        return false;
    }

    // Get min-cut
    get_min_cut_arcs(pi, pi0);

    // Double check that point violates constraint (to avoid precision errors)
    return  is_constraint_violated(pi, pi0);

}

void bdd_chance::get_active_nodes(){

    if (!active_nodes.empty()) return;

    active_nodes.resize(num_vars + 1, vector<int>() );

    int n = 0;
    int count = 0;

    for (int l = 0; l < num_vars + 1; l++ ) {
        n = nodes[l].size();
        count = 0;

        active_nodes[l].resize(nodes[l].size(), -1);

        for (int i = 0; i < n; i++) {
            if (nodes[l][i]->is_free) continue;

            active_nodes[l][count] = i;
            count++;
        }

        active_nodes[l].resize(count);
    }

}

void bdd_chance::set_values_max_flow() {

    // Iterate over nodes and outgoing arcs
    int num_nodes = 0;
    int node_id = 0;
    bdd_node* node = nullptr;

    for (int layer = 0; layer <= num_vars; layer++) {
        num_nodes = active_nodes[layer].size();

        for (int i = 0; i < num_nodes; i++) {
            node_id = active_nodes[layer][i];
            node = nodes[layer][node_id];

            node-> visited = false;
            node->parent.first = -1;
            node->parent.second = -1;

            if (layer == num_vars) break;

            // Update outgoing arcs
            for (int j = 0; j < 2; j++) {
                if (arcs[layer][2*node_id +j]->target == -1) continue;

                if (j == 0) arcs[layer][2*node_id +j]->res_value = 1 - cost[layer];
                else        arcs[layer][2*node_id +j]->res_value = cost[layer];

                arcs[layer][2*node_id +j]->res_rev_value = 0.0;
            }
        }
    }
}

//Ford-Fulkerson algorithm for max_flow
void bdd_chance::compute_max_flow(double& max_flow) {

    double path_flow = 0;
    double edge_flow = 0;

    int layer = -1;
    int node_id = -1;
    bool found_path = true;
    bool found_edge = true;

    bdd_node* node = nullptr;

    max_flow = 0;

    while(found_path) {
        // Find augmentation path
        found_path = find_augmenting_path();

        if(!found_path) break;

        //Get flow in the path
        path_flow = 2;
        node = nodes[num_vars][0]; //terminal node
        layer = num_vars;
        node_id = 0;

        while ( !(layer == 0 && node_id == 0) ) {
            found_edge = false;
            if (node->parent.first == layer -1) {
                //First cases: layer of parent is previous layer
                // Check which outgoing arc of the parent points to the node_id
                for (int j =0; j < 2; j++) {
                    edge_flow = arcs[layer - 1][2*node->parent.second +j]->res_value;
                    if (arcs[layer - 1][2*node->parent.second +j]->target == node_id && edge_flow >= epsilon ) {
                        path_flow = min(path_flow, edge_flow);
                        found_edge = true;
                        break;
                    }
                }
            }
            else {
                //Second case: layer of parent is next layer
                assert(node->parent.first == layer + 1);

                // Check which outgoing arc of the parent (new_node) points to new node_id
                for (int j =0; j < 2; j++) {
                    edge_flow = arcs[layer][2*node_id +j]->res_rev_value;
                    if (arcs[layer][2*node_id +j]->target == node->parent.second && edge_flow >= epsilon) {
                        path_flow = min(path_flow, edge_flow);
                        found_edge = true;
                        break;
                    }
                }
            }

            assert(found_edge);

            //Update values of next node
            layer = node->parent.first;
            node_id = node->parent.second;
            node= nodes[layer][node_id];
        }

        //Update residual values
        node = nodes[num_vars][0]; //terminal node
        layer = num_vars;

        while ( !(layer == 0 && node_id == 0) ) {
            if (node->parent.first == layer -1) {
                //First cases: layer of parent is previous layer
                // Check which outgoing arc of the parent points to the node_id
                for (int j =0; j < 2; j++) {
                    edge_flow = arcs[layer - 1][2*node->parent.second +j]->res_value;
                    if (arcs[layer - 1][2*node->parent.second +j]->target == node_id && edge_flow >= epsilon ) {
                        arcs[layer - 1][2*node->parent.second +j]->res_value -= path_flow;
                        arcs[layer - 1][2*node->parent.second +j]->res_rev_value += path_flow;
                        break;
                    }
                }
            }
            else {
                //Second case: layer of parent is next layer
                // Check which outgoing arc of the current node that points to new  node_id
                for (int j =0; j < 2; j++) {
                    edge_flow = arcs[layer][2*node_id +j]->res_rev_value;
                    if (arcs[layer][2*node_id +j]->target == node->parent.second && edge_flow >= epsilon) {
                        arcs[layer][2*node_id +j]->res_rev_value -= path_flow;
                        arcs[layer][2*node_id +j]->res_value += path_flow;
                        break;
                    }
                }
            }

            //Update values of next node
            layer = node->parent.first;
            node_id = node->parent.second;
            node= nodes[layer][node_id];

        }

        //Update max_flow
        max_flow += path_flow;
    }

    if (debug) {
        cout << "\t max_flow = "<<  max_flow << endl;
    }

}

bool bdd_chance::find_augmenting_path() {

    // Reset visited nodes
    reset_visited_nodes();

    queue <pair<int,int> > q_nodes; //(layer, node_id);
    q_nodes.push(make_pair(0,0));

    nodes[0][0]->visited = true;

    int layer = -1;
    int node_id = -1;

    int target = -1;
    int source = -1;

    bool found_terminal = false;

    while(!q_nodes.empty()) {
        layer = q_nodes.front().first;
        node_id = q_nodes.front().second;

        q_nodes.pop();

        // Iterate over outgoing arcs
        for (int j = 0; j< 2; j++) {
            target = arcs[layer][node_id*2 + j]->target;
            if (target == -1) continue;
            if (arcs[layer][node_id*2 + j]->res_value <= epsilon) continue; // residual value is "0"
            if (nodes[layer + 1][target]->visited) continue; //target already visited

            //Add node to the queue
            nodes[layer + 1][target]->visited = true;
            nodes[layer + 1][target]->parent.first = layer;
            nodes[layer + 1][target]->parent.second = node_id;

            q_nodes.push(make_pair(layer + 1, target));

            if (layer == num_vars - 1){
                found_terminal = true;
                break;
            }
        }

        if (found_terminal) break;

        // Iterate over incoming arcs (for reversible edges)
        for (int j : nodes[layer][node_id]->in_arcs) {
            if (arcs[layer - 1][j]->res_rev_value <= epsilon) continue; // residual value is "0"

            source = arcs[layer - 1][j]->source;

            if (nodes[layer - 1][source]->visited) continue; //source already visited

            //Add node to the queue
            nodes[layer - 1][source]->visited = true;
            nodes[layer - 1][source]->parent.first = layer;
            nodes[layer - 1][source]->parent.second = node_id;

            q_nodes.push(make_pair(layer - 1, source));
        }

    }

    if (found_terminal) return true;

    return false;
}

void bdd_chance::reset_visited_nodes() {

    int num_nodes = 0;
    int node_id = 0;
    bdd_node* node = nullptr;

    for (int layer = 0; layer <= num_vars; layer++) {
        num_nodes = active_nodes[layer].size();

        for (int i = 0; i < num_nodes; i++) {
            node_id = active_nodes[layer][i];
            node = nodes[layer][node_id];

            node->visited = false;
            node->parent.first = -1;
            node->parent.second = -1;
        }
    }

}

void bdd_chance::get_min_cut_arcs(double* pi, double& pi0) {

    // Reset visited nodes
    reset_visited_nodes();

    // Compute Visited nodes in residual graph
    compute_visited_nodes();

    //Get min-cut edges in original graph
    int num_nodes = 0;
    bdd_node* node = nullptr;
    double min_cut = 0;

    pi0 = -1;

    for (int layer = 0; layer < num_vars; layer++) {
        num_nodes = nodes[layer].size();

        for (int i = 0; i < num_nodes; i++) {
            node = nodes[layer][i];
            if (node->is_free || !node->visited) continue;

            //Iterate over outgoing arcs
            for (int j=0; j < 2; j++) {
                if (arcs[layer][2*i+j]->target == -1) continue;
                if (nodes[layer + 1][arcs[layer][2*i+j]->target]->visited ) continue;

                // Found arc that is part of the min_cut!
                if (j == 0) {
                    min_cut += 1 - cost[layer];
                    pi[layer] += 1;
                    pi0 +=1;
                }
                else   {
                    min_cut += cost[layer];
                    pi[layer] -= 1;
                }
            }
        }
    }

    if (debug) {
        cout << "\t min_cut  = "<<  min_cut << endl;
    }
}

void bdd_chance::compute_visited_nodes() {
    queue <pair<int,int> > q_nodes; //(layer, node_id);
    q_nodes.push(make_pair(0,0));

    nodes[0][0]->visited = true;

    int layer = -1;
    int node_id = -1;

    int target = -1;
    int source = -1;

    while(!q_nodes.empty()) {
        layer = q_nodes.front().first;
        node_id = q_nodes.front().second;

        q_nodes.pop();

        // Iterate over outgoing arcs
        for (int j = 0; j< 2; j++) {
            target = arcs[layer][node_id*2 + j]->target;
            if (target == -1) continue;
            if (arcs[layer][node_id*2 + j]->res_value <= epsilon) continue; // residual value is "0"
            if (nodes[layer + 1][target]->visited) continue; //target already visited

            //Add node to the queue
            nodes[layer + 1][target]->visited = true;

            q_nodes.push(make_pair(layer + 1, target));
        }

        // Iterate over incoming arcs (for reversible edges)
        for (int j : nodes[layer][node_id]->in_arcs) {
            if (arcs[layer - 1][j]->res_rev_value <= epsilon) continue; // residual value is "0"

            source = arcs[layer - 1][j]->source;

            if (nodes[layer - 1][source]->visited) continue; //source already visited

            //Add node to the queue
            nodes[layer - 1][source]->visited = true;

            q_nodes.push(make_pair(layer - 1, source));
        }
    }
}

bool bdd_chance::is_constraint_violated(double* pi, double& pi0) {

    double lhs = 0;

    for (int i = 0; i < num_vars; i++) {
        lhs +=pi[i]*cost[i];
    }

    return (lhs - pi0 >= epsilon);
}


// =============== MIN-CUT SEPARATION FUNCTIONS - END

// =============== TARGET CUTS FUNCTIONS - BEGIN
void bdd_chance::compute_center_point() {

    if (overflow) return;

    //Initialize structure
    count_paths_top.resize(num_vars + 1, vector<long long>() );
    count_paths_bottom.resize(num_vars + 1, vector<long long>() );

    for (int l = 0; l < num_vars + 1; l++ ) {
        count_paths_top[l].resize(nodes[l].size(), 0);
        count_paths_bottom[l].resize(nodes[l].size(), 0);
    }


    //Count paths top & bottom
    count_path_top();
    if (overflow) return;
    count_path_bottom();

    //Calculate center point
    center_point.resize(num_vars, 0);

    int num_nodes = 0;
    double num_zero = 0;
    double num_one = 0;
    int arc_id = 0;
    bdd_node* node = nullptr;

    for (int l = 0; l < num_vars; l ++) {
        num_nodes = nodes[l].size();
        num_zero = 0;
        num_one = 0;

        //Count number of zero and one paths in each layer
        for (int i = 0; i < num_nodes; i++) {
            if (nodes[l][i]->is_free) continue;
            node = nodes[l][i];

            //Iterate over outgoing arcs
            for (int j = 0; j< 2; j++) {
                arc_id = 2*i + j;
                if(arcs[l][arc_id]->target == -1) continue;

                if (j == 0) {
                    num_zero += count_paths_top[l][i] * count_paths_bottom[l+1][arcs[l][arc_id]->target];
                }
                else{
                    num_one += count_paths_top[l][i] * count_paths_bottom[l+1][arcs[l][arc_id]->target];
                }

            }
        }
        center_point[l] = num_one/count_paths_bottom[0][0];
    }
}

void bdd_chance::count_path_top() {

    int num_nodes = 0;
    bdd_node* node = nullptr;

    count_paths_top[0][0] = 1;

    for (int l = 1; l < num_vars + 1; l++ ) {
        num_nodes = nodes[l].size();

        for (int i = 0; i < num_nodes; i++) {
            node = nodes[l][i];

            if (node->is_free) continue;

            //Iterate over incoming arcs
            for (int arc_id : node->in_arcs) {
                count_paths_top[l][i] += count_paths_top[l-1][ arcs[l-1][arc_id]->source ];

                if (count_paths_top[l][i] < 0 || count_paths_top[l][i] < count_paths_top[l-1][ arcs[l-1][arc_id]->source ] ) {
                    cout << " Long overflow - cannot compute center point" << endl;
                    overflow = true;
                    return;
                }
            }
        }
    }
}

void bdd_chance::count_path_bottom() {

    int arc_id = 0;
    int num_nodes = 0;
    bdd_node* node = nullptr;

    count_paths_bottom[num_vars][0] = 1;

    for (int l = num_vars - 1; l >= 0; l--) {
        num_nodes = nodes[l].size();

        for (int i = 0; i < num_nodes; i++) {
            node = nodes[l][i];

            if (node->is_free) continue;

            //Iterate over outgoing arcs
            for (int j = 0; j< 2; j++) {
                arc_id = 2*i + j;

                if(arcs[l][arc_id]->target == -1) continue;

                count_paths_bottom[l][i] += count_paths_bottom[l+1][ arcs[l][arc_id]->target ];
            }
        }
    }

}

// =============== TARGET CUTS FUNCTIONS - ENDS

// =============== SLACK FUNCTIONS - BEGIN

bool bdd_chance::get_bdd_slacks(double* _cost, double* _slack, double& rhs){

    if (!is_constructed){
        cout << "Error: Need to construct BDD before getting slacks. Use method construct_bdd()" << endl;
        exit(1);
    }

    if (debug) {
        cout << "= Computing BDD slacks =" << endl;
    }

    // Set values
    cost = _cost;
    slack = _slack;

    // Compute top-down cost
    compute_cost_top();

    //Check if rhs is correct or not
    if (rhs > nodes[num_vars][0]->cost_top) {
        if(debug) cout << "Improve rhs of the inequality form " << rhs << " to " <<nodes[num_vars][0]->cost_top <<endl;

        rhs = nodes[num_vars][0]->cost_top;
    }
    else if(rhs < nodes[num_vars][0]->cost_top) {
        if(debug) cout << "Can't improve inequality with BDD " << bdd_id << endl;

        return false;
    }

    // Compute bottom-up cost
    compute_cost_bottom();

    // Compute slacks
    compute_slacks();


    if (debug){
        cout << "slacks = ";
        for (int i = 0; i < num_vars; i++) {
            cout << slack[i] << " ";
        }
        cout << "\n";
    }


    return true;
}

void bdd_chance::compute_cost_top() {

//    if (debug) {
//        cout << "\n== Computing BDD cost - top ==" <<endl;
//    }

    int size = 0;

    // Iterate over all node layers
    for (int l = 1; l < num_vars + 1; l++) {
        size = nodes[l].size();

        for (int i = 0; i < size; i++){
            if (nodes[l][i]->is_free) continue;

            compute_cost_node_top(l, i);
        }
    }
}

void bdd_chance::compute_cost_bottom() {

//    if (debug) {
//        cout << "\n== Computing BDD cost - bottom ==" <<endl;
//    }

    int size = 0;

    // Iterate over all node layers
    for (int l = num_vars - 1; l >= 0; l--) {
        size = nodes[l].size();

        for (int i = 0; i < size; i++){
            if (nodes[l][i]->is_free) continue;

            compute_cost_node_bottom(l, i);
        }
    }
}

void bdd_chance::compute_slacks() {

    int size = 0;
    double value = 0;

    // Iterate over all the edge layers
    for (int l = 0; l < num_vars; l++) {
        size = arcs[l].size();

        //Iterate over 0-edges (even) first and then 1-edges (odd)
        for (int j = 0; j < 2; j++) {

            cost_edge[j] = int_min;

            for (int i = j; i < size; i += 2) {
                if (arcs[l][i]->target == -1) continue;

                value = nodes[l][arcs[l][i]->source]->cost_top + j*cost[l]
                        + nodes[l+1][arcs[l][i]->target]->cost_bottom;

//                cout << "edge["<<l<<"]["<<i<<"] = " << value << endl;

                if (cost_edge[j] < value)
                    cost_edge[j] = value;
            }
//            cout << "- cost_edge["<<l<<"]["<<j<<"] = " << cost_edge[j] <<endl;
        }

        if (cost_edge[0] == int_min || cost_edge[1] == int_min){
            slack[l] = 0;
        }
        else{
            slack[l] = cost_edge[0] - cost_edge[1];
        }


//        cout << "- slack["<<l<<"] = " << slack[l] << endl;
    }


}

void bdd_chance::compute_cost_node_top(int layer, int id) {

    bdd_node* node = nodes[layer][id];
    double value = 0;

    node->cost_top = int_min;

    // Iterate over incoming arcs
    for (int i : node->in_arcs) {
        value = nodes[layer - 1][arcs[layer - 1][i]->source]->cost_top
                + arcs[layer - 1][i]->value*cost[layer - 1];

        if (node->cost_top < value){
            node->cost_top = value;
            node->cost_arc_id  = i;
        }
    }

//    if (debug) {
//        cout << "\nNode " << id << " in layer " << layer << " -- top-down" << endl;
//        cout << "\tin_arcs  = ";
//        for (int i : node->in_arcs) {
//            cout << i << " ";
//        }
//        cout<<"\n\tcost = " << node->cost_top << endl;
//        cout<<"\tarc = " << node->cost_arc_id << endl;
//    }
}

void bdd_chance::compute_cost_node_bottom(int layer, int id) {

    bdd_node* node = nodes[layer][id];
    double value = 0;

    node->cost_bottom = int_min;

    // Iterate over outgoing arcs
    for (int i = 0; i < 2; i++) {
        if (arcs[layer][2*id + i]->target == -1) continue;

        value = nodes[layer +1][arcs[layer][2*id + i]->target]->cost_bottom
                + arcs[layer][2*id + i]->value*cost[layer];

        if (node->cost_bottom < value) {
            node->cost_bottom = value;
        }
    }

//    if (debug){
//        cout << "\nNode " << id << " in layer " << layer << " -- bottom-up" << endl;
//        cout << "\tout_arcs  = ";
//        for (int i = 0; i < 2; i++) {
//            if (arcs[layer][2*id + i]->target == -1) continue;
//            cout << 2*id + i << " ";
//        }
//        cout<<"\n\tcost = " << node->cost_bottom << endl;
//    }

}



// =============== SLACK FUNCTIONS - END

// =============== LONGEST PATH - BEGIN

void bdd_chance::compute_longest_path(double* _cost, int* path) {

    if (!is_constructed){
        cout << "Error: Need to construct BDD before getting slacks. Use method construct_bdd()" << endl;
        exit(1);
    }

    if (debug) {
        cout << "= Computing Shortest path =" << endl;
    }

    // Set values
    cost = _cost;

    // Compute top-down cost
    compute_cost_top();

    // Extract Longest path
    extract_longest_path(path);
}

void bdd_chance::extract_longest_path(int* path) {

    bdd_arc* arc = nullptr;
    bdd_node* current_node = nodes[num_vars][0]; // start with target node

    for (int l = num_vars; l > 0; l--) {
        assert(current_node->cost_arc_id >= 0);

        arc = arcs[l-1][current_node->cost_arc_id];
        path[l-1] = arc->value;

        current_node = nodes[l-1][arc->source];
    }

    if (debug){
        cout << "Longest path = ";
        for (int i =0 ; i < num_vars; i++) {
            cout << path[i] << " ";
        }
        cout << endl;
    }

}

// =============== LONGEST PATH - END


// Export BDD graph
void bdd_chance::export_graph(const char* filename) {

    cout << "Exporting graph in gml format " << endl;
    cout << "\t max_width = " << max_width << "\t,\tmax_layer = "<< num_vars + 1 << endl;

    //create file
    ofstream bdd_file(filename);
    bdd_file << "graph [\n";
    bdd_file << "\t directed 1\n";
    bdd_file << "\t hierarchic 1\n";

    //Print nodes
    for (int l = 0; l < num_vars + 1; ++l) {
        for (int i = 0; i < nodes[l].size(); i++) {
            if (nodes[l][i]->is_free) continue;

            bdd_file << "node [ \n";
            bdd_file << "\t id " << l * max_width + i << "\n";
            bdd_file << "\t label \"" << l << "," << i << "\n";
            bdd_file << "\n\tlin = [" << nodes[l][i]->lin_min_top << " , " << nodes[l][i]->lin_max_top << "]"<< endl;
            for (int j = 0; j < num_vars; j++) {
                bdd_file << "\ts_"<<j<<" = [" << nodes[l][i]->sigma_min_top[j] << " , " << nodes[l][i]->sigma_max_top[j] << "]"<< endl;
            }
            bdd_file << "\"\n";
            bdd_file << "\t graphics [ \n";
            bdd_file << "\t\t type \"rectangle\"\n";
            if (nodes[l][i]->is_exact) {
                bdd_file << "\t\t fill \"#cce5ff\" \n"; // light blue -> exact node
            }
            else if(nodes[l][i]->is_feasible) {
                bdd_file << "\t\t fill \"#d8ffcc\" \n"; // light green -> feasible node
            }
            else{
                bdd_file << "\t\t hasFill 0 \n";
            }
            bdd_file << "\t\tw 110.0";
            bdd_file << "\t\th 120.0";
            bdd_file << "\t\t ] \n";
            bdd_file << "\t ]\n";

        }
    }

    //Print Edges
    for (int l = 0; l < num_vars; ++l) {
        for (int i = 0; i < arcs[l].size(); i++) {
            if(arcs[l][i]->target == -1) continue;
            bdd_file << "edge [ \n";
            bdd_file << "\t source " << l*max_width + arcs[l][i]->source << "\n";
            bdd_file << "\t target " << (l+1)*max_width + arcs[l][i]->target << "\n";
            bdd_file << "\t label \"" << arcs[l][i]->value;
            if (arcs[l][i]->value == 1)
                bdd_file<< " (" << cqc->vars[l] <<") ";

            bdd_file << "\"\n";
            bdd_file << "\t graphics [ \n";
            bdd_file << "\t \t fill \"#000000\" ";    //give color to an edge
            bdd_file << "\t\ttargetArrow \"diamond\"";
            bdd_file << "\t \t ]\n";
            bdd_file << "\t ]\n";
        }
    }
}


bdd_chance::~bdd_chance() {

    for (int l = 0; l < num_vars; ++l) {
        for (auto & i : arcs[l]) {
            delete i;
        }
    }

    for (int l = 0; l < num_vars + 1; ++l) {
        for (auto & i : nodes[l]) {
            delete i;
        }
    }

}