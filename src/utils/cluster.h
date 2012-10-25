#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "main/node.h"
#include "utils/model_factory.h"
#include "utils/settings_handle.h"
#include "utils/log_output.h"

using namespace std;
using namespace ppa;

class Cluster {
    Node* root;
    Node* prev_root;
    vector<int> cluster_members;
    double id_threshold;
    bool fast_cluster;
    
    double calc_identity(Node* node, Node* read);
    void finalise();
    
 public:
    Cluster(Node* n, int read_index) :
        root(n),
        prev_root(NULL),
        cluster_members(),
        id_threshold(Settings_handle::st.get("cluster-identity").as<float>()),
        fast_cluster(Settings_handle::st.is("cluster-fast")) {
    
        cluster_members.push_back(read_index);
    }
    
    ~Cluster() {
        if(root)
            delete root;
        if(prev_root)
            delete prev_root;
    }
    
    double add_to_cluster(Node* n, Model_factory *mf, int read_index);
    void remove_previous_read();
    
    Node* get_root() {
        return root;
    }
    
    vector<int>* get_members() {
        return &cluster_members;
    }
};

#endif
