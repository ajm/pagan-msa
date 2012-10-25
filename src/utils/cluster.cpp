
#include "utils/cluster.h"

double Cluster::add_to_cluster(Node* new_read, Model_factory *mf, int read_index) {
    
    finalise();
    
    stringstream ss;
    ss << "#" << read_index << "#";
    
    Node* new_root = new Node();
    new_root->set_name(ss.str());
    new_root->add_left_child(root);
    new_root->add_right_child(new_read);
    
    new_root->set_nhx_tid(root->get_nhx_tid());    
    new_read->set_nhx_tid(root->get_nhx_tid());
    
    root->set_distance_to_parent(0.001);
    
    
    new_root->align_sequences_this_node(mf, true, false);
    
    
    double read_identity = calc_identity(new_root, new_read);
    bool cluster_success = read_identity > id_threshold;
    
    new_root->has_left_child(false);
    new_root->has_right_child(false);
    
    
    if(cluster_success) {
        // clustering was successful, so add the read to list of members
        // and replace the old root with the new one
        if(fast_cluster) {
            delete root;
        }
        else {
            prev_root = root;
        }
    
        root = new_root;
        cluster_members.push_back(read_index);
    }
    else {
        // not added to cluster, so ignore old root and read and
        // delete the failed new root
        delete new_root;
    }
    
    return read_identity;
}

// this should not be called when user selects
// --cluster-fast, but then prev_root will always be NULL
void Cluster::remove_previous_read() {
    if(prev_root != NULL) {
        delete root;
        root = prev_root;
        prev_root = NULL;
        cluster_members.pop_back();
    }
}

void Cluster::finalise() {
    if(prev_root) {
        delete prev_root;
        prev_root = NULL;
    }
}

double Cluster::calc_identity(Node* node, Node* read) {
    Sequence* node_sequence = node->get_sequence();
    Sequence* read_sequence = read->get_sequence();
    Sequence* prev_sequence = root->get_sequence();
    
    int read_length = 0;
    int ref_length = 0;
    int aligned = 0;
    int total_sites = 0;
    
    int last_read_index = read_sequence->sites_length() - 1;
    bool in_read = false;
    
    for(int j = 0; j < node_sequence->sites_length(); ++j) {
        Site_children* offspring = node_sequence->get_site_at(j)->get_children();
        
        if(offspring->right_index == 1) {
            in_read = true;
        }
        
        if(offspring->right_index == last_read_index) {
            in_read = false;
            break;
        }
        
        if(not in_read) {
            continue;
        }
        
        bool read_has_site = offspring->right_index != -1;
        bool any_other_has_site = offspring->left_index != -1;
        
        if(read_has_site) {
            read_length++;
            
            if(any_other_has_site) {
                aligned++;
            }
        }
        
        ref_length++;
        
        if(any_other_has_site) {
            int number_sites = prev_sequence->get_site_at(offspring->left_index)->get_sumDNA();
            
            if(number_sites == 0) {
                number_sites = 1;
            }
            
            total_sites += number_sites;
        }
    }
    
    double average_length = total_sites / double(cluster_members.size());
    
    
    // aligned_proportion_cluster prevents deletions
    // aligned_proportion_read prevents insertions
    
    //fprintf(stderr, "aligned=%d read_length=%d average_length=%.3f prop_read=%.3f prop_clust=%.3f\n",
    //                aligned, read_length, average_length, aligned_proportion_read, aligned_proportion_cluster);
    
    //return (aligned / float(read_length)) > min_overlap;
    //return (aligned / min(float(read_length), average_length)) > min_overlap;
    
    // version 3 (identical to version 1, duh)
    return aligned / max(double(read_length), average_length);
}
