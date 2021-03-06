/***************************************************************************
 *   Copyright (C) 2010 by Ari Loytynoja                                   *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef READS_ALIGNER_H
#define READS_ALIGNER_H

#include <vector>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/model_factory.h"
#include "utils/fasta_entry.h"
#include "utils/fasta_reader.h"
#include "main/node.h"
#include "main/sequence.h"
#include "main/viterbi_alignment.h"

using namespace std;

namespace ppa
{


class Reads_aligner
{
    Node *global_root;
    map<string,string> codon_to_aa;
    map<string,string> aa_to_codon;

    void loop_default_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    void loop_simple_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    void loop_two_strand_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    void loop_translated_placement(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);

    void create_clusters(Node *root, vector<Fasta_entry> *reads, Model_factory *mf, int count);
    string create_edit_string(Node* node);
    void get_lowmem_alignment(vector<Fasta_entry>* aligned_sequences, vector<Fasta_entry>* reads, vector<int>* read_order);
    void tokenise(vector<Edit_operation>& tokens, string& editstr);
    void find_and_pad(string& mask, string& out, int& index);
    bool non_recursive_alignment_overlap(Node* node, string read_name, string ref_node_name);

    void find_orfs(Fasta_entry *read,vector<Orf> *open_frames);
    void define_translation_tables();
    string reverse_complement(string dna);

    void find_nodes_for_reads(Node *root, vector<Fasta_entry> *reads, Model_factory *mf);
    void find_nodes_for_all_reads(Node *root, vector<Fasta_entry> *reads, Model_factory *mf);
    void find_nodes_for_all_reads_together(Node *root, vector<Fasta_entry> *reads, Model_factory *mf);
    void remove_overlapping_reads(vector<Fasta_entry> *reads, Model_factory *mf);
    void remove_target_overlapping_identical_reads(vector<Fasta_entry> *reads_for_this, Model_factory *mf);
    void remove_target_overlapping_reads(vector<Fasta_entry> *reads_for_this);
    void align_two_reads(Node *node, Fasta_entry *ri1, Fasta_entry *ri2, Model_factory *mf);
    int reads_pairwise_matching_sites(Node *node);

    double read_match_score(Node *node, Fasta_entry *read, Model_factory *mf, float best_score);
    void read_alignment_scores(Node * node, string read_name, string ref_node_name, float *overlap, float *identity);
    bool read_alignment_overlaps(Node * node, string read_name, string ref_node_name);
    float read_alignment_overlap(Node * node, string read_name, string ref_node_name);
    void add_trimming_comment(vector<Fasta_entry> *reads);
    void merge_paired_reads(vector<Fasta_entry> *reads, Model_factory *mf);
    void find_paired_reads(vector<Fasta_entry> *reads);
    void copy_node_details(Node *reads_node,Fasta_entry *read,bool turn_revcomp = false);
    //void copy_orf_details(Node *reads_node,Fasta_entry *read,Orf *orf,bool turn_revcomp);
    bool correct_sites_index(Node *current_root, string ref_node_name, int alignments_done, map<string,Node*> *nodes_map);

    static bool better_score(const Fasta_entry& a,const Fasta_entry& b)
    {
        return  (a.node_score>b.node_score);
    }

    void sort_reads_vector(vector<Fasta_entry> *reads)
    {
        stable_sort(reads->begin(),reads->end(),Reads_aligner::better_score);
    }

    static bool longer_sequence(const Fasta_entry& a,const Fasta_entry& b)
    {
        return  (a.sequence.length()>b.sequence.length());
    }

    void sort_reads_vector_by_length(vector<Fasta_entry> *reads)
    {
        stable_sort(reads->begin(),reads->end(),Reads_aligner::longer_sequence);
    }

    static bool nodeIsSmaller(const string& l,const string& r)
    {   char a,b,c,d;
        int vl=0; int vr=0;
        stringstream ls(l);
        stringstream rs(r);
        ls>>a>>vl>>b;
        rs>>c>>vr>>d;

        if(a=='#' && b=='#' && c=='#' && d=='#' && vl>0 && vr>0)
            return (vl<vr);

        if(a=='#' && b=='#'&& vl>0)
            return false;

        if(c=='#' && d=='#'&& vr>0)
            return true;

        return (l<r);

    }
public:
    Reads_aligner();
    void align(Node *root, Model_factory *mf,int count);

    void merge_reads_only();

    Node *get_global_root() { return global_root; }
};
}

#endif // READS_ALIGNER_H
