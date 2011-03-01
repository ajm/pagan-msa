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

#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include <string>
#include <vector>
#include "utils/db_matrix.h"
#include "utils/int_matrix.h"
#include "utils/evol_model.h"
#include "utils/settings.h"

namespace ppa{

struct Char_symbol
{
    int index;
    char symbol;
    int n_residues;
    int first_residue;
    int second_residue;
    std::vector<char> residues;
};

class Model_factory
{

    std::string char_alphabet;
    std::string full_char_alphabet;
    int char_as;  // alphabet size
    int char_fas; // full_alphabet size
    int sequence_data_type;

    Db_matrix *charPi;
    float char_ins_rate;
    float char_del_rate;
    float char_ext_prob;
    float char_end_ext_prob;
    float char_break_ext_prob;

    Int_matrix *parsimony_table;
    Int_matrix *child_parsimony_table;

    Db_matrix * charU;
    Db_matrix * charV;
    Db_matrix * charRoot;

//    int i;
//    int j;
//    int k;
//    int l;

    std::vector<Char_symbol> char_letters;

    void define_dna_alphabet();
    void define_protein_alphabet();
    void define_protein_alphabet_old();
    void define_codon_alphabet();

    void print_char_alphabet();

    void build_model(int s,Db_matrix *pi,Db_matrix *q,Db_matrix *wU,Db_matrix *wV,Db_matrix *wRoot);
    void print_char_q_matrices(Db_matrix *charQ);
//    void print_char_p_matrices(Evol_model &model);

public:
    Model_factory(int sequence_data_type);
    ~Model_factory();

    enum Data_type {dna,protein,codon};

    void dna_model(float *char_pi,Settings *st);
    void dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob, float end_ext_prob, float break_ext_prob);

    void protein_model(Settings *st);
    void protein_model(float ins_rate,float del_rate, float ext_prob, float end_ext_prob);

    void codon_model(Settings *st);
    void codon_model(float ins_rate,float del_rate, float ext_prob);

    Evol_model alignment_model(double distance, bool is_local_alignment=false);
//    Evol_model char_alignment_model(double distance);

    void print_int_matrix(Int_matrix *m);
    void print_char_p_matrices(Evol_model &model);

    std::string get_char_alphabet() { return char_alphabet; }
    std::string get_full_char_alphabet() { return full_char_alphabet; }

    int get_child_parsimony_state(int parent_state,int child_state) { return child_parsimony_table->g(parent_state,child_state);}
};

}
#endif // MODEL_FACTORY_H
