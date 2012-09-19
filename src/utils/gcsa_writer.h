#ifndef GCSA_WRITER_H
#define GCSA_WRITER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <cstdio>
#include <tr1/cstdint>

#include "main/node.h"
#include "main/sequence.h"
#include "utils/log_output.h"
#include "utils/model_factory.h"

//#define GCSA_DEBUG 1

using namespace std;
using namespace ppa;

class GCSA_writer {

    char write_buf[BUFSIZ];
    uint32_t* buf_head;
    ofstream f;
    int index;
    
    
    bool _buf_full() {
        return ((char*) buf_head) == (write_buf + BUFSIZ);
    }
    
    void _reset_buf() {
        buf_head = (uint32_t*) write_buf;
    }
    
    void _flush_buf() {
        //f.write(write_buf, BUFSIZ);
        if(f.write(write_buf, (char*)buf_head - write_buf).fail()) {
            Log_output::write_out("Err: write failed flushing buffer to GCSA automata file", 0);
            abort();
        }
        _reset_buf();
    }
    
    void _write_pair(uint32_t a, uint32_t b) {
        *buf_head++ = a;
        *buf_head++ = b;
        
        if(_buf_full()) {
            _flush_buf();
        }
    }
    
    /*
     * WARNING: this only works for proteins, IUPAC codes would need to be translated
     *          into two separate sites and then the edge indices would not be properly
     *          synced
     */
    /*
     start = 21
     A = 1
     R = 2
     N = 3
     D = 4
     C = 5
     Q = 6
     E = 7
     G = 8
     H = 9
     I = 10
     L = 11
     K = 12
     M = 13
     F = 14
     P = 15
     S = 16
     T = 17
     W = 18
     Y = 19
     V = 20
     end = 0
     */
    void _write_sequence(Sequence* seq) {
        vector<Edge>* edges = seq->get_edges();
        vector<Site>* sites = seq->get_sites();
        
        //string alpha = seq->get_full_alphabet();
        string alpha = seq->get_data_type() == Model_factory::protein ? \
                        Model_factory::get_protein_char_alphabet() :
                        Model_factory::get_dna_char_alphabet();
        
        // kill if not protein
        if(seq->get_data_type() != Model_factory::protein) {
            Log_output::write_out("Err: GCSA files can only be written for protein alignments at the moment\n", 0);
            abort();
        }
        
        // metadata
        _write_pair(sites->size(), edges->size()-1);
        
#ifdef GCSA_DEBUG
        cerr << sites->size() << " nodes\n";
        cerr << edges->size()-1 << " edges\n";
#endif
        
        // nodes, start and ends nodes needs special labels
        _write_pair(alpha.size()+1, 0);
#ifdef GCSA_DEBUG
        cerr << "node: " << alpha.size()+1 << ":" << 0 << "\n";
#endif
        
        for(int i = 1; i < int(sites->size() - 1); ++i) {
            _write_pair((*sites)[i].get_state()+1, i);
        
#ifdef GCSA_DEBUG
            cerr << "node: " << (*sites)[i].get_state()+1 << ":" << i << "\n";
#endif
        }
        
        _write_pair(0, sites->size() - 1);
#ifdef GCSA_DEBUG
        cerr << "node: " << 0 << ":" << sites->size()-1 << "\n";
#endif
        
        
        // edges
        for(int i = 1; i < int(edges->size()); ++i) {
            _write_pair((*edges)[i].get_start_site_index(), (*edges)[i].get_end_site_index());
            
#ifdef GCSA_DEBUG
            cerr << "edge: " << (*edges)[i].get_start_site_index() << ":" << (*edges)[i].get_end_site_index() << "\n";
#endif
        }
        
        _flush_buf();
    }
    
public:
    GCSA_writer() :
        buf_head(0),
        f(),
        index(0) {
    
        _reset_buf();
    }
    
    void write(string filebasename, Node* n) {
        
        if(not n->is_leaf()) {
            write(filebasename, n->get_left_child());
            write(filebasename, n->get_right_child());
        }
        
        stringstream ss;
        ss << filebasename << "." << index++;
        
#ifdef GCSA_DEBUG
        //cerr << ss.str() << "\n";
#endif
        
        f.open(ss.str().c_str());
        
        _write_sequence(n->get_sequence());
        
        f.close();
        
#ifdef GCSA_DEBUG
        exit(0);
#endif
    }
};

#endif
