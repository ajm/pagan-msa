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
        f.write(write_buf, BUFSIZ);
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
    // TODO: double check with Ari that the first and last nodes are not part of the
    //       actual sequence
    void _write_sequence(Sequence* seq) {
        string alpha = seq->get_full_alphabet();
        vector<Edge>* edges = seq->get_edges();
        vector<Site>* sites = seq->get_sites();
        
        // kill if not protein
        if(seq->get_data_type() != Model_factory::protein) {
            Log_output::write_out("Err: GCSA files can only be written for protein alignments at the moment\n", 0);
            abort();
        }
        
        // metadata
        _write_pair(sites->size(), edges->size());
        
        
        // nodes, start and ends nodes needs special labels
        _write_pair(alpha.size()+1, 0);
        
        for(int i = 1; i < int(sites->size() - 1); ++i) {
            _write_pair((*sites)[i].get_state()+1, i);
        }
        
        _write_pair(0, sites->size() - 1);
        
        
        // edges
        for(int i = 0; i < int(edges->size()); ++i) {
            _write_pair((*edges)[i].get_start_site_index(), (*edges)[i].get_end_site_index());
        }
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
        
        f.open(ss.str().c_str());
        
        _write_sequence(n->get_sequence());
        
        f.close();
    }
};

#endif
