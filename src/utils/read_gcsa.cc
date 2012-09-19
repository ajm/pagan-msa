#include <iostream>
#include <fstream>

#include <tr1/cstdint>

using namespace std;

void read_int(ifstream& f, uint32_t* a) {
    if(f.read((char*)a, sizeof(*a)).fail()) {
        cerr << "read fail\n";
        exit(-1);
    }
}

int main(int argc, char** argv) {
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << " <gcsa file>\n";
        exit(-1);
    }

    ifstream f;
    f.open(argv[1]);

    if(not f.good()) {
        cerr << "could not read '" << argv[1] << "'\n";
        exit(-1);
    }

    uint32_t a, b, numnodes, numedges;
    
    read_int(f, &numnodes);
    read_int(f, &numedges);
    
    cout << numnodes << " nodes\n";
    cout << numedges << " edges\n";

    for(int i = 0; i < int(numnodes); ++i) {
        read_int(f, &a);
        read_int(f, &b);
        cout << "node: " << a << ":" << b << "\n";
    }
    
    for(int i = 0; i < int(numedges); ++i) {
        read_int(f, &a);
        read_int(f, &b);
        cout << "edge: " << a << ":" << b << "\n";
    }
    
    f.close();

    return EXIT_SUCCESS;
}

