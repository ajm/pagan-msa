#include <sys/mman.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <arpa/inet.h>
#include <ctype.h>

#include "sff_reader.h"

extern int errno;


uint32_t pad(uint32_t i) {
    return (i % 8) == 0 ? i : i + (8 - (i % 8));
}

void get_read(char* baseaddr, uint16_t flow_len, uint32_t nbases, uint16_t trim_left, uint16_t trim_right) {
    
    uint16_t* flowgram  = (uint16_t*) baseaddr;
    uint8_t* flow_index = (uint8_t*)  baseaddr + (flow_len * 2);
    char* bases         = (char*)     baseaddr + (flow_len * 2) + nbases;
    uint8_t* quality    = (uint8_t*)  baseaddr + (flow_len * 2) + (nbases * 2);
    
    
    for(int i = 0; i < nbases; ++i) {
        printf("%c", (i < trim_left) || (i > trim_right) ? tolower(bases[i]) : toupper(bases[i]));
    }
    printf("\n");
    
    /*
    for(int i = 0; i < flow_len; ++i) {
        printf("%.2f ", ntohs(flowgram[i]) / 100.0); 
    }
    printf("\n");
    */
    
    uint16_t total = 0;
    for(int i = 0; i < nbases; ++i) {
        total += flow_index[i];
        printf("%.2f ", ntohs(flowgram[total-1]) / 100.0);
    }
    printf("\n");
    
    for(int i = 0; i < nbases; ++i) {
        printf("%d ", quality[i]);
    }
    printf("\n");
}

void read_sff_file(char* addr) {
    sff_common_header* head = (sff_common_header*) addr;
    
    if(ntohl(head->magic) != SFF_MAGIC) {
        fprintf(stderr, "Warning: Bad magic number\n");
    }

    uint32_t nreads     = ntohl(head->nreads);
    uint16_t header_len = ntohs(head->header_len);
    uint16_t flow_len   = ntohs(head->flow_len);
    char* current_addr  = addr + header_len;
        
    for(int i = 0; i < int(nreads); ++i) {
        sff_read_header* read = (sff_read_header*) current_addr;
        uint16_t rheader_len = ntohs(read->header_len);
        uint32_t nbases = ntohl(read->nbases);
        uint16_t name_len = ntohs(read->name_len);
        uint16_t trim_left = ntohs(read->clip_qual_left) - 1; // these appear to be indexed from 1
        uint16_t trim_right = ntohs(read->clip_qual_right) - 1;
        
        /*
        printf( "header_len = %d\n"
                "name_len = %d\n"
                "nbases = %d\n\n",
                ntohs(read->header_len),
                ntohs(read->name_len),
                ntohl(read->nbases));
        */
        
        get_read(current_addr + rheader_len, flow_len, nbases, trim_left, trim_right);
        
        //exit(EXIT_SUCCESS);
                
        current_addr = current_addr + rheader_len + pad((flow_len * 2) + (nbases * 3));
    }
}

int main(int argc, char** argv) {
    int fd;
    struct stat sb;
    char* addr;
    
    fd = open(argv[1], O_RDONLY);
    if(fd == -1) {
        fprintf(stderr, "Error opening '%s': %s\n", argv[1], strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    if(fstat(fd, &sb) == -1) {
        fprintf(stderr, "Error fstat-ing '%s': %s\n", argv[1], strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    addr = (char*) mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if(addr == MAP_FAILED) {
        fprintf(stderr, "Error mmap-ing '%s': %s\n", argv[1], strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    read_sff_file(addr);
    
    close(fd);
    
    return EXIT_SUCCESS;
}
