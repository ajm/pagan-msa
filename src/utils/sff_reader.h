#ifndef SFF_FILE_H
#define SFF_FILE_H

#include <stdint.h>


#define SFF_MAGIC   0x2e736666 /* ".sff" */

// structs stolen from sff2fasta sff.h
// https://github.com/indraniel/sff2fastq/blob/master/sff.h

typedef struct {
    uint32_t  magic;
    char      version[4];
    uint64_t  index_offset;
    uint32_t  index_len;
    uint32_t  nreads;
    uint16_t  header_len;
    uint16_t  key_len;
    uint16_t  flow_len;
    uint8_t   flowgram_format;
    char     *flow;
    char     *key;
    
} sff_common_header;

/*
 * There is a read_header per "reading" in the SFF archive.
 * It too is padded to an 8-byte boundary.
 */
typedef struct {
    uint16_t  header_len;
    uint16_t  name_len;
    uint32_t  nbases;
    uint16_t  clip_qual_left;
    uint16_t  clip_qual_right;
    uint16_t  clip_adapter_left;
    uint16_t  clip_adapter_right;
    char     *name;
    
} sff_read_header;

/*
 * The read_data section per reading, following the read_header.
 * It is padded to an 8-byte boundary.
 */
typedef struct {
    uint16_t *flowgram;   /* x 100.0 */
    uint8_t  *flow_index; /* relative to last */
    char     *bases;
    uint8_t  *quality;
    
} sff_read_data;



#endif
