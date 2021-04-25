#ifndef _SCQR
#define _SCQR
#include <stdint.h>
#include <string>
#include <vector>

#define USE_PHMAP
#ifndef USE_PHMAP
#include <set>
#include <map>
#define scqr_set std::set
#define scqr_map std::map
#else
#include <parallel_hashmap/phmap.h>
#include <parallel_hashmap/btree.h>
#define scqr_set phmap::btree_set
#define scqr_map phmap::btree_map

#endif

#define _RUNMODE_ 1
#if ((_RUNMODE_))
const uint32_t barcode = 0x10;
const uint32_t umi     = 0xC;
#else
const uint32_t barcode = 0x8;
const uint32_t umi     = 0x6;
#endif

const uint32_t polyT      = 0x8;
const uint32_t polyTMask  = 0xFFFF;
const uint32_t polyTMask1 = 0xFFF0;
#define jumpbit 0

typedef uint64_t tag_t;
typedef uint32_t umi_t;
typedef uint64_t barcode_t;

typedef struct
{
        umi_t    umi;
        uint32_t sample; //barcode
} rOneRead;

typedef struct
{
        barcode_t sample;
        uint32_t  count;
        barcode_t _sample;
        uint32_t  _count;
} barcodeCount_t;

typedef struct
{
        barcode_t   sample;
        std::string seq;
        //std::string q; // quality
} rSecondRead;

typedef struct
{
        barcode_t bc;
        umi_t     umitag;
} read_t;

int test(){
        scqr_map<barcode_t,uint32_t> a;
        scqr_set<uint32_t > b;
        return 0;
}


class ss{
 public:
        scqr_set<uint32_t> s;
};

#endif