#ifndef _SCQR
#define _SCQR
#include <stdint.h>
#include <string>

#define USE_PHMAP
#ifndef USE_PHMAP
#include <set>
#include <map>
#define scqr_set std::set
#define scqr_map std::map
#else
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/btree.h"
#define scqr_set phmap::btree_set
#define scqr_map phmap::btree_map

#endif


#ifdef __USE_ZSTD
#define ZWRAP_USE_ZSTD 1
#include "zstd_zlibwrapper.h"
#else
#include <zlib.h>
#endif

#define BARCODE_NULL 0xFFFFFFFF


//#define __LHJ

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
        uint32_t sample; //barcode id,a number for each barcode
        uint32_t id;
}rOneRead;

typedef struct
{
        barcode_t sample;
        uint32_t  count;
        barcode_t _sample;
        uint32_t  _count;
} barcodeCount_t;

typedef struct
{
        //uint32_t   sample;
        std::string seq;
        //std::string q; // quality
} rSecondRead;

typedef struct
{
        barcode_t bc;
        umi_t     umitag;
} read_t;

#endif
