#include "../include/mkdb.h"
#include <stdint.h>

#include <immintrin.h>
#include <popcntintrin.h>

#define _RUNMODE_ 0
#if ((_RUNMODE_))
const uint32_t barcode = 0x10;
const uint32_t umi     = 0xC;
#else
const uint32_t barcode = 0x8;
const uint32_t umi     = 0x6;
#endif

const uint32_t polyT = 0x8;
const uint32_t polyTMask = 0xFFFF;
const uint32_t polyTMask1 = 0xFFF0;
#define jumpbit 0

typedef uint64_t tag_t;
typedef uint32_t umi_t;
typedef uint64_t barcode_t;

typedef struct
{
        tag_t tag;
        //uint32_t id;
        umi_t umi;
        barcode_t sample; //barcode
        //std::string pair1;
} rOneRead;

class rFirst {
       public:
        /*a unique identifier, tag=[barcode]+[umi]*/
        std::vector<rOneRead> reads;
        std::set<barcode_t>    barcodeSet;
        /*assign 0-based id for each sample*/
        std::map<barcode_t, uint32_t> barcodeMapList;
        /*get origin sampleId(barcode) with index*/
        std::vector<barcode_t> barcodeVector;
        /*store pair1 string for R2 search for the paired reads R1 readId*/
        std::map<uint64_t, uint32_t> readsNameTable;

        std::vector<barcode_t> barcodeDataBase;
        rFirst(const std::string r1gz);
        std::hash<std::string> stringHash;
        void barcodeCorrect();
};
/*
template<class T>
inline bool med(T in,T target){
        uint32_t bitWise = sizeof(T)*8;
        uint32_t p = in xor target ;
        
}
*/

inline barcode_t readBarcode(const char *seq)
{
        //return (barcode_t)baseToBinaryForward_Barcode(seq + jumpbit, barcode);
        return (barcode_t)baseToBinaryForward(seq + jumpbit, barcode);
}

inline umi_t readUmi(const char *seq)
{
        return (umi_t)baseToBinaryForward(seq + barcode + jumpbit , umi);
}

inline tag_t readTag(const char *seq)
{
        return (tag_t)baseToBinaryForward(seq + jumpbit , umi + barcode);
}

inline uint64_t  checkPolyT(const char *seq){
        uint64_t T = baseToBinaryForward(seq + jumpbit + umi +barcode, polyT);
        T |= 0ULL;
        uint32_t p = (T xor polyTMask);
        uint32_t shift = _mm_popcnt_u64(p);
        uint32_t q = 32;
        uint32_t count =0;
        while (q != 0 ){
        }
        while (q < polyT)
        {
                if (*(seq + jumpbit + umi + barcode) != 'T'){
                        
                }
        }
        return 0;
}

/*
void rFirst::barcodeCorrect(){
        auto dbscan = DBSCAN<uint64_t, uint64_t>();

        auto data = std::vector<uint64_t>{barcodeSet.begin(),barcodeSet.end()};

        dbscan.Run(&data, 1, 1, 4000);
        auto noise = dbscan.Noise;
        auto clusters = dbscan.Clusters;
}
*/

/*
  filter unique tag(barcode + umi) and get these reads from R1
  save out these sample(barcode) id and transfer into 0-based index.
  How to find pair-end reads?
  construct readsNameTable with pair of (name.s , readId)
*/
rFirst::rFirst(const std::string r1gz)
{
        gzFile  fp;
        kseq_t *seq;
        int     l;
        fp              = gzopen(r1gz.c_str(), "r");
        seq             = kseq_init(fp);
        uint32_t readId = 0;

        std::vector<std::pair<uint64_t,uint32_t>> readsNameVector;
        std::time_t t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tStart of R1 " << std::endl;

        while ((l = kseq_read(seq)) >= 0)
        {
                barcode_t thisbarcode = readBarcode(seq->seq.s);
                umi_t thisumi     = readUmi(seq->seq.s);
                //uint32_t thistag     = ((uint32_t)thisbarcode << 12) || (0x0u || thisumi);
                {
                        rOneRead tmpRead;
                        tmpRead.umi = thisumi;
                        //tmpRead.id  = readId;
                        //tmpRead.sample = thisbarcode;
                        //tmpRead.pair1 = seq->name.s;
                        reads.push_back(tmpRead);
                        barcodeSet.insert(thisbarcode);
                        /*
                        readsNameTable.insert(
                                std::pair<std::string, uint32_t>(seq->name.s, readId));
                        */
                        //readsNameVector.push_back(std::make_pair(seq->name.s, readId));
                        //readsNameVector.push_back(std::make_pair(std::hash<std::string>{}(seq->name.s),readId));
                        readsNameVector.push_back(std::make_pair(stringHash(seq->name.s),readId));
                        readId++;
                }
                //readId++;
        }
        
        readsNameTable = std::map<uint64_t,uint32_t>(readsNameVector.begin(),readsNameVector.end());

        readsNameVector.resize(0);
        readsNameVector.clear();



        std::cout << reads.size() << std::endl;
        std::cout << barcodeSet.size() << std::endl;
        // reads.size() >> seqId.size()
        //assert(reads.size() == seqId.size());

        auto     it = barcodeSet.begin();
        barcode_t count = 0;
        while (it != barcodeSet.end())
        {
                barcodeMapList.insert(std::pair<barcode_t, uint32_t>(*it, count));
                barcodeVector.push_back(*it);
                it++;
                count++;
        }

        //for (auto &read : reads)
        //{
        //        read.sample = sampleList.find(read.tag >> 12)->second;
        //}

        #pragma omp parallel for
        for (uint32_t i=0;i< reads.size();i++){
                reads[i].sample = barcodeMapList.find(reads[i].tag >> 12)->second;
        }

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tEnd of R1 " << std::endl;


}
