#include "../include/mkdb.h"
#include <stdint.h>

#define _RUNMODE_ 0
#if ((_RUNMODE_))
const uint32_t barcode = 0x10;
const uint32_t umi     = 0xC;
#else
const uint32_t barcode = 0x8;
const uint32_t umi     = 0x6;
#endif

#define jumpbit 1

typedef uint64_t tag_t;
typedef uint32_t umi_t;
typedef uint32_t sample_t;

typedef struct
{
        tag_t tag;
        //uint32_t id;
        sample_t sample; //barcode
        //std::string pair1;
} rOneRead;

class rFirst {
       public:
        /*a unique identifier, tag=[barcode]+[umi]*/
        std::set<tag_t>    seqId;
        std::vector<rOneRead> reads;
        std::set<sample_t>    sampleSet;
        /*assign 0-based id for each sample*/
        std::map<sample_t, uint32_t> sampleList;
        /*get origin sampleId(barcode) with index*/
        std::vector<sample_t> sampleId;
        /*store pair1 string for R2 search for the paired reads R1 readId*/
        std::map<std::string, uint32_t> readsNameTable;

        rFirst(const std::string r1gz);
};

inline sample_t readBarcode(const char *seq)
{
        return (sample_t)baseToBinaryForward(seq + jumpbit, barcode);
}
inline umi_t readUmi(const char *seq)
{
        return (umi_t)baseToBinaryForward(seq + barcode + jumpbit , umi);
}

inline tag_t readTag(const char *seq)
{
        return (tag_t)baseToBinaryForward(seq + jumpbit , umi + barcode);
}
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

        std::vector<std::pair<std::string,uint32_t>> readsNameVector;
        std::time_t t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tStart of R1 " << std::endl;

        while ((l = kseq_read(seq)) >= 0)
        {
                sample_t thisbarcode = readBarcode(seq->seq.s);
                umi_t thisumi     = readUmi(seq->seq.s);
                //uint32_t thistag     = ((uint32_t)thisbarcode << 12) || (0x0u || thisumi);
                tag_t thistag = readTag(seq->seq.s);
                if (seqId.insert(thistag).first != seqId.end())
                {
                        rOneRead tmpRead;
                        tmpRead.tag = thistag;
                        //tmpRead.id  = readId;
                        //tmpRead.sample = thisbarcode;
                        //tmpRead.pair1 = seq->name.s;
                        reads.push_back(tmpRead);
                        sampleSet.insert(thisbarcode);
                        /*
                        readsNameTable.insert(
                                std::pair<std::string, uint32_t>(seq->name.s, readId));
                        */
                        readsNameVector.push_back(std::make_pair(seq->name.s, readId));
                        readId++;
                }
                //readId++;
        }

        readsNameTable = std::map<std::string,uint32_t>(readsNameVector.begin(),readsNameVector.end());

        readsNameVector.resize(0);
        readsNameVector.clear();

        std::cout << reads.size() << std::endl;
        std::cout << seqId.size() << std::endl;
        // reads.size() >> seqId.size()
        //assert(reads.size() == seqId.size());

        auto     it    = sampleSet.begin();
        sample_t count = 0;
        while (it != sampleSet.end())
        {
                sampleList.insert(std::pair<uint16_t, uint16_t>(*it, count));
                sampleId.push_back(*it);
                it++;
                count++;
        }

        //for (auto &read : reads)
        //{
        //        read.sample = sampleList.find(read.tag >> 12)->second;
        //}

        #pragma omp parallel for
        for (uint32_t i=0;i< reads.size();i++){
                reads[i].sample = sampleList.find(reads[i].tag >> 12)->second;
        }

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tEnd of R1 " << std::endl;


}
