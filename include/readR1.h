#include "../include/mkdb.h"
#include <stdint.h>

const uint32_t barcode = 0x10;
const uint32_t umi     = 0xC;
typedef struct
{
        uint32_t tag;
        uint32_t id;
        uint16_t sample; //barcode
        //std::string pair1;
} rOneRead;

class rFirst {
       public:
        /*a unique identifier, tag=[barcode]+[umi]*/
        std::set<uint32_t>    seqId;
        std::vector<rOneRead> reads;
        std::set<uint16_t>    sampleSet;
        /*assign 0-based id for each sample*/
        std::map<uint16_t, uint16_t> sampleList;
        /*get origin sampleId(barcode) with index*/
        std::vector<uint16_t> sampleId;
        /*store pair1 string for R2 search for the paired reads R1 readId*/
        std::map<std::string, uint32_t> readsNameTable;

        rFirst(const std::string r1gz);
};

inline uint16_t readBarcode(const char *seq)
{
        return (uint16_t)baseToBinaryForward(seq, barcode);
}
inline uint16_t readUmi(const char *seq)
{
        return (uint16_t)baseToBinaryForward(seq + barcode, umi);
}

inline uint32_t readTag(const char *seq)
{
        return (uint32_t)baseToBinaryForward(seq, umi + barcode);
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
        while ((l = kseq_read(seq)) >= 0)
        {
                uint16_t thisbarcode = readBarcode(seq->seq.s);
                uint16_t thisumi     = readUmi(seq->seq.s);
                //uint32_t thistag     = ((uint32_t)thisbarcode << 12) || (0x0u || thisumi);
                uint32_t thistag = readTag(seq->seq.s);
                if (seqId.insert(thistag).first != seqId.end())
                {
                        rOneRead tmpRead;
                        tmpRead.tag = thistag;
                        tmpRead.id  = readId;
                        //tmpRead.sample = thisbarcode;
                        //tmpRead.pair1 = seq->name.s;
                        reads.push_back(tmpRead);
                        sampleSet.insert(thisbarcode);
                        readsNameTable.insert(
                                std::pair<std::string, uint32_t>(seq->name.s, readId));
                        readId++;
                }
                //readId++;
        }
        std::cout << reads.size() << std::endl;
        std::cout << seqId.size() << std::endl;
        assert(reads.size() == seqId.size());

        auto     it    = sampleSet.begin();
        uint16_t count = 0;
        while (it != sampleSet.end())
        {
                sampleList.insert(std::pair<uint16_t, uint16_t>(*it, count));
                sampleId.push_back(*it);
                it++;
                count++;
        }

        for (auto &read : reads)
        {
                read.sample = sampleList.find(read.tag >> 12)->second;
        }
}
