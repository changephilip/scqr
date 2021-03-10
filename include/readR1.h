#include "../include/mkdb.h"
#include <stdint.h>

const uint32_t barcode = 0x10;
const uint32_t umi     = 0xC;
typedef struct
{
        uint32_t tag;
        uint32_t id;
        uint16_t sample;
        std::string pair1;
} rOneRead;

class rFirst {
       public:
        std::set<uint32_t>       seqId;
        std::vector<rOneRead>    reads;
        std::set<uint16_t> sample;
        std::vector<uint16_t> sampleV;
        std::map<uint16_t, uint16_t> sampleList;
        std::map<std::string , uint32_t> readsNameTable;
        rFirst(const std::string r1gz);
};

inline uint16_t readBarcode(const char *seq)
{
        return (uint16_t)baseToBinaryForward(seq, barcode);
}
inline uint16_t readUmi(const char *seq)
{
        return (uint16_t)baseToBinaryForward(seq+barcode, umi);
}

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
                uint32_t thistag     = ((uint32_t)thisbarcode << 12) || thisumi;
                if (seqId.insert(thistag).first != seqId.end())
                {
                        rOneRead tmpRead;
                        tmpRead.tag = thistag;
                        tmpRead.id  = readId;
                        //tmpRead.sample = thisbarcode;
                        tmpRead.pair1 = seq->name.s;
                        reads.push_back(tmpRead);
                        sample.insert(thisbarcode);
                        readsNameTable.insert(std::pair<std::string, uint32_t>(seq->name.s,readId));
                }
                readId++;
        }
        std::cout << reads.size() << std::endl;
        std::cout << seqId.size() << std::endl;
        assert(reads.size() == seqId.size());

        auto it = sample.begin();
        uint16_t count=0;
        while (it!=sample.end()){
                sampleList.insert(std::pair<uint16_t,uint16_t>(*it,count));
                count++;
                sampleV.push_back(*it);
        }

        for (auto &read: reads){
                read.sample = sampleList.find(read.tag>>12)->second;
        }

}


