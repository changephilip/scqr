#include "../include/readR1.h"


typedef struct
{
        //std::string pair2;
        //uint32_t readR1Id;
        uint16_t sample;
        std::string seq;
        //std::string q; // quality
} rSecondRead;

class rSecond {
 public:
        std::vector<rSecondRead> reads;
        std::vector<std::vector<rSecondRead>> readsPack;
        rSecond(const std::string &r2gz, const rFirst & rfirst);
};

rSecond::rSecond(const std::string &r2gz, const rFirst & rfirst){
        gzFile fp;
        kseq_t *seq;
        int l;
        fp = gzopen(r2gz.c_str(),"r");
        seq = kseq_init(fp);
        uint32_t readId = 0;

        std::time_t t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tStart of R2 " << std::endl;

        readsPack.resize(rfirst.sampleId.size());
        while (( l = kseq_read(seq)) >=0){
                std::string thisName = seq->name.s;
                if (rfirst.readsNameTable.find(thisName)!= rfirst.readsNameTable.end()){
                        uint32_t seqId = rfirst.readsNameTable.find(thisName)->second;
                        rSecondRead tmp;
                        tmp.seq = seq->seq.s;
                        //tmp.readR1Id = seqId;
                        tmp.sample = rfirst.reads[seqId].sample;
                        //tmp.sample = rfirst.sampleList.find(seqId)->second;
                        //tmp.q = seq->qual.s;
                        //assert(rfirst.readsNameTable.find(thisName).second);
                        //reads.push_back(tmp);
                        readsPack[tmp.sample].push_back(tmp);
                }
                readId ++;
        }

        std::cout << std::asctime(std::localtime(&t)) << "\tEnd of R2 " << std::endl;

}
