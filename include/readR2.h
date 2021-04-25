#ifndef __READR2
#define __READR2
#include "scqr.h"
#include "readR1.h"



class rSecond {
 public:
        std::vector<rSecondRead> reads;
        std::vector<std::vector<rSecondRead>> readsPack;
        rSecond(const std::string &r2gz, const rFirst & rfirst);
        rSecond();
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

        readsPack.resize(rfirst.correctedBarcodeVector.size());
        while (( l = kseq_read(seq)) >=0){
                //auto nameHash = rfirst.stringHash(seq->name.s);
                auto nameHash = rfirst.readsNameTable.find(seq->name.s);
                if (nameHash != rfirst.readsNameTable.end()){
                        uint32_t seqId = nameHash->second;
                        rSecondRead tmp;
                        tmp.seq = seq->seq.s;
                        //tmp.readR1Id = seqId;
                        if (rfirst.reads[seqId].sample != 0xFFFFFFFF)
                        {
                                tmp.sample = rfirst.reads[seqId].sample;
                                //tmp.sample = rfirst.sampleList.find(seqId)->second;
                                //tmp.q = seq->qual.s;
                                //assert(rfirst.readsNameTable.find(thisName).second);
                                //reads.push_back(tmp);
                                readsPack[tmp.sample].push_back(tmp);
                        }
                }
                readId ++;
        }

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tEnd of R2 " << std::endl;

}

#endif
