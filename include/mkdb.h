#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <stdio.h>
#include <cassert>

#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include <zlib.h>
#include <omp.h>
#include <parallel/algorithm>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

const uint32_t rNATailLength        = 128;
const uint32_t kmerLength           = 32;
const uint64_t conversionTable[128] = {4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       //32
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       //64 A=65,C=67,G=71,T=84
                                       4,
                                       0,
                                       4,
                                       1,
                                       4,
                                       4,
                                       4,
                                       2,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       3,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       //96 a=97,c=99,g=103,t=116
                                       4,
                                       0,
                                       4,
                                       1,
                                       4,
                                       4,
                                       4,
                                       2,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       3,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4,
                                       4};

enum strainType
{
        forward,
        backward,
        forward_supply,
        backward_supply
};
enum baseType
{
        A,
        C,
        G,
        T
};

typedef struct
{
        uint32_t    transId;
        uint32_t    length = 0;
        char*       seq;
        std::string description;
} transFa;

inline uint64_t baseToBinaryForward(const char* seq, uint64_t length)
{
        uint64_t p        = 0;
        uint64_t bitwidth = (length - 1) * 2;
        uint64_t tmp      = 0ULL;

        while (p < length)
        {
                tmp |= (conversionTable[seq[p]] << bitwidth);
                p++;
                bitwidth -= 2;
        }

        return tmp;
}

inline uint64_t baseToBinaryBackward(const char* seq, uint64_t length)
{
        int64_t  p        = length - 1;
        uint64_t bitwidth = (length - 1) * 2;
        uint64_t tmp      = 0ULL;
        while (p >= 0)
        {
                tmp |= (conversionTable[seq[p]] << bitwidth);
                p--;
                bitwidth -= 2;
        }
        return tmp;
}

inline uint64_t baseToBinary(const char* seq, uint64_t length, strainType strain)
{
        if (strain == forward)
                return baseToBinaryForward(seq, length);
        else
        {
                return baseToBinaryBackward(seq, length);
        }
}
inline bool forwardOrBackward(const char* seq, uint64_t length)
{
        auto forward  = baseToBinaryForward(seq, length);
        auto backward = baseToBinaryBackward(seq, length);
        return forward > backward;
}
typedef struct
{
        uint32_t transId;
        uint64_t kmer;
} kmer_t;

inline void sort_kmer(std::vector<kmer_t> &table, const int32_t thread)
{
        omp_set_dynamic(false);
        omp_set_num_threads(thread);
        std::sort(table.begin(),table.end(), [](kmer_t left,kmer_t right){
                        return left.kmer < right.kmer;
                });
}

class transcriptomeFa {
       public:
        std::vector<transFa>  originFa;
        uint32_t              transNum = 0;
        std::vector<uint32_t> dbId;
        std::vector<uint64_t> dbEntryValue;
        std::vector<uint32_t> hashedDBId;
        std::vector<kmer_t>   kmerTable;
        std::vector<kmer_t> contact;
        uint32_t thread=1;
        
        transcriptomeFa(const char* fa, uint32_t thread_);
        void addKmer(transFa& rna);
        void mkdb();
};

transcriptomeFa::transcriptomeFa(const char* fa, uint32_t thread_)
{
        gzFile  fp;
        kseq_t* seq;
        int     l;
        fp  = gzopen(fa, "r");
        seq = kseq_init(fp);
        while ((l = kseq_read(seq)) >= 0)
        {
                transFa tmp;
                tmp.transId     = transNum;
                tmp.description = seq->name.s;
                tmp.length      = strlen(seq->seq.s);
                tmp.seq         = new char[tmp.length + 1];
                memcpy(tmp.seq, seq->seq.s, sizeof(char) * (tmp.length + 1));
                tmp.seq[tmp.length] = '\0';
                this->originFa.push_back(tmp);
                transNum++;
        }
        this->thread = thread_;
        mkdb();
}

void transcriptomeFa::addKmer(transFa& rna)
{
        std::vector<kmer_t> kmer;
        uint32_t            cutLength = rNATailLength;
        if (rna.length < cutLength)
        {
                cutLength = rna.length;
        }

        char* p         = NULL;
        p               = rna.seq + rna.length - cutLength;
        uint32_t n_kmer = cutLength - kmerLength + 1;
        for (uint64_t i = 0; i < n_kmer; i++)
        {
                kmer_t     tmpK;
                strainType thisStrainType;
                if (forwardOrBackward(p, kmerLength))
                {
                        thisStrainType = forward;
                }
                else
                {
                        thisStrainType = backward;
                }

                tmpK.kmer    = baseToBinary(p, kmerLength, thisStrainType);
                tmpK.transId = rna.transId;
                kmer.emplace_back(tmpK);
                p++;
        }

        //supply
        p = rna.seq + rna.length - cutLength;
        char* supplyRNA;
        supplyRNA = new char[cutLength + 1];
        for (uint64_t i = 0; i < cutLength; i++)
        {
                supplyRNA[i] = *p;
                p++;
        }
        supplyRNA[cutLength] = '\0';
        p                    = supplyRNA;
        for (uint64_t i = 0; i < n_kmer; i++)
        {
                kmer_t     tmpK;
                strainType thisStrainType;
                if (forwardOrBackward(p, kmerLength))
                {
                        thisStrainType = forward;
                }
                else
                {
                        thisStrainType = backward;
                }

                tmpK.kmer    = baseToBinary(p, kmerLength, thisStrainType);
                tmpK.transId = rna.transId;
                kmer.emplace_back(tmpK);
                p++;
        }
        
        for (auto item : kmer){
                kmerTable.emplace_back(item); 
        }
}

void transcriptomeFa::mkdb()
{
        //check if each transcriptome is longer than 256
        uint32_t sum      = 0;
        uint32_t count150 = 0;
        uint32_t count256 = 0;
        /*
        for (auto iter : this->originFa)
        {
                sum += iter.length;
                if (iter.length < 150)
                {
                        count150++;
                }
                if (iter.length < 256)
                {
                        count256++;
                }

        }
        std::cout << count150 << "\t" << count256 << std::endl;
        std::cout << sum / float(this->originFa.size()) << std::endl;
        */
        for (auto rna: this->originFa){
                addKmer(rna);
        }

        sort_kmer(kmerTable, thread);
        auto iter = kmerTable.begin()++;
        auto prev = kmerTable.begin();
        while (iter != kmerTable.end()){
                if (iter->kmer != prev->kmer) {
                        if (iter - prev ==1){
                                contact.push_back(*prev);
                                prev=iter;
                        }else{
                                prev=iter;
                        }
                }
                iter++;
        }
        iter = contact.begin()++;
        while (iter!=contact.end()){
                assert(iter->kmer != (iter-1)->kmer);
        }
        std::cout << contact.size() << std::endl;
}
