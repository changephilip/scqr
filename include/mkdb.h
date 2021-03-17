#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <stdint.h>
#include <stdio.h>
#include <cassert>
#include <ctime>

#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>


//#include <fmt/ranges.h>
#include <zlib.h>
#include <omp.h>
#include <parallel/algorithm>

#include <regex>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

const uint32_t rNATailLength        = 150;
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
const char supplyTable[] = "ACGT";

typedef struct
{
        uint32_t    transId;
        uint32_t    length = 0;
        char*       seq;
        std::string description;
        std::string geneName;
        uint32_t    geneId;
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
        uint64_t kmer;
        uint32_t geneId;
} kmer_t;

inline void sort_kmer(std::vector<kmer_t>& table, const int32_t thread)
{
        omp_set_dynamic(false);
        omp_set_num_threads(thread);
        __gnu_parallel::sort(table.begin(), table.end(), [](kmer_t left, kmer_t right) {
                return left.kmer < right.kmer;
        });
}

typedef struct
{
        uint32_t index_start;
        uint32_t index_end;
} index_t;

class transcriptomeFa {
       public:
        std::vector<transFa>               originFa;
        std::vector<std::string>           geneList;
        uint32_t                           transNum = 0;
        std::vector<uint32_t>              dbId;
        std::vector<uint64_t>              dbEntryValue;
        std::vector<uint32_t>              hashedDBId;
        std::vector<std::vector<uint64_t>> gene_T;
        std::vector<kmer_t>                geneKmerTable;
        std::vector<kmer_t>                contact;
        uint32_t                           index[256];
        uint32_t                           thread = 1;

        transcriptomeFa(const char* fa, uint32_t thread_);
        void addKmer(transFa& rna);
        void mkdb();
        void mkindex();
        void writedb();
        void clear();
};

transcriptomeFa::transcriptomeFa(const char* fa, uint32_t thread_)
{
        gzFile                fp;
        kseq_t*               seq;
        int                   l;
        std::set<std::string> geneSet;
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
                tmp.seq[tmp.length]     = '\0';
                std::string geneid      = seq->comment.s;
                uint32_t    leftbracket = geneid.find_last_of("(");
                //if (geneid.find_first_of("(") != leftbracket){
                //        std::cout << tmp.description << seq->comment.s << std::endl;
                //}
                uint32_t rightbracket = geneid.find_last_of(")");
                tmp.geneName =
                        geneid.substr(leftbracket + 1, rightbracket - leftbracket - 1);
                geneSet.insert(tmp.geneName);
                this->originFa.push_back(tmp);
                transNum++;
        }
        for (auto item : geneSet)
        {
                geneList.push_back(item);
        }

        std::map<std::string, uint32_t> geneTable;

        for (auto item : geneSet)
        {
                //std::cout << item << std::endl;
        }

        for (uint64_t i = 0; i < geneList.size(); i++)
        {
                geneTable.insert({geneList[i], i});
        }
        for (auto& rna : originFa)
        {
                //rna.geneId = geneTable[rna.geneName];
                auto iter = geneTable.find(rna.geneName);
                assert(iter != geneTable.end() and iter->second < geneTable.size());
                rna.geneId = geneTable.find(rna.geneName)->second;
                //std::cout << geneTable.find(rna.geneName)->second << std::endl;
        }
        gene_T.resize(geneSet.size());
        this->thread = thread_;
        mkdb();
        mkindex();
        writedb();
        clear();
}

void transcriptomeFa::addKmer(transFa& rna)
{
        std::vector<kmer_t> pkmer;
        std::vector<kmer_t> rkmer;
        uint32_t            cutLength = 0;
        if (rna.length >= rNATailLength)
        {
                cutLength = rNATailLength;
        }
        else if (rna.length < cutLength and rna.length > kmerLength)
        {
                cutLength = rna.length;
        }
        else
        {
                return;
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

                tmpK.kmer = baseToBinary(p, kmerLength, thisStrainType);
                //tmpK.transId = rna.transId;
                tmpK.geneId = rna.geneId;
                pkmer.push_back(tmpK);
                p++;
        }

        //supply

        p = rna.seq + rna.length - cutLength;
        char* supplyRNA;
        supplyRNA = new char[cutLength + 1];
        for (uint64_t i = 0; i < cutLength; i++)
        {
                supplyRNA[i] = supplyTable[3 - conversionTable[*p]];
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

                tmpK.kmer = baseToBinary(p, kmerLength, thisStrainType);
                //tmpK.transId = rna.transId;
                rkmer.emplace_back(tmpK);
                p++;
        }

        for (auto item : pkmer)
        {
                gene_T[rna.geneId].push_back(item.kmer);
        }
        for (auto item : rkmer)
        {
                gene_T[rna.geneId].push_back(item.kmer);
        }
}

void redupKmer(std::vector<kmer_t>& in, std::vector<kmer_t>& out)
{
        auto iter = in.begin();
        iter++;
        auto prev = in.begin();
        while (iter != in.end())
        {
                if (iter->kmer != prev->kmer)
                {
                        if (std::distance(prev, iter) == 1)
                        {
                                out.push_back(*prev);
                                prev = iter;
                        }
                        else
                        {
                                prev = iter;
                        }
                }
                iter++;
        }
}
// T for simple type
template <class T> void reduplicate(std::vector<T>& in, std::vector<T>& out)
{
        std::sort(in.begin(), in.end());
        std::set<T> tmp;
        for (auto item : in)
        {
                tmp.insert(item);
        }
        out.clear();
        out.resize(0);
        for (auto item : tmp)
        {
                out.push_back(item);
        }
        std::sort(out.begin(), out.end());
}

void transcriptomeFa::mkdb()
{
        //check if each transcriptome is longer than 256
        uint32_t sum      = 0;
        uint32_t count150 = 0;
        uint32_t count256 = 0;

        for (auto iter : this->originFa)
        {
                sum += iter.length;
                if (iter.length < 128)
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

        for (auto rna : this->originFa)
        {
                addKmer(rna);
        }

        for (uint64_t i = 0; i < gene_T.size(); i++)
        {
                std::vector<uint64_t> tmp;
                reduplicate<uint64_t>(gene_T[i], tmp);
                for (auto item : tmp)
                {
                        kmer_t tmpK;
                        tmpK.kmer   = item;
                        tmpK.geneId = i;
                        geneKmerTable.push_back(tmpK);
                }
        }
        sort_kmer(geneKmerTable, thread);
        redupKmer(geneKmerTable, contact);
        std::cout << contact.size() << std::endl;
}

void transcriptomeFa::mkindex()
{
        uint32_t shiftSize = (sizeof(uint64_t) - sizeof(uint8_t)) * 8;
        uint32_t index_stash[256];
        for (uint32_t i=0;i<256;i++){
                index_stash[i]=0;
                index[i]=0;
        }
        for (uint32_t i = 0; i < contact.size(); i++)
        {
                uint64_t s = contact[i].kmer;
                s          = s >> shiftSize;
                index_stash[s]++;
        }
        //do prefix sum
        uint32_t index_presum[256];
        index[0] = 0;
        //index[1] = index_stash[0];
        for (uint32_t i = 1; i < 256; i++)
        {
                //index[i] += index_stash[i - 1];
                index[i] = index[i-1] + index_stash[i-1];
        }
}

void transcriptomeFa::writedb()
{
        FILE* fdb;
        FILE* findex;
        FILE* fgene;

        const std::string f_db    = "SCQRDB.DATA";
        const std::string f_index = "SCQRDB.INDEX";
        const std::string f_gene  = "SCQRDB.GENE";

        fdb = fopen(f_db.c_str(), "wb");
        std::fwrite(contact.data(), sizeof(kmer_t), contact.size(), fdb);
        std::fclose(fdb);
        findex = fopen(f_index.c_str(), "wb");
        std::fwrite(index, sizeof(index[0]), 256, findex);
        std::fclose(findex);

        fgene           = fopen(f_gene.c_str(), "w");
        uint32_t geneId = 0;
        for (auto geneName : geneList)
        {
                std::fprintf(fgene, "%d\t%s\n", geneId, geneName.c_str());
                geneId++;
        }
        std::fclose(fgene);
}

void transcriptomeFa::clear()
{
        for (auto item : originFa)
        {
                delete[] item.seq;
        }
}

typedef struct
{
        std::vector<kmer_t>      db;
        uint32_t                 index[256];
        std::vector<std::string> gene;
} loadedDB;

void read_DB_Index(const std::string         f_db,
                   const std::string         f_index,
                   const std::string         f_gene,
                   loadedDB & ldDB)
{
        FILE* fdb;
        FILE* findex;
        FILE* fgene;

        fdb                = std::fopen(f_db.c_str(), "rb");
        uint32_t db_length = 0;
        std::fseek(fdb, 0, SEEK_END);
        db_length = std::ftell(fdb) / sizeof(kmer_t);
        ldDB.db.resize(db_length);
        std::fseek(fdb, 0, SEEK_SET);
        std::fread(ldDB.db.data(), sizeof(kmer_t), db_length, fdb);
        std::fclose(fdb);

        findex                = std::fopen(f_index.c_str(), "rb");
        uint32_t index_length = 256;
        std::fread(ldDB.index, sizeof(uint32_t), index_length, findex);
        std::fclose(findex);

        fgene = std::fopen(f_gene.c_str(), "r");
        std::vector<std::string> _gene;
        while (true)
        {
                char     geneName_c[128];
                if (fscanf(fgene, "%s\n", geneName_c) != EOF)
                {
                        _gene.push_back(std::string(geneName_c));
                        //ldDB.gene.push_back(std::string(geneName_c));
                }
                else{
                        break;
                }
        }
        for (uint32_t i=0;i<_gene.size()/2;i++){
                ldDB.gene.push_back(_gene[2*i+1]);
        }
        //fmt::print("{}\n", ldDB.gene);
        std::fclose(fgene);
};
