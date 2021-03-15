#include "../include/readR2.h"
#include <mutex>

class rQuery {
       public:
        rQuery(const rFirst &r1, const rSecond &r2, const loadedDB &lddb, uint32_t _thread, const std::string &result);
        std::vector<uint32_t> rnaQ;
        loadedDB              mydb;
        uint32_t **           quantRNAMatrix;
        uint32_t              align(const rSecondRead &read);
        uint32_t              search(uint64_t key, const loadedDB &lddb);
};

rQuery::rQuery(const rFirst &r1, const rSecond &r2, const loadedDB &lddb, uint32_t _thread, const std::string &result)
{
        mydb           = lddb;
        quantRNAMatrix = new uint32_t *[r1.sampleId.size()];
        for (uint16_t i = 0; i < r1.sampleId.size(); i++)
        {
                quantRNAMatrix[i] = new uint32_t[lddb.gene.size()];
        }
        omp_set_dynamic(false);
        omp_set_num_threads(_thread);
        #pragma omp parallel for
        for (uint16_t i = 0; i < r1.sampleId.size(); i++)
        {
                for (uint32_t j = 0; j < r2.readsPack[i].size(); j++)
                {
                        uint32_t q = align(r2.readsPack[i][j]);
                        quantRNAMatrix[i][q]++;
                }
        }


        //print result
        FILE* fresult;
        fresult=std::fopen(result.c_str(),"w");
        std::fprintf(fresult,"\t");
        //print only sample 1
        fprintf(fresult, "%d\t",r1.sampleId[0]);
        auto thissampleid = r1.sampleId[0];
        for (uint32_t i =0 ;i< lddb.gene.size();i++){
                fprintf(fresult, "%d\n",quantRNAMatrix[0][i]);
        }
 }

template <class T1, class T2> inline T1 findMapMax(std::map<T1, T2> &in)
{
        auto p = in.begin();
        auto m = in.begin();
        p++;
        while (p != in.end())
        {
                if (m->second < p->second)
                {
                        m = p;
                }
                p++;
        }
        return m->first;
}

inline uint32_t rQuery::align(const rSecondRead &read)
{
        uint32_t                     stride = 1;
        uint32_t                     n_kmer = read.seq.size() - kmerLength + 1;
        std::map<uint32_t, uint16_t> score;
        auto                         p        = read.seq.begin();
        uint32_t                     distance = 0;
        while (p != read.seq.end())
        {
                strainType thisStrainType;
                if (forwardOrBackward(p.base(), kmerLength))
                {
                        thisStrainType = forward;
                }
                else
                {
                        thisStrainType = backward;
                }
                uint64_t kmer = baseToBinary(p.base(), kmerLength, thisStrainType);
                uint64_t key  = search(kmer, mydb);
                if (key != 0xFFFFFFFF and score.find(key) != score.end())
                {
                        score.find(key)->second++;
                }
                else
                {
                        score.insert(std::pair<uint32_t, uint16_t>(key, 0));
                }
                distance++;
                p++;
                auto p = score.begin();
        }

        return findMapMax(score);
}

inline uint32_t rQuery::search(uint64_t key, const loadedDB &lddb)
{
        uint32_t start;
        uint32_t end;
        start = lddb.index[key >> (64 - 8)];
        if (start != 0xFF)
        {
                end = lddb.index[(key >> (64 - 8)) + 1];
        }
        else
        {
                end = lddb.db.size();
        }

        uint32_t ps = start;
        uint32_t pe = end;
        uint32_t mid = ps + (pe-ps)/2;
       
        while (ps != pe)
        {
                if (lddb.db[mid].kmer < key)
                {
                        ps  = mid;
                        mid = ps + (pe - ps) / 2;
                }
                else if (lddb.db[mid].kmer > key)
                {
                        pe  = mid;
                        mid = ps + (pe - ps) / 2;
                }
                else
                {
                        return lddb.db[mid].geneId;
                }
        }
        return 0xFFFFFFFF;
}
