#ifndef __QUERYDB
#define __QUERYDB

#include "read.h"
#include "readR2.h"
#include <mutex>
#include <ctime>
typedef uint32_t gene_t;
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
        quantRNAMatrix = new uint32_t *[r1.correctedBarcodeVector.size()];
        for (uint32_t i = 0; i < r1.correctedBarcodeVector.size(); i++)
        {
                quantRNAMatrix[i] = new uint32_t[lddb.gene.size()];
                for (uint32_t j=0;j< lddb.gene.size();j++){
                        quantRNAMatrix[i][j]=0;
                }
        }
        std::cout << "GENE LIST LENGTH is " << lddb.gene.size() << std::endl;
        std::cout << "Num of Barcode(sample) is " << r1.correctedBarcodeVector.size() << std::endl;
        omp_set_dynamic(false);
        omp_set_num_threads(_thread);
        #pragma omp parallel for
        for (uint32_t i = 0; i < r1.correctedBarcodeVector.size(); i++)
        {
                for (uint32_t j = 0; j < r2.readsPack[i].size(); j++)
                {
                        uint32_t q = align(r2.readsPack[i][j]);
                        if (q < lddb.gene.size())
                        {
                                quantRNAMatrix[i][q]++;
                        }
                }
                if ( i % 32 ==0 and true){
                        std::time_t t = std::time(nullptr);
                        std::cout << std::asctime(std::localtime(&t)) << "\tStep is " << i << std::endl;
                }
        }


        //print result
        FILE* fresult;
        fresult=std::fopen(result.c_str(),"w");

        //header
        std::fprintf(fresult,"G/S");
        for (uint32_t i=0;i< r1.correctedBarcodeVector.size();i++){
                std::fprintf(fresult, "\t%lx",r1.correctedBarcodeVector[i]);
        }
        std::fprintf(fresult,"\n");
        std::fprintf(fresult,"RC");
        for (uint32_t i=0;i< r1.correctedBarcodeVector.size();i++){
                std::fprintf(fresult, "\t%d",r1.cellReadsCount[i]);
        }
        std::fprintf(fresult, "\n");

        for (uint32_t i=0;i < lddb.gene.size();i++){
                std::fprintf(fresult, "%s",lddb.gene[i].c_str());
                for (uint32_t j=0;j< r1.correctedBarcodeVector.size();j++){
                        std::fprintf(fresult, "\t%d",quantRNAMatrix[j][i]);
                }
                std::fprintf(fresult, "\n");
        }
        //print only sample 1
        /*
        fprintf(fresult, "%d\t",r1.sampleId[0]);
        auto thissampleid = r1.sampleId[0];
        for (uint32_t i =0 ;i< lddb.gene.size();i++){
                fprintf(fresult, "%d\n",quantRNAMatrix[0][i]);
        }
        */
        std::fclose(fresult);

        for (uint32_t i=0;i< r1.correctedBarcodeVector.size();i++){
                delete[] quantRNAMatrix[i];
        }
        delete[] quantRNAMatrix;
 }

template <class T1, class T2> inline T1 findMapMax(scqr_map<T1, T2> &in)
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
        if (m->first < 32){return 0xFFFFFFF;}
        return m->first;
}

inline uint32_t rQuery::align(const rSecondRead &read)
{
        uint32_t                     stride = 1;
        if (read.seq.size() < kmerLength) return 0xFFFFFFFF;
        uint32_t                     n_kmer = read.seq.size() - kmerLength + 1;
        //score : [geneId, count]
        scqr_map<uint32_t, uint16_t> score;
        //auto                         p        = read.seq.begin();
        const char *t = read.seq.c_str();
        uint32_t                     distance = 0;
        for (uint32_t i =0;i < n_kmer ;i++)
        {
                strainType thisStrainType;
                if (forwardOrBackward(t, kmerLength))
                {
                        thisStrainType = forward;
                }
                else
                {
                        thisStrainType = backward;
                }
                uint64_t kmer = baseToBinary(t, kmerLength, thisStrainType);
                uint64_t key  = search(kmer, mydb);
                t++;
                if (key == 0xFFFFFFFF) continue;
                if (score.find(key) != score.end())
                {
                        score.find(key)->second++;
                }
                else
                {
                        score.insert(std::pair<uint32_t, uint16_t>(key, 0));
                }
                distance++;
                //auto p = score.begin();
        }
        if (score.size()==0){
                return 0xFFFFFFFF;
        }
        return findMapMax(score);
}

inline uint32_t rQuery::search(uint64_t key, const loadedDB &lddb)
{
        uint32_t start;
        uint32_t end;
        //indexSize = 65536 => 64 - 16
        //indexSize = 256 => 64 -8
        start = lddb.index[key >> (64 - 16)];
        if ((key >> 48) != 0xFFFF)
        {
                end = lddb.index[(key >> (48)) + 1];
        }
        else
        {
                end = lddb.db.size();
        }

        uint32_t ps = start;
        uint32_t pe = end;
        uint32_t mid = ps + (pe-ps)/2;
       
        while ( ps + 1 < pe)
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
        if ( lddb.db[ps].kmer == key){
                return ps;
        }
        if (lddb.db[pe].kmer == key){
                return pe;
        }
        return 0xFFFFFFFF;
}

#endif
