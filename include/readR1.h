#ifndef  __READR1
#define __READR1
#include "scqr.h"
#include "mkdb.h"

scqr_map<uint32_t,uint32_t> ss;
scqr_set<uint32_t> consT;

class rFirst {
       public:
        uint32_t thread;
        /*a unique identifier, tag=[barcode]+[umi]*/
        std::vector<rOneRead> reads;
        scqr_set<barcode_t>   barcodeSet;
        /*assign 0-based id for each sample*/
        scqr_map<barcode_t, uint32_t> correctedBarcodeMapList;
        /*get origin sampleId(barcode) with index*/
        std::vector<barcode_t>        correctedBarcodeVector;
        scqr_map<barcode_t, uint32_t> barcodeCount;
        std::vector<uint32_t>         cellReadsCount;
        /*store pair1 string for R2 search for the paired reads R1 readId*/
        scqr_map<std::string, uint32_t> readsNameTable;
        std::vector<barcode_t>       barcodeOfRead;
        rFirst(const std::string r1gz, uint32_t _labels, uint32_t _thread);
        rFirst();
        std::hash<std::string> stringHash;
        uint32_t               labels;
        void
             exchange(std::vector<barcodeCount_t> &in, std::vector<uint32_t> &top, uint32_t i);
        void transform(std::vector<barcodeCount_t> &in,
                       std::vector<uint32_t> &      src,
                       std::vector<uint32_t> &      dst,
                       uint32_t                     i);
        void removeDuplicateUMI(std::vector<rOneRead> &in,uint32_t start,uint32_t end);
        void barcodeCorrect();
        void indelCorrect();
        void umiDuplicate();
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
        return (umi_t)baseToBinaryForward(seq + barcode + jumpbit, umi);
}

inline tag_t readTag(const char *seq)
{
        return (tag_t)baseToBinaryForward(seq + jumpbit, umi + barcode);
}

/*
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
*/

/*
void rFirst::barcodeCorrect(){
        auto dbscan = DBSCAN<uint64_t, uint64_t>();

        auto data = std::vector<uint64_t>{barcodeSet.begin(),barcodeSet.end()};

        dbscan.Run(&data, 1, 1, 4000);
        auto noise = dbscan.Noise;
        auto clusters = dbscan.Clusters;
}
*/
inline uint64_t shift2bit(barcode_t &in, uint32_t x)
{
        return ((in >> (62 - x * 2)) & 3U);
}
inline bool compareBarcode(barcode_t &a, barcode_t &b, uint32_t r)
{
        uint32_t d = 0;
        for (uint32_t i = 0; i < 32; i++)
        {
                //if (_mm_popcnt_u64(shift2bit(a, i) xor shift2bit(b, i)) > 0)
                if ((shift2bit(a, i) xor shift2bit(b, i)) != 0)
                {
                        d++;
                        if (d > r)
                        {
                                return false;
                        }
                }
        }
        if (d > r)
        {
                return false;
        }
        return true;
}

void rFirst::exchange(std::vector<barcodeCount_t> &in,
                      std::vector<uint32_t> &      top,
                      uint32_t                     i)
{
        for (uint32_t j = top.size() - 1; j > i + 1; j--)
        {
                if (compareBarcode(in[top[i]]._sample, in[top[j]]._sample,4))
                {
                        {
                                if (in[top[i]]._count <= in[top[j]]._count)
                                {
                                        in[top[i]]._sample = in[top[j]]._sample;
                                        in[top[i]]._count  = in[top[j]]._count;
                                        return;
                                }
                        }
                }
        }
}

void rFirst::transform(std::vector<barcodeCount_t> &in,
                       std::vector<uint32_t> &      src,
                       std::vector<uint32_t> &      dst,
                       uint32_t                     i)
{
        //for (uint32_t j = dst.size() - 1; j > i + 1; j--)
        //for (uint32_t j = dst.size() - 1; j >= 0; j--)
        for (uint32_t j = 0; j < dst.size(); j++)
        {
                if (compareBarcode(in[src[i]]._sample, in[dst[j]]._sample,4))
                {
                        in[src[i]]._sample = in[dst[j]]._sample;
                        in[src[i]]._count  = in[dst[j]]._count;
                        return;
                }
        }
        in[src[i]]._count = 0;
}
/*
  probability of error barcode is very low and real barcode has most of reads count.
*/
void rFirst::barcodeCorrect()
{
        omp_set_dynamic(false);
        omp_set_num_threads(this->thread);
        /*
        int                                  N = barcodeSet.size();
        std::vector<std::array<uint64_t, 1>> data;
        for (auto item : barcodeSet)
        {
                data.push_back(std::array<uint64_t, 1>{item});
        }
        */
        std::vector<barcodeCount_t> barcodeCountVector;
        for (auto item : barcodeCount)
        {
                barcodeCount_t t;
                t.sample  = item.first;
                t._sample = t.sample;
                t.count   = item.second;
                t._count  = t.count;
                barcodeCountVector.emplace_back(t);
        }
        __gnu_parallel::sort(barcodeCountVector.begin(),
                             barcodeCountVector.end(),
                             [](barcodeCount_t left, barcodeCount_t right) {
                                     return left.count < right.count;
                             });
        /*
        std::sort(barcodeCountVector.begin(),
                  barcodeCountVector.end(),
                  [](barcodeCount_t &left, barcodeCount_t &right) {
                          return left.count < right.count;
                  });
        FILE *f;
        f = std::fopen("R1_Statistic.txt", "w");
        for (auto item : barcodeCountVector)
        {
                std::fprintf(f, "%lu\t%d\n", item.sample, item.count);
        }
        std::fclose(f);
        */

        std::vector<uint32_t> top;
        std::vector<uint32_t> buttom;
        uint32_t              threshold = 128;
        // split all barcode into top(>128) and buttom (<128)
        // after sort barcodeCountVector, if i1 < i2 then i1._count i2._count
        for (uint32_t i = 0; i < barcodeCountVector.size(); i++)
        {
                if (barcodeCountVector[i].count > threshold)
                {
                        top.push_back(i);
                }
                else
                {
                        buttom.push_back(i);
                }
        }

        //find duplicate in top,not recommended to run in parallel
        //#pragma omp parallel for
        for (uint32_t i = 0; i < top.size(); i++)
        {
                exchange(barcodeCountVector, top, i);
        }

        //merge duplicate into compact,step 1
        std::vector<uint32_t> compact;
        for (uint32_t i = 0; i < top.size(); i++)
        {
                if (barcodeCountVector[top[i]]._sample ==
                    barcodeCountVector[top[i]].sample)
                {
                        compact.push_back(top[i]);
                }
        }
        // step 2
        for (uint32_t i=0;i < top.size(); i++){
                if (barcodeCountVector[top[i]]._sample != barcodeCountVector[top[i]].sample){
                        transform(barcodeCountVector, top, compact, i);
                }
        }
        // check
        
        //find and merge duplicate from buttom into compact
        //rewrite its _sample with _sample from compact
        #pragma omp parallel for
        for (uint32_t i = 0; i < buttom.size(); i++)
        {
                transform(barcodeCountVector, buttom, compact, i);
        }
        //check homeless barcode
        std::vector<barcode_t> cellLess;
        uint32_t               cellLessCount = 0;
        for (uint32_t i = 0; i < buttom.size(); i++)
        {
                if (barcodeCountVector[buttom[i]]._count == 0)
                {
                        cellLess.push_back(barcodeCountVector[buttom[i]].sample);
                        cellLessCount += barcodeCountVector[buttom[i]].count;
                }
        }

        //generate the map from origin sample to corrected sample
        std::map<barcode_t, barcode_t> reMapToBarcode;
        for (auto item : barcodeCountVector)
        {
                if (item._count != 0)
                {
                        reMapToBarcode.insert(std::pair<barcode_t, barcode_t>(
                                item.sample, item._sample));
                }
        }
        // generate cell set
        scqr_set<barcode_t> compactSet;
        for (auto i : compact)
        {
                compactSet.insert(barcodeCountVector[i].sample);
        }

        // assign a cell id for each cell start from 0
        auto     it     = compactSet.begin();
        uint32_t cellId = 0;
        while (it != compactSet.end())
        {
                correctedBarcodeMapList.insert(
                        std::pair<barcode_t, uint32_t>(*it, cellId));
                it++;
                cellId++;
                correctedBarcodeVector.push_back(*it);
        }
        //assgin barcode id for each read
#pragma omp parallel for
        for (uint32_t i = 0; i < reads.size(); i++)
        {
                auto it = reMapToBarcode.find(barcodeOfRead[i]);
                if (it != reMapToBarcode.end())
                {
                        reads[i].sample =
                                correctedBarcodeMapList.find(it->second)->second;
                }
                else{
                        reads[i].sample = 0xFFFFFFFF;  // cellLess read
                }
        }

        //calculate reads in every cell
        //cellReadsCount.resize(compact.size());
        cellReadsCount = std::vector<uint32_t>(compact.size(), 0);
        for (barcodeCount_t item : barcodeCountVector)
        {
                if (item._count != 0)
                {
                        cellReadsCount[correctedBarcodeMapList.find(item._sample)
                                               ->second] += item._count;
                }
        }
}

void rFirst::indelCorrect(){
        std::vector<barcode_t>  target;
        for (auto item: correctedBarcodeMapList){
                target.push_back(item.first);
        }

        scqr_set<uint64_t> miniKmer;
        uint32_t miniKmerSize= 8;

        for (auto item: correctedBarcodeMapList){
                
        }
}

void rFirst::removeDuplicateUMI(std::vector<rOneRead> &in, uint32_t start, uint32_t end){
        scqr_set<umi_t> umiSet;
        for (uint32_t i=start;i<end;i++){
                auto it = umiSet.find(in[i].umi);
                if (it ==umiSet.end()){
                        umiSet.insert(i);
                }
                else{
                        in[i].sample = 0xFFFFFFFF;
                }
        }
}

void rFirst::umiDuplicate(){
        __gnu_parallel::sort(reads.begin(), reads.end(),[](rOneRead left, rOneRead right){return left.sample < right.sample;});
        auto p = reads.begin();
        auto n = reads.begin()++;
        std::vector<uint32_t> s;
        std::vector<uint32_t> e;
        while (n!= reads.end()){
                if (n->sample == p->sample){
                        n++;
                }
                else{
                        s.push_back(std::distance(reads.begin(), p));
                        e.push_back(std::distance(reads.begin(), n));
                        p = n;
                        n++;
                }
        }
        omp_set_dynamic(false);
        omp_set_num_threads(this->thread);
        #pragma  omp parallel for
        for (uint32_t i =0 ;i< s.size();i++){
                //__gnu_parallel::sort(reads.begin() + s[i], reads.begin()+e[i],[](rOneRead &left,rOneRead &right){return left.umi < right.umi;});
                removeDuplicateUMI(this->reads, s[i], e[i]);
        }

}

/*
  filter unique tag(barcode + umi) and get these reads from R1
  save out these sample(barcode) id and transfer into 0-based index.
  How to find pair-end reads?
  construct readsNameTable with pair of (name.s , readId)
*/
rFirst::rFirst(const std::string r1gz, uint32_t _labels, uint32_t _thread)
{
        this->thread = _thread;
        this->labels = _labels;
        gzFile  fp;
        kseq_t *seq;
        int     l;
        fp              = gzopen(r1gz.c_str(), "r");
        seq             = kseq_init(fp);
        uint32_t readId = 0;

        std::vector<std::pair<std::string, uint32_t>> readsNameVector;
        std::time_t                                t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tStart of R1 " << std::endl;

        while ((l = kseq_read(seq)) >= 0)
        {
                barcode_t thisbarcode = readBarcode(seq->seq.s);
                umi_t     thisumi     = readUmi(seq->seq.s);
                auto it = barcodeCount.find(thisbarcode);
                if (it != barcodeCount.end())
                {
                        it->second++;
                }
                else
                {
                        barcodeCount.insert(
                                std::pair<barcode_t, uint32_t>(thisbarcode, 1));
                }
                //uint32_t thistag     = ((uint32_t)thisbarcode << 12) || (0x0u || thisumi);
                {
                        rOneRead tmpRead;
                        tmpRead.umi = thisumi;
                        //tmpRead.sample = thisbarcode;
                        barcodeOfRead.push_back(thisbarcode);
                        //tmpRead.id  = readId;
                        //tmpRead.sample = thisbarcode;
                        //tmpRead.pair1 = seq->name.s;
                        reads.push_back(tmpRead);
                        barcodeSet.insert(thisbarcode);
                        /*
                        readsNameTable.insert(
                                std::pair<std::string, uint32_t>(seq->name.s, readId));
                        */
                        readsNameVector.push_back(std::make_pair(seq->name.s, readId));
                        //readsNameVector.push_back(std::make_pair(std::hash<std::string>{}(seq->name.s),readId));
                        readsNameVector.push_back(
                                std::make_pair(seq->name.s, readId));
                        readId++;
                }
                //readId++;
        }

        readsNameTable = scqr_map<std::string, uint32_t>(readsNameVector.begin(),
                                                      readsNameVector.end());

        readsNameVector.resize(0);
        readsNameVector.clear();

        std::cout << reads.size() << std::endl;
        std::cout << barcodeSet.size() << std::endl;
        // reads.size() >> seqId.size()
        //assert(reads.size() == seqId.size());
        barcodeCorrect();
        barcodeOfRead.resize(0);

        umiDuplicate();
        //for (auto &read : reads)
        //{
        //        read.sample = sampleList.find(read.tag >> 12)->second;
        //}

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "\tEnd of R1 " << std::endl;
        scqr_set<uint32_t> mm;
}
#endif
