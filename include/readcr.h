#ifndef __readcr
#define __readcr
#include "scqr.h"
#include "mkdb.h"
#include "readR1.h"
#include "readR2.h"

class rCellRanger {
       public:
        uint32_t thread;
        rCellRanger(const std::string barcgz,
                    const std::string umigz,
                    const std::string readgz,
                    uint32_t          _thread);
        rFirst  rf;
        rSecond rs;
};

rCellRanger::rCellRanger(const std::string barcgz,
                         const std::string umigz,
                         const std::string readgz,
                         uint32_t          _thread)
{
        this->thread = _thread;
        gzFile fpBarc;
        gzFile fpUmi;
        gzFile fpRead2;

        kseq_t *seqbarc;
        kseq_t *seqUmi;
        kseq_t *seqRead2;

        rf.thread = _thread;

        fpBarc  = gzopen(barcgz.c_str(), "r");
        fpUmi   = gzopen(umigz.c_str(), "r");
        fpRead2 = gzopen(readgz.c_str(), "r");

        kseq_init(fpBarc);
        kseq_init(fpUmi);
        kseq_init(fpRead2);

        uint32_t readId = 0;

        std::time_t t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Loading Reads Start "
                  << std::endl;

        while (kseq_read(seqbarc) >= 0 and kseq_read(seqUmi) >= 0)
        {
                barcode_t thisbarcode = readBarcode(seqbarc->seq.s);
                umi_t     thisumi     = readUmi(seqUmi->seq.s);
                rf.barcodeVector.push_back(thisbarcode);
                rOneRead tmpRead;
                tmpRead.umi = thisumi;
                tmpRead.id  = readId;
                rf.reads.push_back(tmpRead);
                readId++;
        }

        rf.barcodeOfRead =
                std::vector<barcode_t>(rf.barcodeVector.begin(), rf.barcodeVector.end());

        __gnu_parallel::sort(
                rf.barcodeVector.begin(),
                rf.barcodeVector.end(),
                [](barcode_t left, barcode_t right) { return left < right; });
        std::vector<std::pair<barcode_t, uint32_t>> uniqueCount;
        uint32_t                                    p = 0;
        uint32_t                                    n = 1;

        while (n <= rf.barcodeVector.size())
        {
                if (rf.barcodeVector[p] != rf.barcodeVector[n])
                {
                        uniqueCount.push_back(std::pair<barcode_t, uint32_t>(
                                rf.barcodeVector[p], n - p));
                        p = n;
                }
                n++;
        }

        rf.barcodeVector.clear();
        rf.barcodeCount =
                scqr_map<barcode_t, uint32_t>(uniqueCount.begin(), uniqueCount.end());
        uniqueCount.clear();

        rf.barcodeCorrect();
        rf.umiDuplicate();

        uint64_t numOftaggedReads = 0;

        numOftaggedReads =
                std::accumulate(rf.cellReadsCount.begin(), rf.cellReadsCount.end(), 0);

        std::cout << "Num Of Tagged Reads is " << numOftaggedReads << std::endl;
        rs.readsPack.resize(rf.correctedBarcodeVector.size());

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Read R2 Started " << std::endl;

        uint32_t readCount = 0;
        while (kseq_read(seqRead2) >= 0)
        {
                if (rf.reads[readCount].sample != BARCODE_NULL)
                {
                        rSecondRead tmp;
                        tmp.seq = seqRead2->seq.s;

                        rs.readsPack[rf.reads[readCount].sample].emplace_back(tmp);
                }
                readCount++;
        }
#ifdef DEBUG
        uint32_t readsPackSum = 0;
        for (uint32_t i = 0; i < rs.readsPack.size(); i++)
        {
                readsPackSum += rs.readsPack[i].size();
        }

        std::cout << "readsPackSum\t" << readsPackSum << std::endl;
#endif

        assert(readCount == rf.barcodeOfRead.size());

        rf.barcodeOfRead.resize(0);

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Generate reads pack Completed "
                  << std::endl;

}

#endif
