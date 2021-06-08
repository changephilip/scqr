#ifndef __READ
#define __READ
#include "scqr.h"
#include "mkdb.h"
#include "readR1.h"
#include "readR2.h"
//#include "../include/readR2.h"

//#include <parallel_hashmap/phmap.h>

class rRead {
       public:
        uint32_t thread;

        rRead(const std::string r1gz, const std::string r2gz, uint32_t _thread);
        void    CountBarcode(barcode_t _bc);
        rFirst  rf;
        rSecond rs;
};

void rRead::CountBarcode(barcode_t _bc)
{
        auto it = rf.barcodeCount.find(_bc);
        if (it == rf.barcodeCount.end())
        {
                rf.barcodeCount.insert(std::pair<barcode_t, uint32_t>(_bc, 1));
        }
        else
        {
                it->second++;
        }
}

rRead::rRead(const std::string r1gz, const std::string r2gz, uint32_t _thread)
{
        this->thread = _thread;
        gzFile  fp1;
        gzFile  fp2;
        kseq_t *seq1;
        kseq_t *seq2;

        rf.thread = _thread;

        std::time_t t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Loading Reads Start "
                  << std::endl;

        fp1             = gzopen(r1gz.c_str(), "r");
        fp2             = gzopen(r2gz.c_str(), "r");
        seq1            = kseq_init(fp1);
        seq2            = kseq_init(fp2);
        uint32_t readId = 0;
        //while (kseq_read(seq1) >= 0 and kseq_read(seq2) >= 0)
        while (kseq_read(seq1) >= 0)
        {
                {
                        //std::string s1 = seq1->name.s;
                        //std::string s2 = seq2->name.s;
                        if (false)
                        //if (s1 != s2)
                        {
                                std::perror("Please use fastqc to check origin file!\n");
                                //std::printf("SEQ1\n%s\n", s1.c_str());
                                //std::printf("SEQ2\n%s\n", s2.c_str());
                                exit(EXIT_FAILURE);
                        }
                        else
                        {
                                barcode_t thisbarcode = readBarcode(seq1->seq.s);
                                umi_t     thisumi     = readUmi(seq1->seq.s);
                                //CountBarcode(thisbarcode);
                                rf.barcodeVector.push_back(thisbarcode);

                                rOneRead tmpRead;
                                tmpRead.umi = thisumi;
                                tmpRead.id  = readId;
                                //rf.barcodeOfRead.push_back(thisbarcode);
                                rf.reads.push_back(tmpRead);
                                //rf.barcodeSet.insert(thisbarcode);
                                /*
                                rSecondRead tmpRead2;
                                tmpRead2.sample = thisbarcode;
                                tmpRead2.seq    = seq2->seq.s;
                                reads.emplace_back(tmpRead2);
                                */
                                readId++;
                        }
                }
        }

        rf.barcodeOfRead =
                std::vector<barcode_t>(rf.barcodeVector.begin(), rf.barcodeVector.end());

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Loading Reads R1 Completed "
                  << std::endl;

        std::cout << "Num Of Reads is " << rf.barcodeOfRead.size() << std::endl;

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

        //std::unique(rf.barcodeVector.begin(), rf.barcodeVector.end());
        //rf.barcodeSet = scqr_set<barcode_t>(rf.barcodeVector.begin(),rf.barcodeVector.end());

        rf.barcodeCorrect();
        rf.umiDuplicate();

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Correcting Barcode Completed "
                  << std::endl;

        uint64_t numOftaggedReads = 0;

        numOftaggedReads =
                std::accumulate(rf.cellReadsCount.begin(), rf.cellReadsCount.end(), 0);

        std::cout << "Num Of Tagged Reads is " << numOftaggedReads << std::endl;

        rs.readsPack.resize(rf.correctedBarcodeVector.size());

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Read R2 Started " << std::endl;

        //gzclose(fp2);
        //fp2 = gzopen(r2gz.c_str(),"r");
        //seq2= kseq_init(fp2);
        uint32_t readCount = 0;
        while (kseq_read(seq2) >= 0)
        {
                if (rf.reads[readCount].sample != BARCODE_NULL)
                {
                        rSecondRead tmp;
                        tmp.seq = seq2->seq.s;

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
        /*
        for (uint32_t i = 0; i < reads.size(); i++)
        {
                if (rf.reads[i].sample != 0xFFFFFFFF)
                {
                        //tmp.sample = rf.reads[i].sample;
                        rs.readsPack[rf.reads[i].sample].emplace_back(reads[i]);
                }
        }
        */

        t = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&t)) << "Generate reads pack Completed "
                  << std::endl;
}
#endif
