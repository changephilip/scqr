#include "../include/read.h"
#include "../include/readcr.h"
#include "../include/querydb.h"
#include <getopt.h>
#include "../include/cxxopts.hpp"
#include <stdlib.h>

int main(int argc, char *argv[])
{
        /*
        if (argv[1] != NULL and std::atoi(argv[2]))
        {
                //        transcriptomeFa t(argv[1], 8);
        }
        else
        {
                exit(EXIT_FAILURE);
        }
        */

        /*
        int mode = std::atoi(argv[1]);
        if (argv[1] != NULL)
        {
                switch (mode)
                {
                case 1: {
                        transcriptomeFa t(argv[2], std::atoi(argv[3]));
                        break;
                }
                case 2: {
                        loadedDB lldb;
                        read_DB_Index("SCQRDB.DATA", "SCQRDB.INDEX", "SCQRDB.GENE", lldb);
                        rFirst  r1(argv[2],4000,32);
                        rSecond r2(argv[3], r1);
                        rQuery  rQ(r1, r2, lldb, 32, argv[4]);
                        break;
                }
                default:
                        exit(EXIT_FAILURE);
                }
        }
        */

        cxxopts::Options options("scqr", "Quantification of RNA in Single Cell");
        options.add_options()("m,mode",
                              "--mode=index for build Index or \n--mode=run for run "
                              "SCQR,\n--mode=pair-hide for sequence-sorted paired reads",
                              cxxopts::value<std::string>())(
                "t,thread", "Thread used", cxxopts::value<uint32_t>())(
                "f,fasta", "RNA fasta for index", cxxopts::value<std::string>())(
                "r,R1", "R1 fasta file ,fastq.gz", cxxopts::value<std::string>())(
                "l,R2", "R2 fasta file ,fastq.gz", cxxopts::value<std::string>())(
                "o,out", "Output expression matrix", cxxopts::value<std::string>())(
                "h,help", "Print usage");
        auto results = options.parse(argc, argv);
        while (results.count("mode") == 1)
        {
                if (results["mode"].as<std::string>() == "index")
                {
                        transcriptomeFa buildIndex(
                                results["fasta"].as<std::string>().c_str(),
                                results["thread"].as<uint32_t>());
                        break;
                }
                else if (results["mode"].as<std::string>() == "run")
                {
                        loadedDB lldb;
                        read_DB_Index("SCQRDB.DATA", "SCQRDB.INDEX", "SCQRDB.GENE", lldb);
                        rFirst  r1(results["R1"].as<std::string>(),
                                  4000,
                                  results["thread"].as<uint32_t>());
                        rSecond r2(results["R2"].as<std::string>(), r1);
                        rQuery  rQ(r1,
                                  r2,
                                  lldb,
                                  results["thread"].as<uint32_t>(),
                                  results["out"].as<std::string>());
                        break;
                }
                else if (results["mode"].as<std::string>() == "pair-hide")
                {
                        loadedDB lldb;
                        read_DB_Index("SCQRDB.DATA", "SCQRDB.INDEX", "SCQRDB.GENE", lldb);
                        rRead  pair_reads(results["R1"].as<std::string>(),
                                         results["R2"].as<std::string>(),
                                         results["thread"].as<uint32_t>());
                        rQuery rQ(pair_reads.rf,
                                  pair_reads.rs,
                                  lldb,
                                  results["thread"].as<uint32_t>(),
                                  results["out"].as<std::string>());
                        break;
                }
                else if (results["mode"].as<std::string>() == "cellRanger")
                {
                        loadedDB lldb;
                        read_DB_Index("SCQRDB.DATA", "SCQRDB.INDEX", "SCQRDB.GENE", lldb);
                        rCellRanger rCR(results["I1"].as<std::string>(),
                                        results["R1"].as<std::string>(),
                                        results["R2"].as<std::string>(),
                                        results["thread"].as<uint32_t>());
                        rQuery      rQ(rCR.rf,
                                  rCR.rs,
                                  lldb,
                                  results["thread"].as<uint32_t>(),
                                  results["out"].as<std::string>());
                        break;
                }
                else
                {
                        std::cout << options.help() << std::endl;
                        std::cout << "Example for index : scqr --mode=index grch38rna.fa "
                                     "--thread=8"
                                  << std::endl;
                        std::cout << "Example for run   : scqr --mode=run "
                                     "--R1=foo.R1.fastq.gz --R2=foo.R2.fastq.gz "
                                     "--out=matrix.txt --thread=8"
                                  << std::endl;
                        std::cout << "Example for run   : scqr --mode=cellRanger "
                                     "--I1=foo.I1.fastq.gz --R1=foo.R1.fastq.gz "
                                     "--R2=foo.R2.fastq.gz "
                                     "--out=matrix.txt --thread=8"
                                  << std::endl;

                        std::cout << "\n" << std::endl;
                        exit(EXIT_FAILURE);
                }
        }
        if (results.count("mode") != 1)
        {
                std::cout << options.help() << std::endl;
                std::cout
                        << "Example for index : scqr --mode=index grch38rna.fa --thread=8"
                        << std::endl;
                std::cout << "Example for run   : scqr --mode=run --R1=foo.R1.fastq.gz "
                             "--R2=foo.R2.fastq.gz --out=matrix.txt --thread=8"
                          << std::endl;
                std::cout << "Example for run   : scqr --mode=cellRanger "
                             "--I1=foo.I1.fastq.gz --R1=foo.R1.fastq.gz "
                             "--R2=foo.R2.fastq.gz "
                             "--out=matrix.txt --thread=8"
                          << std::endl;

                std::cout << "\n" << std::endl;
        }
        //std::vector<kmer_t> db;
        //uint32_t index[256];
        //std::vector<std::string> gene;
        return 0;
}
