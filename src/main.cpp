#include "../include/querydb.h"
#include <getopt.h>
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
                        rFirst  r1(argv[2],4000);
                        rSecond r2(argv[3], r1);
                        rQuery  rQ(r1, r2, lldb, 32, argv[4]);
                        break;
                }
                default:
                        exit(EXIT_FAILURE);
                }
        }
        //std::vector<kmer_t> db;
        //uint32_t index[256];
        //std::vector<std::string> gene;
        return 0;
}
