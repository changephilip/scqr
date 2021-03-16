#include "../include/querydb.h"
#include <getopt.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
        if (argv[1] != NULL and std::atoi(argv[2]))
        {
                //        transcriptomeFa t(argv[1], 8);
        }
        else
        {
                exit(EXIT_FAILURE);
        }

        //std::vector<kmer_t> db;
        //uint32_t index[256];
        //std::vector<std::string> gene;
        loadedDB lldb;
        read_DB_Index("SCQRDB.DATA", "SCQRDB.INDEX", "SCQRDB.GENE", lldb);
        rFirst r1(argv[3]);
        rSecond r2(argv[4],r1);
        rQuery rQ(r1,r2,lldb,32,argv[5]);
        return 0;
}
