#include "../include/mkdb.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
        if (argv[1] != NULL and std::atoi(argv[2]))
        {
                transcriptomeFa t(argv[1], 8);
        }
        else
        {
                exit(EXIT_FAILURE);
        }

        std::vector<kmer_t> db;
        uint32_t index[256];
        std::vector<std::string> gene;
        read_DB_Index("SCQRDB.DATA", "SCQRDB.INDEX", "SCQRDB.GENE", db, index, gene);
        return 0;
}
