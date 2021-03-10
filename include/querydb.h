#include "../include/readR2.h"

class rQuery{
 public:
        rQuery(const rFirst& r1, const rSecond & r2, const loadedDB & lddb);
        std::vector<uint32_t> rnaQ;
};

rQuery::rQuery(const rFirst &r1,const rSecond &r2, const loadedDB & lddb ){
        rnaQ.resize(0);
        for (rSecondRead r2read:r2.reads){

        }
}
