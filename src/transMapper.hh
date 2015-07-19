/*
 * transmap projection of annotations
 */
#ifndef transMapper_hh
#define transMapper_hh
extern "C" {
#include "pslTransMap.h"
}
#include <string>
using namespace std;

/*
 * transmap via alignment chains
 */
class TransMapper {
    private:
    struct genomeRangeTree* fMapAlns;  // mapping alingments
    void mapAlnsAdd( struct psl *mapPsl);
    struct psl* chainToPsl(struct chain *ch,
                           bool swapMap);
    void loadMapChains(const string& chainFile,
                       bool swapMap);
    struct psl* mapPslPair(struct psl *inPsl, struct psl *mapPsl);

    public:
    /* constructor, loading chains */
    TransMapper(const string& chainFile,
                bool swapMap);

    /* destructor */
    ~TransMapper() {
        // just memory leak and let exit() clean up
    }

    /* Map a single input PSL and return a list of resulting mappings */
    struct psl* mapPsl(struct psl* inPsl);

};

#endif
