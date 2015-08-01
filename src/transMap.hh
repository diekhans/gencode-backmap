/*
 * transmap projection of annotations
 */
#ifndef transMap_hh
#define transMap_hh
extern "C" {
#include "pslTransMap.h"
}
#include <string>
#include <map>
#include <vector>
using namespace std;

/*
 * vector of psls
 */
typedef vector<struct psl*> PslVector;

/*
 * transmap via alignment chains
 */
class TransMap {
    private:
    struct genomeRangeTree* fMapAlns;  // mapping alingments
    typedef map<const string, int> SizeMap;
    
    SizeMap fQuerySizes;   // query sequence sizes
    SizeMap fTargetSizes;  // target sequence sizes
    
    void mapAlnsAdd( struct psl *mapPsl);
    struct psl* chainToPsl(struct chain *ch,
                           bool swapMap);
    void loadMapChains(const string& chainFile,
                       bool swapMap);
    void addSeqSize(const string& seqName,
                    int seqSize,
                    SizeMap& sizeMap);
    struct psl* mapPslPair(struct psl *inPsl, struct psl *mapPsl) const;

    public:
    /* constructor, loading chains */
    TransMap(const string& chainFile,
             bool swapMap);

    /* destructor */
    ~TransMap() {
        // just memory leak and let exit() clean up
    }

    /* get the size of the a mapping query sequence */
    int getQuerySeqSize(const string& qName) const {
        return fQuerySizes.at(qName);
    }
    
    /* get the size of the a mapping target sequence */
    int getTargetSeqSize(const string& tName) const {
        return fTargetSizes.at(tName);
    }
    
    /* Map a single input PSL and return a list of resulting mappings */
    PslVector mapPsl(struct psl* inPsl) const;

};

#endif
