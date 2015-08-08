/*
 * transmap projection of annotations
 */
#ifndef transMap_hh
#define transMap_hh
#include "jkinclude.hh"
#include <string>
#include <map>
#include <vector>
using namespace std;

/*
 * vector of psls
 */
typedef vector<struct psl*> PslVector;

/* set of mapped alignments */
class PslMapping {
    public:
    struct psl* fSrcPsl;
    PslVector fMappedPsls;  // top one is best scoring
    int fScore;  // mapping score of top psl; 0 is perfect

    private:
    static int numPslAligned(struct psl* psl);
    void sortMappedPsls();

    public:
    /* constructor, sort mapped PSLs */
    PslMapping(struct psl* srcPsl,
               PslVector& mappedPsls);

    /* free up psls */
    ~PslMapping();

    /* Compute a mapping score between the src and mapped psl.  A perfect
     * mapping is a zero score.  Extra inserts count against the score. */
    static int calcPslMappingScore(struct psl* srcPsl,
                                   struct psl* mappedPsl);
};


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
    
    /* Map a single input PSL and return a list of resulting mappings.
     * Keep PSL in the same order, even if it creates a `-' on the target. */
    PslMapping* mapPsl(struct psl* inPsl) const;

};

#endif
