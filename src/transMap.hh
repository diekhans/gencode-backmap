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
    static int numPslAligned(const struct psl* psl);
    void sortMappedPsls();

    public:
    /* constructor, sort mapped PSLs */
    PslMapping(struct psl* srcPsl,
               PslVector& mappedPsls);

    /* free up psls */
    ~PslMapping();

    /* Compute a mapping score between the src and mapped psl.  A perfect
     * mapping is a zero score.  Extra inserts count against the score. */
    static int calcPslMappingScore(const struct psl* srcPsl,
                                   const struct psl* mappedPsl);
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
    ~TransMap();

    /* do we have a mapping query sequence */
    bool haveQuerySeq(const string& qName) const {
        return fQuerySizes.find(qName) != fQuerySizes.end();
    }
    
    /* do we have a mapping target sequence */
    bool haveTargetSeq(const string& tName) const {
        return fTargetSizes.find(tName) != fTargetSizes.end();
    }
    
    /* get the size of a mapping query sequence */
    int getQuerySeqSize(const string& qName) const {
        return fQuerySizes.at(qName);
    }
    
    /* get the size of a mapping target sequence */
    int getTargetSeqSize(const string& tName) const {
        return fTargetSizes.at(tName);
    }
    
    /* Map a single input PSL and return a list of resulting mappings.
     * Keep PSL in the same order, even if it creates a `-' on the target. */
    PslMapping* mapPsl(struct psl* inPsl) const;

};

#endif
