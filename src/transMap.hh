/*
 * transmap projection of annotations
 */
#ifndef transMap_hh
#define transMap_hh
#include "jkinclude.hh"
#include <string>
#include <map>
#include "pslOps.hh"
using namespace std;

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

    /* constructor */
    TransMap();

    public:
    /* factory from a chain file */
    static TransMap* factoryFromChainFile(const string& chainFile,
                                      bool swapMap);
    /* factory from a list of psls */
    static TransMap* factoryFromPsls(struct psl* psls,
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
    PslVector mapPsl(struct psl* inPsl) const;
};

/* Vector of transmap objects.  Doesn't own them. */
class  TransMapVector: public vector<const TransMap*> {
    public:
    /* free all objects in the vector */
    void free() {
        for (int i = 0; i < size(); i++) {
            delete (*this)[i];
        }
        clear();
    }
};
#endif
