/*
 * transmap projection of annotations
 */
#ifndef transMap_hh
#define transMap_hh
#include <string>
#include <map>
#include "pslOps.hh"
using namespace std;


class GenomeSizeMap: public map<const string, int> {
    public:
    /* add a size if we don't have it */
    void add(const string& name,
             int size) {
        if (find(name) == end()) {
            (*this)[name] = size;
        }
    }

    /* do we have a sequence query sequence */
    bool have(const string& name) const {
        return find(name) != end();
    }

    /* get the size of a sequence */
    int get(const string& name) const {
        return at(name);
    }
};

/*
 * transmap via alignment chains
 */
class TransMap {
    private:
    struct genomeRangeTree* fMapAlns;  // mapping alingments

    public:
    GenomeSizeMap fQuerySizes;   // query sequence sizes
    GenomeSizeMap fTargetSizes;  // target sequence sizes
   
    private:
    void mapAlnsAdd(struct psl *mapPsl);
    void mapPslPair(struct psl *inPsl,
                    struct psl *mapPsl,
                    PslVector& allMappedPsls) const;

    /* is a mapping alignment file a chain or psl? */
    static bool isChainMappingAlign(const string& fileName) {
        if (stringEndsWith(fileName, ".chain") or stringEndsWith(fileName, ".chain.gz")) {
            return true;
        } else if (stringEndsWith(fileName, ".psl") or stringEndsWith(fileName, ".psl.gz")) {
            return false;
        } else {
            errAbort(toCharStr("Error: expected mapping alignments file with an extension of .chain, .chain.gz, .psl, or .psl.gz: " + fileName));
            return false;
        }
    }

    /* constructor */
    TransMap();

    public:
    /* consumes PSLs */
    static TransMap* factoryFromPsls(struct psl** psls,
                                     bool swapMap);

    /* clones PSL */
    static TransMap* factoryFromPsl(struct psl* psl,
                                    bool swapMap);

    /* factory from a chain file */
    static TransMap* factoryFromChainFile(const string& chainFile,
                                      bool swapMap);
    /* factory from a  psl file */
    static TransMap* factoryFromPslFile(const string& pslFile,
                                        bool swapMap);
    
    /* factory from a chain or psl file */
    static TransMap* factoryFromFile(const string& fileName,
                                     bool swapMap) {
        if (isChainMappingAlign(fileName)) {
            return factoryFromChainFile(fileName, swapMap);
        } else {
            return factoryFromPslFile(fileName, swapMap);
        }
    }
    
    /* destructor */
    ~TransMap();

    /* do we have a mapping query sequence */
    bool haveQuerySeq(const string& qName) const {
        return fQuerySizes.have(qName);
    }
    
    /* do we have a mapping target sequence */
    bool haveTargetSeq(const string& tName) const {
        return fTargetSizes.have(tName);
    }
    
    /* get the size of a mapping query sequence */
    int getQuerySeqSize(const string& qName) const {
        return fQuerySizes.get(qName);
    }
    
    /* get the size of a mapping target sequence */
    int getTargetSeqSize(const string& tName) const {
        return fTargetSizes.get(tName);
    }
    
    /* Map a single input PSL and return a list of resulting mappings.  Keep
     * PSL in the same query order, even if it creates a `-' on the target. */
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
