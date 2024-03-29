/*
 * Mapping of GxF features using transmap
 */
#ifndef featureTransMap_hh
#define featureTransMap_hh
#include "transMap.hh"
#include "featureTree.hh"
struct psl;
class PslMapping;

/* conversion of a list of features to a PSL */
class FeaturesToPsl {
    private:
    static bool checkFeatureOrder(const FeatureNodeVector& features);
    static int sumFeatureSizes(const FeatureNodeVector& features);
    static void makePslBlocks(struct psl* psl,
                              const FeatureNodeVector& features);
    static struct psl* makeFeaturesPsl(const string& qName,
                                       int qSize, int tStart, int tEnd, int tSize,
                                       const FeatureNodeVector& features);
    public:
    /* create a psl from a list of features. Assumes features are sorter in
     * ascending order.  */
    static struct psl* toPsl(const string& qName,
                             int tSize,
                             const FeatureNodeVector& features);
    
};

/* Mapping of GxF features using transmap. */
class FeatureTransMap {
    private:
    TransMapVector fTransMaps;  // object(s) to performance mappings. NOT OWNED!
    
    void mapPslVector(const PslVector& srcPsls,
                      int iTransMap,
                      PslVector& mappedPsls) const ;
    PslVector recursiveMapPsl(struct psl* srcPsl,
                              int iTransMap) const;

    static TransMapVector mkVector(const TransMap* transMap) {
        TransMapVector transMaps;
        transMaps.push_back(transMap);
        return transMaps;
    }
    public:
    /* Constructor with single transmap.  Doesn't own TransMap objects */
    FeatureTransMap(const TransMap* transMap):
        fTransMaps(mkVector(transMap)) {
    }

    /* Constructor with mulitple transmap object, in order of projection
     * Doesn't own TransMap objects */
    FeatureTransMap(const TransMapVector& transMaps):
        fTransMaps(transMaps) {
    }

    /* destructor */
    ~FeatureTransMap() {
    }

    /* Convert features to PSL and map via mapping alignment(s).  PSL will be
     * kept in the input order of the features, even if this creates a
     * negative target strand. */
    PslMapping* mapFeatures(const string& qName,
                            const FeatureNodeVector& features) const;

    /* map a single feature */
    PslMapping* mapFeature(const string& qName,
                           const FeatureNode* feature) const;

};



#endif
