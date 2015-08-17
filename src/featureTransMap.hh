/*
 * Mapping of GxF features using transmap
 */
#ifndef featureTransMap_hh
#define featureTransMap_hh
#include "transMap.hh"
#include "gxf.hh"
struct psl;
class PslMapping;

/* Mapping of GxF features using transmap  */
class FeatureTransMap {
    private:
    TransMapVector fTransMaps;  // object(s) to performance mappings
    
    int sumFeatureSizes(const GxfFeatureVector& features) const;
    void makePslBlocks(struct psl* psl,
                       const GxfFeatureVector& features) const;
    struct psl* makeFeaturesPsl(const string& qName,
                                int qSize, int tStart, int tEnd, int tSize,
                                const GxfFeatureVector& features) const;
    struct psl* featuresToPsl(const string& qName,
                              const GxfFeatureVector& features) const;
    bool checkFeatureOrder(const GxfFeatureVector& features) const;

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
    /* Constructor with single transmap */
    FeatureTransMap(const TransMap* transMap):
        fTransMaps(mkVector(transMap)) {
    }
    /* Constructor with mulitple transmap object, in order of projection */
    FeatureTransMap(const TransMapVector transMaps):
        fTransMaps(transMaps) {
    }

    /* Convert features to PSL and map via mapping alignment(s).  PSL will be
     * kept in the input order of the features, even if this creates a
     * negative target strand. */
    PslMapping* mapFeatures(const string& qName,
                            const GxfFeatureVector& features) const;
};



#endif
