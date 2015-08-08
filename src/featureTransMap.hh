/*
 * Mapping of GxF features using transmap
 */
#ifndef featureTransMap_hh
#define featureTransMap_hh
#include "transMap.hh"
#include "gxf.hh"
struct psl;

/* Mapping of GxF features using transmap  */
class FeatureTransMap {
    private:
    const TransMap* fTransMap;  // object to performance mappings
    
    int sumFeatureSizes(const GxfFeatureVector& features) const;
    void makePslBlocks(struct psl* psl,
                       const GxfFeatureVector& features) const;
    struct psl* featuresToPsl(const string& qName,
                              const GxfFeatureVector& features) const;
    bool checkFeatureOrder(const GxfFeatureVector& features) const;

    public:
    /* Constructor */
    FeatureTransMap(const TransMap* transMapper):
        fTransMap(transMapper) {
    }

    /* Map features to list of PSLs.  PSL will be kept in the input order
     * of the features, even if this creates a negative target strand. */
    PslMapping* mapFeatures(const string& qName,
                            const GxfFeatureVector& features) const;
};



#endif
