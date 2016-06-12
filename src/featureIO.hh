/*
 * Paerser
 */
#ifndef featureIO_hh
#define featureIO_hh
#include <assert.h>
#include "feature.hh"
class GxfParser;

/**
 * Parser to group genes records together in a tree.
 */
class FeatureParser {
    private:
    GxfParser* fGxfParser;
    Feature* fNextGene;  // hold pending feature
    
    Feature* nextFeature();
    Feature* nextGeneFeature();
    Feature* findGff3Parent(Feature* geneLeaf,
                            const Feature* feature);
    const string& getGtfParentType(const string& featureType);
    Feature* findGtfParent(Feature* geneLeaf,
                           const Feature* feature);
    Feature* findParent(Feature* geneLeaf,
                        const Feature* feature);
    void addGeneFeature(Feature* gene,
                        Feature*& geneLeaf,
                        Feature* feature);
    void loadGeneChildren(Feature* gene);
    static void removeTransAttrsOnGenes(Feature* gene);
    static GxfFeature* featureFactory(const string& seqid, const string& source, const string& type,
                                      int start, int end, const string& score, const string& strand,
                                      const string& phase, const AttrVals& attrs);
    
    public:
    /* Constructor */
    FeatureParser(const string& gxfFile);

    /* Destructor */
    ~FeatureParser();
        
    /* load next gene */
    Feature* nextGene();
};

#endif
