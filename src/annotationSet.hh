/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef annotationSet_hh
#define annotationSet_hh

#include "gxf.hh"
#include <map>
#include <stdexcept>
#include "featureTree.hh"
struct genomeRangeTree;

/* stored in range tree to link target location to
 * tree. */
struct LocationLink {
    struct LocationLink *next;
    FeatureNode* featureNode;
};
 
/*
 * Locations in target genome of old transcripts, by base id
 */
class AnnotationSet {
    private:
       
    // map of gene or transcripts to features. Keep up to two for PAR
    typedef map<const string, FeatureNodeVector> FeatureMap;
    typedef FeatureMap::iterator FeatureMapIter;
    typedef FeatureMap::const_iterator FeatureMapConstIter;

    // map by base id of genes and transcripts
    FeatureMap fIdFeatureMap;

    // map by names of genes and transcripts
    FeatureMap fNameFeatureMap;

    // list of all gene features found
    FeatureNodeVector fGenes;
    
    // map of location used to find mapping to other locis
    struct genomeRangeTree* fLocationMap;
    
    void addFeature(FeatureNode* featureNode);
    void addGene(FeatureNode* geneTree);
    void processRecord(GxfParser *gxfParser,
                       GxfRecord* gxfRecord);
    void addLocationMap(FeatureNode* featureNode);
    void buildLocationMap();
    void freeLocationMap();
    FeatureNode* getFeatureNodeByKey(const string& key,
                                     const FeatureMap& featureMap,
                                     const string& seqIdForParCheck) const;

    public:
    /* constructor, load gene and transcript objects from a GxF */
    AnnotationSet(const string& gxfFile);

    /* destructor */
    ~AnnotationSet();

    /* get a target gene or transcript node with same base id or NULL.
     * special handling for PARs. Getting node is used if you need whole tree. */
    FeatureNode* getFeatureNodeById(const string& id,
                                    const string& seqIdForParCheck) const;

    /* get a target gene or transcript node with same name or NULL.
     * special handling for PARs. Getting node is used if you need whole tree. */
    FeatureNode* getFeatureNodeByName(const string& name,
                                      const string& seqIdForParCheck) const;

    /* get a target gene or transcript with same base or NULL.
     * special handling for PARs. */
    GxfFeature* getFeatureById(const string& id,
                               const string& seqIdForParCheck) const;

    /* find overlapping features */
    FeatureNodeVector findOverlappingFeatures(const string& seqid,
                                             int start,
                                             int end);

    /* get list of all gene features */
    const FeatureNodeVector& getGenes() const {
        return fGenes;
    }
};

#endif
