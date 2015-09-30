/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef targetAnnotations_hh
#define targetAnnotations_hh

#include "gxf.hh"
#include <map>
#include <stdexcept>
#include "featureTree.hh"
struct genomeRangeTree;

/* stored in range tree to link target location to
 * tree. */
struct TargetLocationLink {
    struct TargetLocationLink *next;
    FeatureNode* featureNode;
};
 
/*
 * Locations in target genome of old transcripts, by base id
 */
class TargetAnnotations {
    private:
       
    // map of gene or transcripts to features. Keep up to two for PAR
    typedef map<const string, FeatureNodeVector> IdFeatureMap;
    typedef IdFeatureMap::iterator IdFeatureMapIter;
    typedef IdFeatureMap::const_iterator IdFeatureMapConstIter;

    // map by base id of genes and transcripts
    IdFeatureMap fIdFeatureMap;

    // list of all gene features found
    FeatureNodeVector fGenes;
    
    // map of location used to find mapping to other locis
    struct genomeRangeTree* fLocationMap;
    
    void loadFeature(FeatureNode* featureNode);
    void loadGene(FeatureNode* geneTree);
    void processRecord(GxfParser *gxfParser,
                       GxfRecord* gxfRecord);
    void freeLocationMap();

    public:
    /* constructor, load gene and transcript objects from a GxF */
    TargetAnnotations(const string& gxfFile);

    /* destructor */
    ~TargetAnnotations();

    /* get a target gene or transcript with same base or NULL.
     * special handling for PARs/ */
    GxfFeature* getFeature(const string& id,
                           const string& seqIdForParCheck) const;

    /* find overlapping features */
    FeatureNodeVector findOverlappingFeatures(const string& seqid,
                                             int start,
                                             int end);
};

#endif
