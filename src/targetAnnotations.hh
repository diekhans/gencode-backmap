/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef targetAnnotations_hh
#define targetAnnotations_hh

#include "gxf.hh"
#include <map>
#include <stdexcept>
struct genomeRangeTree;

/* stored in range tree to link target location to
 * tree. */
struct TargetLocationLink {
    struct TargetLocationLink *next;
    GxfFeature* feature;
};

/*
 * Locations in target genome of old transcripts, by base id
 */
class TargetAnnotations {
    private:
    // Keep up to two for PAR
    typedef map<const string, GxfFeatureVector> IdFeatureMap;
    typedef IdFeatureMap::iterator IdFeatureMapIter;
    typedef IdFeatureMap::const_iterator IdFeatureMapConstIter;

    // map by base id of genes and transcripts
    IdFeatureMap fIdFeatureMap;

    // map of location used to find mapping to other locis
    struct genomeRangeTree* fLocationMap;
    
    void loadFeature(GxfFeature* gxfFeature);
    void processRecord(GxfRecord* gxfRecord);

    public:
    /* constructor, load gene and transcript objects from a GxF */
    TargetAnnotations(const string& gxfFile);

    /* destructor */
    ~TargetAnnotations();

    /* get a target gene or transcript with same base or NULL.
     * special handling for PARs/ */
    GxfFeature* get(const string& id,
                    const string& seqIdForParCheck) const;

    /* find overlapping features */
    GxfFeatureVector findOverlappingFeatures(const string& seqid,
                                             int start,
                                             int end);
};

#endif
