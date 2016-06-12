/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef annotationSet_hh
#define annotationSet_hh

#include "gxfRecord.hh"
#include <map>
#include <stdexcept>
#include "feature.hh"
struct genomeRangeTree;
class GenomeSizeMap;
class GxfWriter;

/*
 * Locations in target genome of old transcripts, by base id
 */
class AnnotationSet {
    private:
    /* stored in range tree to link target location to
     * tree. */
    struct LocationLink {
        struct LocationLink *next;
        Feature* feature;
    };
 
       
    // map of gene or transcripts to features. Keep up to two for PAR
    typedef map<const string, FeatureVector> FeatureMap;
    typedef FeatureMap::iterator FeatureMapIter;
    typedef FeatureMap::const_iterator FeatureMapConstIter;

    // map by base id of genes and transcripts (not exons).
    FeatureMap fIdFeatureMap;

    // map by base id of exons
    FeatureMap fIdExonMap;

    // map by names of genes and transcripts
    FeatureMap fNameFeatureMap;

    // list of all gene features found
    FeatureVector fGenes;

    // map of location to feature
    struct genomeRangeTree* fLocationMap;

    // mapped sequence ids that have been written
    StringSet fSeqRegionsWritten;

    // optional table of chromosome sequence sizes
    const GenomeSizeMap* fGenomeSizes;
    
    void addFeature(Feature* feature);
    void addLocationMap(Feature* feature);
    void buildLocationMap();
    void freeLocationMap();
    bool isOverlappingGene(const Feature* gene,
                           const Feature* overlappingFeature,
                           float minSimilarity,
                           bool manualOnlyTranscripts);

    Feature* getFeatureByKey(const string& key,
                                     const FeatureMap& featureMap,
                                     const string& seqIdForParCheck) const;

    /* check if a seqregion for seqid has been written, if so, return true,
     * otherwise record it and return false.  */
    bool checkRecordSeqRegionWritten(const string& seqid) {
        if (fSeqRegionsWritten.find(seqid) == fSeqRegionsWritten.end()) {
            fSeqRegionsWritten.insert(seqid);
            return false;
        } else {
            return true;
        }
    }
    void outputSeqRegion(const string& seqId,
                         int size,
                         GxfWriter& gxfFh);
    void outputMappedSeqRegionIfNeed(const Feature* gene,
                                     GxfWriter& mappedGxfFh);
    void outputFeature(const Feature* feature,
                       GxfWriter& gxfFh) const;

    public:
    /* constructor, load gene and transcript objects from a GxF */
    AnnotationSet(const string& gxfFile,
                  const GenomeSizeMap* genomeSizes=NULL);

    /* constructor, empty set */
    AnnotationSet(const GenomeSizeMap* genomeSizes=NULL):
        fLocationMap(NULL),
        fGenomeSizes(genomeSizes) {
    }

    /* destructor */
    ~AnnotationSet();

    /* add a gene the maps */
    void addGene(Feature* gene);

    /* get a gene or transcript with same base id or NULL.  special
     * handling for PARs. */
    Feature* getFeatureById(const string& id,
                            const string& seqIdForParCheck) const;

    /* get a gene or transcript with same name or NULL.  special handling
     * for PARs. */
    Feature* getFeatureByName(const string& name,
                              const string& seqIdForParCheck) const;

    /* get exon features by base id */
    FeatureVector getExonsById(const string& exonId) const;

    /* find overlapping features */
    FeatureVector findOverlappingFeatures(const string& seqid,
                                          int start,
                                          int end);
    
    /* find overlapping genes with minimum similarity at the transcript level */
    FeatureVector findOverlappingGenes(const Feature* gene,
                                       float minSimilarity,
                                       bool manualOnlyTranscripts);

    /* get list of all gene features */
    const FeatureVector& getGenes() const {
        return fGenes;
    }

    /* print for debugging */
    void dump(ostream& fh) const;

    /* output genes */
    void write(GxfWriter& gxfFh);
};

#endif
