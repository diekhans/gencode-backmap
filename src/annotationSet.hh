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
class GenomeSizeMap;
class GxfWriter;

/*
 * Locations in target genome of old transcripts, by base id
 */
class AnnotationSet {
    public:
    // Maps of gene or transcripts id or name feature object.  A list is kept
    // because gene names can have multiple associated gene ids.  including
    // chrom handles pre-V44 releases where the same ids are used for both PAR
    // regions
    typedef map<const string, FeatureNodeVector> FeatureMap;
    typedef FeatureMap::iterator FeatureMapIter;
    typedef FeatureMap::const_iterator FeatureMapConstIter;

    private:

    /* stored in range tree to link target location to tree. */
    struct LocationLink {
        struct LocationLink *next;
        FeatureNode* feature;
    };

    // Warn if an id is added to the set were there is already a different version
    bool fWarnIdDiffVersions;
    
    // map by base id of genes and transcripts (not exons).
    FeatureMap fIdFeatureChromMap;
   
    // map by names of genes; small non-coding RNAs are know to have non-unique names and are not included.
    FeatureMap fNameFeatureChromMap;

    // map of base genes and transcript ids. Used to detect multiple version
    // of a transcript getting added, which has happened with bugs in target
    // substitutions.
    FeatureMap fIdFeatureMap;
    
    // list of all gene features found
    FeatureNodeVector fGenes;

    // map of location to feature
    struct genomeRangeTree* fLocationMap;

    // mapped sequence ids that have been written
    StringSet fSeqRegionsWritten;

    // optional table of chromosome sequence sizes
    const GenomeSizeMap* fGenomeSizes;

    string mkFeatureIdKey(const string& typeId,
                          const string& chrom) const;
    void insertInFeatureMap(const string& key,
                            FeatureNode* feature,
                            FeatureMap& featureMap);
    void addFeature(FeatureNode* feature);
    void processRecord(GxfParser *gxfParser,
                       GxfRecord* gxfRecord);
    void addLocationMap(FeatureNode* feature);
    void buildLocationMap();
    void freeLocationMap();
    void dumpFeatureMap(const FeatureMap& featureMap,
                        const string& label,
                        ostream& fh) const;
    bool isOverlappingGene(const FeatureNode* gene,
                           const FeatureNode* overlappingFeature,
                           float minSimilarity,
                           bool manualOnlyTranscripts);
    FeatureNode* getFeatureByKey(const string& baseKey,
                                 const string& chrom,
                                 const FeatureMap& featureMap) const;
    void checkForMultiIdVersions(FeatureNode* feature) const;

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
    void outputMappedSeqRegionIfNeed(const FeatureNode* gene,
                                     GxfWriter& mappedGxfFh);
    void outputFeature(const FeatureNode* feature,
                       GxfWriter& gxfFh) const;

    public:
    /* constructor, load gene and transcript objects from a GxF */
    AnnotationSet(const string& gxfFile,
                  const GenomeSizeMap* genomeSizes=NULL,
                  bool warnIdDiffVersions=false);

    /* constructor, empty set */
    AnnotationSet(const GenomeSizeMap* genomeSizes=NULL,
                  bool warnIdDiffVersions=false):
        fWarnIdDiffVersions(warnIdDiffVersions),
        fLocationMap(NULL),
        fGenomeSizes(genomeSizes) {
    }

    /* destructor */
    ~AnnotationSet();

    /* add a gene the maps */
    void addGene(FeatureNode* gene);

    /* get list of features with the base id or empty if none */
    const FeatureNodeVector& getFeaturesById(const string& baseId) const {
        static const FeatureNodeVector emptyVector;
        FeatureMapConstIter it = fIdFeatureMap.find(baseId);
        if (it == fIdFeatureMap.end()) {
            return emptyVector;
        } else {
            return it->second;
        }
    }

    /* get a target gene or transcript node with same base id or NULL.
     * special handling for PARs. Getting node is used if you need whole tree. */
    FeatureNode* getFeatureByIdChrom(const string& id,
                                     const string& chrom) const {
        return getFeatureByKey(getBaseId(id), chrom, fIdFeatureChromMap);
    }

    /* get a target gene or transcript node with same name or NULL.  Chrom is used
     * for for PARs. Getting node is used if you need whole tree. */
    FeatureNode* getFeatureByNameChrom(const string& name,
                                       const string& chrom) const {
        return getFeatureByKey(name, chrom, fNameFeatureChromMap);
    }

    /* find overlapping features */
    FeatureNodeVector findOverlappingFeatures(const string& seqid,
                                              int start,
                                          int end);
    
    /* find overlapping genes with minimum similarity at the transcript level */
    FeatureNodeVector findOverlappingGenes(const FeatureNode* gene,
                                           float minSimilarity,
                                           bool manualOnlyTranscripts);

    /* get list of all gene features */
    const FeatureNodeVector& getGenes() const {
        return fGenes;
    }

    /* sort genes */
    void sortGenes() {
        fGenes.sortChrom();
    }

    /* sort into standard gencode order */
    void sortGencode();
    
    /* print genes for debugging */
    void dump(ostream& fh) const;

    /* print id maps for debugging */
    void dumpIdMaps(ostream& fh) const;

    /* output genes */
    void write(GxfWriter& gxfFh);
};

#endif
