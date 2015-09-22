/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef targetAnnotations_hh
#define targetAnnotations_hh

#include "gxf.hh"
#include <map>
#include <stdexcept>

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
                    const string& seqIdForParCheck) const {
        string baseId = getBaseId(id);
        IdFeatureMapConstIter it = fIdFeatureMap.find(baseId);
        if (it == fIdFeatureMap.end()) {
            return NULL;
        } else if (it->second.size() == 2) {
            if (it->second[0]->fSeqid == seqIdForParCheck) {
                return it->second[0];
            } else if (it->second[1]->fSeqid == seqIdForParCheck) {
                return it->second[1];
            } else {
                throw logic_error("PAR target feature hack confused: " + baseId);
            }
        } else {
            return it->second[0];
        }
    }
};

#endif
