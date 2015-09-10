/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef targetAnnotations_hh
#define targetAnnotations_hh

#include "gxf.hh"
#include <map>

/*
 * Locations in target genome of old transcripts, by base id
 */
class TargetAnnotations {
    private:
    typedef map<const string, GxfFeature*> IdFeatureMap;
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

    /* get a target gene or transcript with same base or NULL */
    GxfFeature* get(const string& id) const {
        string baseId = getBaseId(id);
        IdFeatureMapConstIter it = fIdFeatureMap.find(baseId);
        if (it == fIdFeatureMap.end()) {
            return NULL;
        } else {
            return it->second;
        }
    }
};

#endif
