#ifndef bedMap_hh
#define bedMap_hh
#include <string>
using namespace std;

struct genomeRangeTree;
struct bed;

/*
 * container to hold BED records to search for overlaps
 */
class BedMap {
    private:
    // map of location to bed
    struct genomeRangeTree* fLocationMap;

    public:
    /* constructor */
    BedMap(const string& bedFile);

    /* destructor */
    ~BedMap();

    /* check for overlap on a range */
    bool anyOverlap(const string& seqid,
                    int start,
                    int end) const;

};

#endif
