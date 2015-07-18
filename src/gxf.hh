/*
 * Extremely naive and specialized GFF3 and GTF parsers.
 *
 * The goal preserve the exact structure of the GFF3/GTF files, including
 * comments, and only update coordinates (and occasionally split lines).
 * Neither parser in the kent tree design for this.
 *
 * These parser assume the ordering of the GENCODE GFF3/GTF files.
 */

#ifndef gxf_hh
#define gxf_hh
#include "typeOps.hh"
#include <queue>
using namespace std;

/* it seems so stupid to need to keep writing one-off GFF/GTF parsers */

class Feature;
class FIOStream;

typedef enum {
    GFF3_FORMAT,
    GTF_FORMAT,
} GxfFormat;

/* vector of feature objects */
typedef vector<Feature*>  FeatureVector;

/*
 * GxF base record type.  Use instanceOf to determine actually type
 */
class GxfRecord {
    public:
};

/*
 * non-feature line.
 */
class GxfLine: public string, public GxfRecord {
    public:
    GxfLine(const string& line):
        string(line) {
    }
};

/*
 * A row parsed from a GTF/GFF file. Immutable object.
 */
class GxfFeature: public GxfRecord {
public:
    // attribute/value pairs
    class AttrVal {
        public:
        const string fAttr;
        const string fVal;
        const bool fQuoted;

        AttrVal(const string& attr, const string& val, bool quoted):
            fAttr(attr), fVal(val), fQuoted(quoted) {
        }
    };
    typedef vector<AttrVal> AttrVals;
    
    // columns parsed from file.
    const string fSeqid;
    const string fSource;
    const string fType;
    const int fStart;
    const int fEnd;
    const string fScore;
    const string fStrand;
    const string fPhase;
    const AttrVals fAttrs;

public:
    /* construct a new feature object */
    GxfFeature(const string& seqid, const string& source, const string& type,
               int start, int end, const string& score, const string& strand,
               const string& phase, const AttrVals& attrs):
        fSeqid(seqid), fSource(source), fType(type),
        fStart(start), fEnd(end),
        fScore(score), fStrand(strand),
        fPhase(phase), fAttrs(attrs) {
    }

    /* destructor */
    virtual ~GxfFeature() {
    }
};

/**
 * gff3 or gtf parser.
 */
class GxfParser {
    private:
    FIOStream* fIn;  // input stream
    GxfFormat fGxfFormat; // format of file
    queue<const GxfRecord*> fPending; // FIFO of pushed records to be read before file

    StringVector splitFeatureLine(const string& line);
    const GxfFeature* createGff3Feature(const StringVector& columns);
    const GxfFeature* createGtfFeature(const StringVector& columns);
    const GxfFeature* createGxfFeature(const StringVector& columns);
    const GxfRecord* read();
    
    public:
    /* constructor that opens file, which maybe compressed */
    GxfParser(const string& fileName,
              GxfFormat gxfFormat);

    /* destructor */
    ~GxfParser();

    /* Read the next record, either queued by push() or from the file , use
     * instanceOf to determine the type.  Return NULL on EOF.
     */
    const GxfRecord* next();

    /* Return a record to be read before the file. */
    void push(const GxfRecord* gxfRecord) {
        fPending.push(gxfRecord);
    }
};

/**
 * Tree container for GxfFeature objects
 */
class GxfFeatureTree {
    public:
    const GxfFeature* fNode;
    vector<const GxfFeature*> fChildren;

    GxfFeatureTree(const GxfFeature* node):
        fNode(node) {
    }

    ~GxfFeatureTree() {
        for (int i = 0; i < fChildren.size(); i++) {
            delete fChildren[i];
        }
    }
};

#endif
