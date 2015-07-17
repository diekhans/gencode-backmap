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
#include <string>
#include <vector>
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
 * A row parsed from a GTF/GFF file. Immutable object.
 */
class Feature {
public:
    // columns parsed from file.
    const string fSeqid;
    const string fSource;
    const string fType;
    const int fStart;
    const int fEnd;
    const string fScore;
    const string fStrand;
    const string fPhase;
    const string fAttrs;

public:
    /* construct a new feature object */
    Feature(const string& seqid, const string& source, const string& type,
            int start, int end, const string& score, const string& strand,
            const string& phase, const string& attrs):
        fSeqid(seqid), fSource(source), fType(type),
        fStart(start), fEnd(end),
        fScore(score), fStrand(strand),
        fPhase(phase), fAttrs(attrs) {
    }

    /* destructor */
    virtual ~Feature() {
    }
};

/* Record returned by the parser, maybe or may not include a feature */
class GxfLine {
    public:
    const string fLine;
    const Feature* fFeature;

    public:
    GxfLine(const string& line,
            const Feature* feature = NULL):
        fLine(line),
        fFeature(feature) {
    }
    ~GxfLine() {
        delete fFeature;
    }
};

/**
 * gff3 or gtf parser.
 */
class GxfParser {
    private:
    FIOStream* fIn;  // input stream
    GxfFormat fGxfFormat; // format of file

    vector<const string> splitFeatureLine(const string& line);
    const Feature* createGff3Feature(const vector<const string>& columns);
    const Feature* createGtfFeature(const vector<const string>& columns);
    const Feature* createGxfFeature(const vector<const string>& columns);
    
    public:
    /* constructor that opens file, which maybe compressed */
    GxfParser(const string& fileName,
              GxfFormat gxfFormat);

    /* destructor */
    ~GxfParser();

    /* Read the next line, parse into a feature if it is one, otherwise just
     * return line. Return NULL on EOF */
    GxfLine* next();
};

#endif
