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
#include <stdexcept>
using namespace std;

/* it seems so stupid to need to keep writing one-off GFF/GTF parsers */

class GxfFeature;
class FIOStream;

typedef enum {
    GFF3_FORMAT,
    GTF_FORMAT,
} GxfFormat;

/*
 * GxF base record type.  Use instanceOf to determine actually type
 */
class GxfRecord {
    public:
    /* destructor */
    virtual ~GxfRecord() {
    }

    /* return record as a string */
    virtual string toString() const = 0;
};
typedef vector<const GxfRecord*> GxfRecordVector;

/*
 * non-feature line.
 */
class GxfLine: public string, public GxfRecord {
    public:
    GxfLine(const string& line):
        string(line) {
    }

    /* destructor */
    virtual ~GxfLine() {
    }

    /* return record as a string */
    virtual string toString() const {
        return *this;
    }
};

/* attribute/value pair */
class AttrVal {
    public:
    const string fName;
    const string fVal;
    const bool fQuoted;

    AttrVal(const string& name, const string& val, bool quoted=false):
        fName(name), fVal(val), fQuoted(quoted) {
        if (stringEmpty(fName)) {
            throw invalid_argument("empty attribute name");
        }
        if (stringEmpty(fVal)) {
            throw invalid_argument("empty attribute value");
        }
    }

    /* copy constructor */
    AttrVal(const AttrVal& src):
        fName(src.fName), fVal(src.fVal), fQuoted(src.fQuoted) {
    }
};

/* list of attributes, derived class of specific types are create for GFF3 and GTF */
class AttrVals: public vector<const AttrVal*> {
    // n.b.  this keeps pointers rather than values due to reallocation if vector changes
    public:
    /* empty constructor */
    AttrVals() {
    }

    /* copy constructor */
    AttrVals(const AttrVals& src) {
        for (size_t i = 0; i < src.size(); i++) {
            push_back(new AttrVal(*(src[i])));
        }
    }

    
    /* destructor */
    ~AttrVals() {
        for (size_t i = 0; i < size(); i++) {
            delete (*this)[i];
        }
    }

    /* find the index of an attribute or -1 if not found */
    int findIdx(const string& name) const {
        for (int i = 0; i < size(); i++) {
            if ((*this)[i]->fName == name) {
                return i;
            }
        }
        return -1;
    }

    
    /* get a attribute, NULL if it doesn't exist */
    const AttrVal* findAttr(const string& name) const {
        int i = findIdx(name);
        if (i < 0) {
            return NULL;
        } else {
            return (*this)[i];
        }
    }
    
    /* get a attribute, error it doesn't exist */
    const AttrVal* getAttr(const string& name) const {
        const AttrVal* attrVal = findAttr(name);
        if (attrVal == NULL) {
            throw invalid_argument("Attribute not found: " + name);
        }
        return attrVal;
    }

    /* add an attribute */
    void add(const AttrVal& attrVal) {
        push_back(new AttrVal(attrVal));
    }

    /* add or replace an attribute */
    void update(const AttrVal& attrVal) {
        int idx = findIdx(attrVal.fName);
        if (idx < 0) {
            add(attrVal);
        } else {
            delete (*this)[idx];
            (*this)[idx] = new AttrVal(attrVal);
        }
    }
};

/*
 * A row parsed from a GTF/GFF file. Immutable object.
 */
class GxfFeature: public GxfRecord {
public:
    // Standard feature names
    static const string GENE;
    static const string TRANSCRIPT;
    static const string EXON;
    static const string CDS;
    static const string START_CODON;
    static const string UTR;
    static const string STOP_CODON;
    static const string STOP_CODON_REDEFINED_AS_SELENOCYSTEINE;

    
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

    protected:
    string baseColumnsAsString() const;
    
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

    /* clone the feature */
    virtual GxfFeature* clone() const = 0;
    
    /* get the format */
    virtual GxfFormat getFormat() const = 0;
    
    /* destructor */
    virtual ~GxfFeature() {
    }

    /* get all attribute */
    const AttrVals& getAttrs() const {
        return fAttrs;
    }


    /* get a attribute, NULL if it doesn't exist */
    const AttrVal* findAttr(const string& name) const {
        return fAttrs.findAttr(name);
    }

    /* get a attribute, error it doesn't exist */
    const AttrVal* getAttr(const string& name) const {
        return fAttrs.getAttr(name);
    }

    /* get a attribute value, error it doesn't exist */
    const string& getAttrValue(const string& name) const {
        return getAttr(name)->fVal;
    }

    /* get the size of the feature */
    int size() const {
        return (fEnd - fStart)+1;
    }
};

/* create the appropriate feature type */
const GxfFeature* gxfFeatureFactory(GxfFormat gxfFormat,
                                    const StringVector& columns);

/* create the appropriate feature type */
const GxfFeature* gxfFeatureFactory(GxfFormat gxfFormat,
                                    const string& seqid, const string& source, const string& type,
                                    int start, int end, const string& score, const string& strand,
                                    const string& phase, const AttrVals& attrs);

/* vector of feature objects, doesn't own features */
class GxfFeatureVector: public vector<const GxfFeature*> {
    public:
    /* free all features in the vector */
    void free() {
        for (int i = 0; i < size(); i++) {
            delete (*this)[i];
        }
        clear();
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

    StringVector splitFeatureLine(const string& line) const;
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

    /* Accessors */
    GxfFormat getGxfFormat() const {
        return fGxfFormat;
    }
};

#endif
