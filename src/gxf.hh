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
#include <algorithm>
#include <cassert>
using namespace std;

/* it seems so stupid to need to keep writing one-off GFF/GTF parsers */

class GxfFeature;
class FIOStream;

typedef enum {
    GXF_UNKNOWN_FORMAT,
    GFF3_FORMAT,
    GTF_FORMAT,
    DEV_NULL_FORMAT,
} GxfFormat;

/*
 * Method used to hack PAR ids to be unique when required by format.
 */
typedef enum {
    PAR_ID_HACK_OLD,  // ENSTR method
    PAR_ID_HACK_NEW   // _PAR_Y method
} ParIdHackMethod;

/* Get format from file name, or error */
GxfFormat gxfFormatFromFileName(const string& fileName);


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
typedef vector<GxfRecord*> GxfRecordVector;

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

/* attribute/value pair.  Maybe multi-valued */
class AttrVal {
    private:
    const string fName;
    StringVector fVals;

    static void checkName(const string& name) {
        if (stringEmpty(name)) {
            throw invalid_argument("empty attribute name");
        }
    }
    static void checkVal(const string& val) {
        if (stringEmpty(val)) {
            throw invalid_argument("empty attribute value");
        }
    }

    public:
    AttrVal(const string& name, const string& val):
        fName(name) {
        checkName(name);
        checkVal(val);
        fVals.push_back(val);
    }

    AttrVal(const string& name, const StringVector& vals):
        fName(name), fVals(vals) {
        checkName(name);
        for (int i = 0; i < vals.size(); i++) {
            checkVal(vals[i]);
        }
    }

    /* add a value */
    void addVal(const string& val) {
        checkVal(val);
        fVals.push_back(val);
    }
    
    /* copy constructor */
    AttrVal(const AttrVal& src):
        fName(src.fName), fVals(src.fVals) {
    }

    const string& getName() const {
        return fName;
    }
    const string& getVal(int iVal=0) const {
        return fVals[iVal];
    }
    const StringVector& getVals() const {
        return fVals;
    }
    int size() const {
        return fVals.size();
    }
};

/* vector of attribute/value pointers */
typedef vector<AttrVal*> AttrValVector;

/* list of attributes,  Multi-valued attributes (tag) are stored as multiple 
 * entries. */
class AttrVals: public AttrValVector {
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

    /* does the attribute exist */
    bool exists(const string& name) const {
        return findIdx(name) >= 0;
    }
    
    /* find the index of the first attribute with name or -1 if not found */
    int findIdx(const string& name) const {
        for (int i = 0; i < size(); i++) {
            if ((*this)[i]->getName() == name) {
                return i;
            }
        }
        return -1;
    }

    
    /* get a attribute, NULL if it doesn't exist */
    const AttrVal* find(const string& name) const {
        int i = findIdx(name);
        if (i < 0) {
            return NULL;
        } else {
            return (*this)[i];
        }
    }
    
    
    /* get a attribute, error it doesn't exist */
    const AttrVal* get(const string& name) const {
        const AttrVal* attrVal = find(name);
        if (attrVal == NULL) {
            throw invalid_argument("Attribute not found: " + name);
        }
        return attrVal;
    }

    /* add an attribute */
    void add(const AttrVal& attrVal) {
        push_back(new AttrVal(attrVal));
    }

    /* insert an attribute at the front */
    void push(const AttrVal& attrVal) {
        insert(begin(), new AttrVal(attrVal));
    }

    /* add or replace an attribute */
    void update(const AttrVal& attrVal) {
        int idx = findIdx(attrVal.getName());
        if (idx < 0) {
            add(attrVal);
        } else {
            delete (*this)[idx];
            (*this)[idx] = new AttrVal(attrVal);
        }
    }

    /* remove an attribute by name, if it exists */
    void remove(const string& attrName) {
        int idx = findIdx(attrName);
        if (idx >= 0) {
            delete (*this)[idx];
            erase(begin()+idx);
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

    // standard attribute names
    static const string ID_ATTR;
    static const string PARENT_ATTR;
    static const string GENE_ID_ATTR;
    static const string GENE_NAME_ATTR;
    static const string GENE_TYPE_ATTR;
    static const string GENE_STATUS_ATTR;
    static const string GENE_HAVANA_ATTR;
    static const string TRANSCRIPT_ID_ATTR;
    static const string TRANSCRIPT_NAME_ATTR;
    static const string TRANSCRIPT_TYPE_ATTR;
    static const string TRANSCRIPT_STATUS_ATTR;
    static const string TRANSCRIPT_HAVANA_ATTR;
    static const string EXON_ID_ATTR;
    static const string EXON_NUMBER_ATTR;
    static const string TAG_ATTR;
    
    /* source names */
    static const string SOURCE_HAVANA;
    static const string SOURCE_ENSEMBL;

    private:
    // columns parsed from file.
    const string fSeqid;
    const string fSource;
    const string fType;
    const int fStart;
    const int fEnd;
    const string fScore;
    const string fStrand;
    const string fPhase;
    AttrVals fAttrs;     // attribute maybe modified

    public:
    /* construct a new feature object */
    GxfFeature(const string& seqid, const string& source, const string& type,
               int start, int end, const string& score, const string& strand,
               const string& phase, const AttrVals& attrs):
        fSeqid(seqid), fSource(source), fType(type),
        fStart(start), fEnd(end),
        fScore(score), fStrand(strand),
        fPhase(phase), fAttrs(attrs) {
        assert(strand.size() == 1);
        assert(phase.size() == 1);
    }

    /* clone the feature */
    GxfFeature* clone() const {
        return new GxfFeature(fSeqid, fSource, fType, fStart, fEnd, fScore, fStrand, fPhase, fAttrs);
    }
    
    /* destructor */
    virtual ~GxfFeature() {
    }

    /* convert all columns, except attributes, to a string */
    string baseColumnsAsString() const;
    
    /* accessors */
    const string& getSeqid() const {
        return fSeqid;
    }
    const string& getSource() const {
        return fSource;
    }
    const string& getType() const {
        return fType;
    }
    int getStart() const {
        return fStart;
    }
    int getEnd() const {
        return fEnd;
    }
    /* get the length of the feature */
    int length() const {
        return (fEnd - fStart)+1;
    }
    const string& getScore() const {
        return fScore;
    }
    const string& getStrand() const {
        return fStrand;
    }
    const string& getPhase() const {
        return fPhase;
    }

    /* get all attribute */
    const AttrVals& getAttrs() const {
        return fAttrs;
    }

    /* get all attribute */
    AttrVals& getAttrs() {
        return fAttrs;
    }
    
    /* does the attribute exist */
    bool hasAttr(const string& name) const {
        return fAttrs.exists(name);
    }

    /* get a attribute, NULL if it doesn't exist */
    const AttrVal* findAttr(const string& name) const {
        return fAttrs.find(name);
    }

    /* get a attribute, error it doesn't exist */
    const AttrVal* getAttr(const string& name) const {
        return fAttrs.get(name);
    }

    /* get a attribute value, error it doesn't exist */
    const string& getAttrValue(const string& name) const {
        return getAttr(name)->getVal();
    }

    /* get a attribute value, default it doesn't exist */
    const string& getAttrValue(const string& name, 
                               const string& defaultVal) const {
        const AttrVal* attrVal = findAttr(name);
        if (attrVal == NULL) {
            return defaultVal;
        } else {
            return attrVal->getVal();
        }
    }

    /* get the id based on feature type, or empty string if it doesn't have an
     * id */
    const string& getTypeId() const;
    
    /* get the id based on feature type, or empty string if it doesn't have an
     * id */
    const string& getHavanaTypeId() const;
    
    /* get the name based on feature type, or empty string if it doesn't have an
     * id */
    const string& getTypeName() const;
    
    /* get the biotype based on feature type, or empty string if it doesn't have an
     * id */
    const string& getTypeBiotype() const;
    
    /* does this feature overlap another */
    bool overlaps(const GxfFeature* other) const {
        if ((fSeqid != other->fSeqid) || (fStrand != other->fStrand)) {
            return false;
        } else if ((fStart > other->fEnd) || (fEnd < other->fStart)) {
            return false;
        } else {
            return true;
        }
    }

    /* return feature as a string */
    virtual string toString() const;
};

/* vector of feature objects, doesn't own features */
class GxfFeatureVector: public vector<GxfFeature*> {
    public:
    /* free all features in the vector */
    void free() {
        for (int i = 0; i < size(); i++) {
            delete (*this)[i];
        }
        clear();
    }

    /* do we contain a particular feature object? */
    bool contains(const GxfFeature* feature) const {
        return std::find(begin(), end(), feature) != end();
    }
    
    /* sort the vector */
    void sort() {
        std::sort(begin(), end(),
                  [](const GxfFeature* a, const GxfFeature* b) -> bool {
                      if (a->getStrand() == "+") {
                          return a->getStart() > b->getStart();
                      } else {
                          return a->getStart() < b->getStart();
                      }   
                  });
    }
};

/**
 * gff3 or gtf parser.
 */
class GxfParser {
    private:
    FIOStream* fIn;  // input stream
    queue<GxfRecord*> fPending; // FIFO of pushed records to be read before file

    StringVector splitFeatureLine(const string& line) const;
    GxfRecord* read();

    protected:
    /* parse a feature */
    virtual GxfFeature* parseFeature(const StringVector& columns) = 0;
    
    /* constructor that opens file */
    GxfParser(const string& fileName);

    public:
    /* destructor */
    virtual ~GxfParser();

    /* get the format being parser */
    virtual GxfFormat getFormat() const = 0;

    /* Factory to create a parser. file maybe compressed.  If gxfFormat is
     * unknown, guess from filename*/
    static GxfParser *factory(const string& fileName,
                              GxfFormat gxfFormat=GXF_UNKNOWN_FORMAT);

    /* Read the next record, either queued by push() or from the file , use
     * instanceOf to determine the type.  Return NULL on EOF.
     */
    GxfRecord* next();

    /* Return a record to be read before the file. */
    void push(GxfRecord* gxfRecord) {
        fPending.push(gxfRecord);
    }
};

/**
 * gff3 or gtf writer.
 */
class GxfWriter {
    private:
    FIOStream* fOut;  // output stream

    protected:
    /* format a feature line */
    virtual string formatFeature(const GxfFeature* feature) = 0;
    
    public:
    /* constructor that opens file */
    GxfWriter(const string& fileName);

    /* destructor */
    virtual ~GxfWriter();

    /* get the format being written */
    virtual GxfFormat getFormat() const = 0;

    static GxfWriter *factory(const string& fileName,
                              ParIdHackMethod parIdHackMethod,
                              GxfFormat gxfFormat=GXF_UNKNOWN_FORMAT);

    /* copy a file to output, normally used for a header */
    void copyFile(const string& inFile);

    /* write one GxF record. */
    void write(const GxfRecord* gxfRecord);

    /* write one GxF line. */
    void write(const string& line);
    
};
#endif
