
#include "gxf.hh"
#include "FIOStream.hh"
#include <typeinfo>
#include <string>
#include <vector>
#include "typeOps.hh"
#include <stdexcept>

/* is a value quotes */
static bool isQuoted(const string& s) {
    return ((s.size() > 1) and (s[0] == '"') and (s[s.size()-1] == '"'));
}

/* strip optional quotes */
static string stripQuotes(const string& s) {
    if (isQuoted(s)) {
        return s.substr(1, s.size()-2);
    } else {
        return s;
    }
}

const string GxfFeature::GENE = "gene";
const string GxfFeature::TRANSCRIPT = "transcript";
const string GxfFeature::EXON = "exon";
const string GxfFeature::CDS = "CDS";
const string GxfFeature::START_CODON = "start_codon";
const string GxfFeature::UTR = "UTR";
const string GxfFeature::STOP_CODON = "stop_codon";
const string GxfFeature::STOP_CODON_REDEFINED_AS_SELENOCYSTEINE = "stop_codon_redefined_as_selenocysteine";


// standard attribute names
const string GxfFeature::ID_ATTR = "ID";
const string GxfFeature::PARENT_ATTR = "Parent";
const string GxfFeature::GENE_ID_ATTR = "gene_id";
const string GxfFeature::GENE_NAME_ATTR = "gene_name";
const string GxfFeature::GENE_TYPE_ID_ATTR = "gene_type";
const string GxfFeature::TRANSCRIPT_ID_ATTR = "transcript_id";
const string GxfFeature::TRANSCRIPT_NAME_ATTR = "transcript_name";
const string GxfFeature::TRANSCRIPT_TYPE_ID_ATTR = "transcript_type";
const string GxfFeature::EXON_ID_ATTR = "exon_id";
    

/* return base columns (excluding attributes) as a string */
string GxfFeature::baseColumnsAsString() const {
    return fSeqid + "\t" + fSource + "\t" + fType + "\t" + to_string(fStart) + "\t"
        + to_string(fEnd) + "\t" + fScore + "\t" + fStrand + "\t" + fPhase + "\t";
}

/* get the id based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getTypeId() const {
    static const string emptyStr;
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_ID_ATTR, emptyStr);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_ID_ATTR, emptyStr);
    } else if (fType == GxfFeature::EXON) {
        return getAttrValue(GxfFeature::EXON_ID_ATTR, emptyStr);
    } else {
        return emptyStr;
    }
}

/* get the name based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getTypeName() const {
    static const string emptyStr;
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_NAME_ATTR, emptyStr);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_NAME_ATTR, emptyStr);
    } else {
        return emptyStr;
    }
}

/* get the biotype based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getTypeBiotype() const {
    static const string emptyStr;
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_TYPE_ID_ATTR, emptyStr);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_TYPE_ID_ATTR, emptyStr);
    } else {
        return emptyStr;
    }
}


/*
 * GFF3 Feature
 */
class Gff3Feature: public GxfFeature {
    public:
    /* get the format */
    virtual GxfFormat getFormat() const {
        return GFF3_FORMAT;
    }
    
    /* parse ID=ENSG00000223972.5 */
    static AttrVal* parseAttr(const string& attrStr) {
        size_t i = attrStr.find('=');
        if (i == string::npos) {
            throw invalid_argument("Invalid GFF3 attribute \"" + attrStr + "\"");
        }
        string name = attrStr.substr(0,i);
        string value = attrStr.substr(i+1);
        return new AttrVal(name, stripQuotes(value), isQuoted(value));
    }
    
    /* parse: ID=ENSG00000223972.5;gene_id=ENSG00000223972.5 */
    static AttrVals parseAttrs(const string& attrsStr) {
        AttrVals attrVals;
        StringVector parts = stringSplit(attrsStr,';');
        // `;' is a separator
        for (size_t i = 0; i < parts.size(); i++) {
            attrVals.push_back(parseAttr(stringTrim(parts[i])));
        }
        return attrVals;
    }
    /* format an attribute */
    string formatAttr(const AttrVal* attrVal) const {
        // n.b. this is not general, doesn't handle embedded quotes
        string strAttr = attrVal->fName + "=";
        if (attrVal->fQuoted) {
            strAttr += "\"" + attrVal->fVal + "\"";
        } else {
            strAttr += attrVal->fVal;
        }
        return strAttr;
    }

    /* format attributes */
    string formatAttrs() const {
        string strAttrs;
        for (size_t i = 0; i < fAttrs.size(); i++) {
            if (i > 0) {
                strAttrs += ";";
            }
            strAttrs += formatAttr(fAttrs[i]);
        }
        return strAttrs;
    }
    
    /* constructor. */
    Gff3Feature(const string& seqid, const string& source, const string& type, int
                start, int end, const string& score, const string& strand, const
                string& phase, const AttrVals& attrs):
        GxfFeature(seqid, source, type, start, end, score, strand, phase, attrs) {
    }

    /* clone the feature */
    virtual GxfFeature* clone() const {
        return new Gff3Feature(fSeqid, fSource, fType, fStart, fEnd, fScore, fStrand, fPhase, fAttrs);
    }
    
    /* return record as a string */
    virtual string toString() const {
        return baseColumnsAsString() + formatAttrs();
    }
};


/*
 * GTF Feature
 */
class GtfFeature: public GxfFeature {
    public:
    /* get the format */
    virtual GxfFormat getFormat() const {
        return GTF_FORMAT;
    }
    
    /* parse ID=ENSG00000223972.5 */
    static AttrVal* parseAttr(const string& attrStr) {
        size_t i = attrStr.find(' ');
        if (i == string::npos) {
            throw invalid_argument("Invalid GTF attribute \"" + attrStr + "\"");
        }
        string name = attrStr.substr(0,i);
        string value = attrStr.substr(i+1);
        return new AttrVal(name, stripQuotes(value), isQuoted(value));
    }

    /* parse: gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene";  */
    static AttrVals parseAttrs(const string& attrsStr) {
        AttrVals attrVals;
        StringVector parts = stringSplit(attrsStr,';');
        // last will be empty, since `;' is a terminator
        for (size_t i = 0; i < parts.size()-1; i++) {
            attrVals.push_back(parseAttr(stringTrim(parts[i])));
        }
        return attrVals;
    }

    /* format an attribute */
    string formatAttr(const AttrVal* attrVal) const {
        // n.b. this is not general, doesn't handle embedded quotes
        string strAttr = attrVal->fName + " ";
        if (attrVal->fQuoted) {
            strAttr += "\"" + attrVal->fVal + "\"";
        } else {
            strAttr += attrVal->fVal;
        }
        return strAttr;
    }

    /* format attribute */
    string formatAttrs() const {
        string strAttrs;
        for (size_t i = 0; i < fAttrs.size(); i++) {
            if (i > 0) {
                strAttrs += " ";  // same formatting as GENCODE
            }
            strAttrs += formatAttr(fAttrs[i]) + ";";
        }
        return strAttrs;
    }
    public:
    /* constructor */
    GtfFeature(const string& seqid, const string& source, const string& type, int
               start, int end, const string& score, const string& strand, const
               string& phase, const AttrVals& attrs):
        GxfFeature(seqid, source, type, start, end, score, strand, phase, attrs) {
    }

    /* clone the feature */
    virtual GxfFeature* clone() const {
        return new GtfFeature(fSeqid, fSource, fType, fStart, fEnd, fScore, fStrand, fPhase, fAttrs);
    }
    
    /* return record as a string */
    virtual string toString() const {
        return baseColumnsAsString() + formatAttrs();
    }
};


/* create the appropriate feature type */
GxfFeature* gxfFeatureFactory(GxfFormat gxfFormat,
                              const StringVector& columns) {
    if (columns.size() != 9) {
        throw invalid_argument("expected 9 columns in GxF file");
    }
    if (gxfFormat == GFF3_FORMAT) {
        return new Gff3Feature(columns[0], columns[1], columns[2],
                               stringToInt(columns[3]), stringToInt(columns[4]),
                               columns[5], columns[6], columns[7], Gff3Feature::parseAttrs(columns[8]));
    } else {
        return new GtfFeature(columns[0], columns[1], columns[2],
                              stringToInt(columns[3]), stringToInt(columns[4]),
                              columns[5], columns[6], columns[7], GtfFeature::parseAttrs(columns[8]));
    }
}

/* create the appropriate feature type */
GxfFeature* gxfFeatureFactory(GxfFormat gxfFormat,
                              const string& seqid, const string& source, const string& type,
                              int start, int end, const string& score, const string& strand, const string& phase,
                              const AttrVals& attrs) {
    if (gxfFormat == GFF3_FORMAT) {
        return new Gff3Feature(seqid, source, type, start, end, score, strand, phase, attrs);
    } else {
        return new GtfFeature(seqid, source, type, start, end, score, strand, phase, attrs);
    }
}

/* split a feature line of GFF3 or GTF */
StringVector GxfParser::splitFeatureLine(const string& line) const {
    StringVector columns = stringSplit(line, '\t');
    if (columns.size() != 9) {
        invalid_argument("invalid row, expected 9 columns: " + line);
    }
    return columns;
}

/* constructor that opens file, which maybe compressed. If gxfFormat is
 * unknown, guess from filename */
GxfParser::GxfParser(const string& fileName,
                     GxfFormat gxfFormat):
    fIn(new FIOStream(fileName)),
    fGxfFormat((gxfFormat==GXF_UNKNOWN_FORMAT) ? formatFromFileName(fileName) : gxfFormat) {
}

/* destructor */
GxfParser::~GxfParser() {
    delete fIn;
}

/* Get format from file name, or error */
GxfFormat GxfParser::formatFromFileName(const string& fileName) {
    if (stringEndsWith(fileName, ".gff3") or stringEndsWith(fileName, ".gff3.gz")) {
        return GFF3_FORMAT;
    } else if (stringEndsWith(fileName, ".gtf") or stringEndsWith(fileName, ".gtf.gz")) {
        return GTF_FORMAT;
    } else {
        errAbort(toCharStr("Error: expected input annotation with an extension of .gff3, .gff3.gz, .gtf, or .gtf.gz: " + fileName));
        return GXF_UNKNOWN_FORMAT;
    }
}

/* Read the next record */
GxfRecord* GxfParser::read() {
    string line;
    if (not fIn->readLine(line)) {
        return NULL;
    } else if ((line.size() > 0) and line[0] != '#') {
        return gxfFeatureFactory(fGxfFormat, splitFeatureLine(line));
    } else {
        return new GxfLine(line);
    }
}

/* Read the next record, either queued by push() or from the file , use
 * instanceOf to determine the type.  Return NULL on EOF.
 */
GxfRecord* GxfParser::next() {
    if (not fPending.empty()) {
        GxfRecord* rec = fPending.front();
        fPending.pop();
        return rec;
    } else {
        return read();
    }
}
