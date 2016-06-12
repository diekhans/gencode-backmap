#include "gxfIO.hh"
#include "FIOStream.hh"
#include <typeinfo>
#include <string>
#include <vector>
#include "typeOps.hh"
#include <stdexcept>
#include "gxfRecord.hh"

/* is a value quotes */
static bool isQuoted(const string& s) {
    return ((s.size() > 1) and (s[0] == '"') and (s[s.size()-1] == '"'));
}


/* is a value a integrate or float */
static bool isNumeric(const string& s) {
    int dotCount = 0;
    for (int i = 0; i < s.size(); i++) {
        if (s[i] == '.') {
            dotCount++;
        } else if (!isdigit(s[i])) {
            return false;
        }
    }
    if (dotCount > 1) {
        return false;
    }
    return true;
}

/* strip optional quotes */
static string stripQuotes(const string& s) {
    if (isQuoted(s)) {
        return s.substr(1, s.size()-2);
    } else {
        return s;
    }
}

/* Get format from file name, or error */
GxfFormat gxfFormatFromFileName(const string& fileName) {
    if (stringEndsWith(fileName, ".gff3") or stringEndsWith(fileName, ".gff3.gz")
        or (fileName == "/dev/null")) {
        return GFF3_FORMAT;
    } else if (stringEndsWith(fileName, ".gtf") or stringEndsWith(fileName, ".gtf.gz")) {
        return GTF_FORMAT;
    } else {
        errAbort(toCharStr("Error: expected input annotation with an extension of .gff3, .gff3.gz, .gtf, or .gtf.gz: " + fileName));
        return GXF_UNKNOWN_FORMAT;
    }
}

/* Parse for GFF 3 */
class Gff3Parser: public GxfParser {
    private:
    /* is this a multi-valued attribute? */
    /* parse ID=ENSG00000223972.5 */
    static void parseAttr(const string& attrStr,
                          AttrVals& attrVals) {
        size_t i = attrStr.find('=');
        if (i == string::npos) {
            throw invalid_argument("Invalid GFF3 attribute \"" + attrStr + "\"");
        }
        string name = attrStr.substr(0,i);
        string value = attrStr.substr(i+1);
        StringVector values = stringSplit(stripQuotes(value), ',');
        AttrVal* attrVal = new AttrVal(name, values[0]);
        attrVals.push_back(attrVal);
        for (int i = 1; i < values.size(); i++) {
            attrVal->addVal(values[i]);
        }
    }

    /* parse: ID=ENSG00000223972.5;gene_id=ENSG00000223972.5 */
    static AttrVals parseAttrs(const string& attrsStr) {
        AttrVals attrVals;
        StringVector parts = stringSplit(attrsStr,';');
        // `;' is a separator
        for (size_t i = 0; i < parts.size(); i++) {
            parseAttr(stringTrim(parts[i]), attrVals);
        }
        return attrVals;
    }

    public:
    /* constructor */
    Gff3Parser(const string& fileName,
               GxfFeatureFactory gxfFeatureFactory):
        GxfParser(fileName, gxfFeatureFactory) {
    }
 
    /* get the format being parser */
    virtual GxfFormat getFormat() const {
        return GFF3_FORMAT;
    }

    /* parse a feature */
    virtual GxfFeature* parseFeature(const StringVector& columns) {
        return fGxfFeatureFactory(columns[0], columns[1], columns[2],
                                  stringToInt(columns[3]), stringToInt(columns[4]),
                                  columns[5], columns[6], columns[7], parseAttrs(columns[8]));
    }
};
    
/* Parse for GTF */
class GtfParser: public GxfParser {
    private:
    /* parse ID=ENSG00000223972.5 */
    static void parseAttr(const string& attrStr,
                          AttrVals& attrVals) {
        size_t i = attrStr.find(' ');
        if (i == string::npos) {
            throw invalid_argument("Invalid GTF attribute \"" + attrStr + "\"");
        }
        string name = attrStr.substr(0,i);
        string value = attrStr.substr(i+1);
        int idx = attrVals.findIdx(name);
        if (idx >= 0) {
            attrVals[idx]->addVal(stripQuotes(value));
        } else {
            attrVals.push_back(new AttrVal(name, stripQuotes(value)));
        }
    }

    /* parse: gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene";  */
    static AttrVals parseAttrs(const string& attrsStr) {
        AttrVals attrVals;
        StringVector parts = stringSplit(attrsStr,';');
        // last will be empty, since `;' is a terminator
        for (size_t i = 0; i < parts.size()-1; i++) {
            parseAttr(stringTrim(parts[i]),  attrVals);
        }
        return attrVals;
    }

    public:
   /* constructor */
    GtfParser(const string& fileName,
              GxfFeatureFactory gxfFeatureFactory):
        GxfParser(fileName, gxfFeatureFactory) {
    }
 
    /* get the format being parser */
    virtual GxfFormat getFormat() const {
        return GTF_FORMAT;
    }

     /* parse a feature */
    virtual GxfFeature* parseFeature(const StringVector& columns) {
        return fGxfFeatureFactory(columns[0], columns[1], columns[2],
                                  stringToInt(columns[3]), stringToInt(columns[4]),
                                  columns[5], columns[6], columns[7], parseAttrs(columns[8]));
    }
};

/* split a feature line of GFF3 or GTF */
StringVector GxfParser::splitFeatureLine(const string& line) const {
    StringVector columns = stringSplit(line, '\t');
    if (columns.size() != 9) {
        invalid_argument("invalid row, expected 9 columns: " + line);
    }
    return columns;
}

/* constructor that opens file, which maybe compressed. */
GxfParser::GxfParser(const string& fileName,
                     GxfFeatureFactory gxfFeatureFactory):
    fIn(new FIOStream(fileName)),
    fGxfFeatureFactory(gxfFeatureFactory) {
}

/* destructor */
GxfParser::~GxfParser() {
    delete fIn;
}

/* Read the next record */
GxfRecord* GxfParser::read() {
    string line;
    if (not fIn->readLine(line)) {
        return NULL;
    } else if ((line.size() > 0) and line[0] != '#') {
        return parseFeature(splitFeatureLine(line));
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

/* Factory to create a parser. file maybe compressed.  If gxfFormat is
 * unknown, guess from filename*/
GxfParser *GxfParser::factory(const string& fileName,
                              GxfFeatureFactory gxfFeatureFactory,
                              GxfFormat gxfFormat) {
    if (gxfFormat==GXF_UNKNOWN_FORMAT) {
        gxfFormat = gxfFormatFromFileName(fileName);
    }
    if (gxfFormat == GFF3_FORMAT) {
        return new Gff3Parser(fileName, gxfFeatureFactory);
    } else {
        return new GtfParser(fileName, gxfFeatureFactory);
    }
}

/* Write for GFF3 */
class Gff3Writer: public GxfWriter {
    public:
    /* format an attribute */
    static string formatAttr(const AttrVal* attrVal) {
        string strAttr = attrVal->getName() + "=";
        for (int i = 0; i < attrVal->getVals().size(); i++) {
            if (i > 0) {
                strAttr += ",";
            }
            strAttr += attrVal->getVals()[i];

        }
        return strAttr;
    }

    /* format attributes */
    static string formatAttrs(const AttrVals& attrVals) {
        string strAttrs;
        for (size_t i = 0; i < attrVals.size(); i++) {
            if (i > 0) {
                strAttrs += ";"; // separator
            }
            strAttrs += formatAttr(attrVals[i]);
        }
        return strAttrs;
    }

    /* constructor */
    Gff3Writer(const string& fileName):
        GxfWriter(fileName) {
        write("##gff-version 3");
    }

    /* get the format being parser */
    virtual GxfFormat getFormat() const {
        return GFF3_FORMAT;
    }

    /* format a feature line */
    virtual string formatFeature(const GxfFeature* feature) {
        return feature->baseColumnsAsString() + formatAttrs(feature->getAttrs());
    }
};

/* Write for GTF */
class GtfWriter: public GxfWriter {
    public:

    /* format an attribute */
    static string formatAttr(const string& name,
                             const string& val) {
        // n.b. this is not general, doesn't handle embedded quotes
        bool numericAttr = isNumeric(val);
        string strAttr = name + " ";
        if (!numericAttr) {
            strAttr += "\"";
        }
        strAttr += val;
        if (!numericAttr) {
            strAttr += "\"";
        }
        return strAttr;
    }

    /* format an attribute and values */
    static string formatAttr(const AttrVal* attrVal) {
        string strAttr;
        for (int i = 0; i < attrVal->getVals().size(); i++) {
            if (i > 0) {
                strAttr += " ";  // same formatting as GENCODE
            }
            strAttr += formatAttr(attrVal->getName(), attrVal->getVals()[i]) +  ";";
        }
        return strAttr;
    }

    /* should this attribute be included */
    static bool includeAttr(const AttrVal* attrVal) {
        // drop GFF3 linkage attributes
        return not ((attrVal->getName() == GxfFeature::ID_ATTR)
                    or (attrVal->getName() == GxfFeature::PARENT_ATTR)
                    or (attrVal->getName() == "remap_original_id"));
    }
    
    /* format attribute */
    static string formatAttrs(const AttrVals& attrVals) {
        string strAttrs;
        for (int i = 0; i < attrVals.size(); i++) {
            if (includeAttr(attrVals[i])) {
                if (strAttrs.size() > 0) {
                    strAttrs += " ";  // same formatting as GENCODE
                }
                strAttrs += formatAttr(attrVals[i]);
            }
        }
        return strAttrs;
    }
    /* constructor */
    GtfWriter(const string& fileName):
        GxfWriter(fileName) {
    }

    /* get the format being parser */
    virtual GxfFormat getFormat() const {
        return GTF_FORMAT;
    }

    /* format a feature line */
    virtual string formatFeature(const GxfFeature* feature) {
        return feature->baseColumnsAsString() + formatAttrs(feature->getAttrs());
    }
};

/* constructor that opens file */
GxfWriter::GxfWriter(const string& fileName):
    fOut(new FIOStream(fileName, ios::out)) {
}

/* destructor */
GxfWriter::~GxfWriter() {
    delete fOut;
}

/* Factory to create a writer. file maybe compressed.  If gxfFormat is
 * unknown, guess from filename*/
GxfWriter *GxfWriter::factory(const string& fileName,
                              GxfFormat gxfFormat) {
    if (gxfFormat==GXF_UNKNOWN_FORMAT) {
        gxfFormat = gxfFormatFromFileName(fileName);
    }
    if (gxfFormat == GFF3_FORMAT) {
        return new Gff3Writer(fileName);
    } else {
        return new GtfWriter(fileName);
    }
}

/* copy a file to output, normally used for a header */
void GxfWriter::copyFile(const string& inFile) {
    FIOStream inFh(inFile);
    string line;
    while (inFh.readLine(line)) {
        write(line);
    }
}

/* write one GxF record. */
void GxfWriter::write(const GxfRecord* gxfRecord) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        *fOut << formatFeature(dynamic_cast<const GxfFeature*>(gxfRecord)) << endl;
    } else {
        *fOut << gxfRecord->toString() << endl;
    }
}

/* write one GxF line. */
void GxfWriter::write(const string& line) {
    *fOut << line << endl;
}

/* return feature as a string */
string GxfFeature::toString() const {
    // just use GFF3 format, this is for debugging, not output
    return baseColumnsAsString() + Gff3Writer::formatAttrs(getAttrs());
}
