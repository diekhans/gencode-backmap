
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
        return s.substr(1, s.size()-1);
    } else {
        return s;
    }
}

/*
 * GFF3 Feature
 */
class Gff3Feature: public GxfFeature {
    private:
    // parse ID=ENSG00000223972.5
    AttrVal parseAttr(const string& attrStr) {
        size_t i = attrStr.find('=');
        if (i == string::npos) {
            throw invalid_argument("Invalid GFF3 attribute \"" + attrStr + "\"");
        }
        string name = attrStr.substr(0,i);
        string value = attrStr.substr(i+1);
        return AttrVal(name, stripQuotes(value), isQuoted(value));
    }
    
    AttrVals parseAttrs(const string& attrsStr) {
        AttrVals attrVals;
        StringVector parts = stringSplit(attrsStr,';');
        // last will be empty, since `;' is a terminator
        for (size_t i = 0; i < parts.size()-1; i++) {
            attrVals.push_back(parseAttr(parts[i]));
        }
        return attrVals;
    }

    public:
    /* constructor */
    Gff3Feature(const string& seqid, const string& source, const string& type, int
                start, int end, const string& score, const string& strand, const
                string& phase, const string& attrs):
        GxfFeature(seqid, source, type, start, end, score, strand, phase, parseAttrs(attrs)) {
    }
};


/*
 * GTF Feature
 */
class GtfFeature: public GxfFeature {
    private:
    // parse ID=ENSG00000223972.5
    AttrVal parseAttr(const string& attrStr) {
        size_t i = attrStr.find(' ');
        if (i == string::npos) {
            throw invalid_argument("Invalid GTF attribute \"" + attrStr + "\"");
        }
        string name = attrStr.substr(0,i);
        string value = attrStr.substr(i+1);
        return AttrVal(name, stripQuotes(value), isQuoted(value));
    }
    
    AttrVals parseAttrs(const string& attrsStr) {
        AttrVals attrVals;
        StringVector parts = stringSplit(attrsStr,';');
        // last will be empty, since `;' is a terminator
        for (size_t i = 0; i < parts.size()-1; i++) {
            attrVals.push_back(parseAttr(parts[i]));
        }
        return attrVals;
    }

    public:
    /* constructor */
    GtfFeature(const string& seqid, const string& source, const string& type, int
               start, int end, const string& score, const string& strand, const
               string& phase, const string& attrs):
        GxfFeature(seqid, source, type, start, end, score, strand, phase, parseAttrs(attrs)) {
    }
};


/* split a feature line of GFF3 or GTF */
StringVector GxfParser::splitFeatureLine(const string& line) {
    StringVector columns = stringSplit(line, '\t');
    if (columns.size() != 9) {
        invalid_argument("invalid row, expected 9 columns: " + line);
    }
    return columns;
}

/* create a gff3 feature */
const GxfFeature* GxfParser::createGff3Feature(const StringVector& columns) {
    return new Gff3Feature(columns[0], columns[1], columns[2],
                           stringToInt(columns[3]), stringToInt(columns[4]),
                           columns[5], columns[6], columns[7], columns[8]);
}

/* create a gtf feature */
const GxfFeature* GxfParser::createGtfFeature(const StringVector& columns) {
    return new GtfFeature(columns[0], columns[1], columns[2],
                          stringToInt(columns[3]), stringToInt(columns[4]),
                          columns[5], columns[6], columns[7], columns[8]);
}

/* create the appropriate feature type */
const GxfFeature* GxfParser::createGxfFeature(const StringVector& columns) {
    if (fGxfFormat == GFF3_FORMAT) {
        return createGff3Feature(columns);
    } else {
        return createGxfFeature(columns);
    }
}

/* constructor that opens file, which maybe compressed */
GxfParser::GxfParser(const string& fileName,
                     GxfFormat gxfFormat):
    fIn(new FIOStream(fileName)),
    fGxfFormat(gxfFormat) {
}

/* destructor */
GxfParser::~GxfParser() {
    delete fIn;
}

/* Read the next record */
const GxfRecord* GxfParser::read() {
    string line;
    if (not fIn->readLine(line)) {
        return NULL;
    } else if ((line.size() > 0) and line[0] != '#') {
        return createGxfFeature(splitFeatureLine(line));
    } else {
        return new GxfLine(line);
    }
}

/* Read the next record, either queued by push() or from the file , use
 * instanceOf to determine the type.  Return NULL on EOF.
 */
const GxfRecord* GxfParser::next() {
    if (not fPending.empty()) {
        const GxfRecord* rec = fPending.front();
        fPending.pop();
        return rec;
    } else {
        return read();
    }
}
