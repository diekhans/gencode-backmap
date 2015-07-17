
#include "gxf.hh"
extern "C" {
#define min jkmin
#define max jkmax
#include "common.h"
#undef min
#undef max
}
#include "FIOStream.hh"
#include <typeinfo>
#include <string>
#include <vector>
#include "typeOps.hh"
#include <stdexcept>


/*
 * GFF3 Feature
 */
class Gff3Feature: public Feature {
    public:
    /* constructor */
    Gff3Feature(const string& seqid, const string& source, const string& type, int
                start, int end, const string& score, const string& strand, const
                string& phase, const string& attrs):
        Feature(seqid, source, type, start, end, score, strand, phase, attrs) {
    }
};


/*
 * GTF Feature
 */
class GtfFeature: public Feature {
    public:
    /* constructor */
    GtfFeature(const string& seqid, const string& source, const string& type, int
               start, int end, const string& score, const string& strand, const
               string& phase, const string& attrs):
        Feature(seqid, source, type, start, end, score, strand, phase, attrs) {
    }
};


/* split a feature line of GFF3 or GTF */
vector<const string> GxfParser::splitFeatureLine(const string& line) {
    vector<const string> columns = stringSplit(line, '\t');
    if (columns.size() != 9) {
        invalid_argument("invalid row, expected 9 columns: " + line);
    }
    return columns;
}

/* create a gff3 feature */
const Feature* GxfParser::createGff3Feature(const vector<const string>& columns) {
    return new Gff3Feature(columns[0], columns[1], columns[2],
                           stringToInt(columns[3]), stringToInt(columns[4]),
                           columns[5], columns[6], columns[7], columns[8]);
}

/* create a gtf feature */
const Feature* GxfParser::createGtfFeature(const vector<const string>& columns) {
    return new GtfFeature(columns[0], columns[1], columns[2],
                          stringToInt(columns[3]), stringToInt(columns[4]),
                          columns[5], columns[6], columns[7], columns[8]);
}

/* create the appropriate feature type */
const Feature* GxfParser::createGxfFeature(const vector<const string>& columns) {
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

/* Read the next line, parse into a feature if it is one, otherwise just
 * return line. Return NULL on EOF */
GxfLine* GxfParser::next() {
    string line;
    if (not fIn->readLine(line)) {
        return NULL;
    } else if ((line.size() > 0) and line[0] != '#') {
        return new GxfLine(line, createGxfFeature(splitFeatureLine(line)));
    } else {
        return new GxfLine(line);
    }
}
