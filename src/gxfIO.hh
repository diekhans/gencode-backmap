/*
 * Extremely naive and specialized GFF3 and GTF parsers.
 *
 * These parser assume the ordering of the GENCODE GFF3/GTF files.
 */

#ifndef gxfIO_hh
#define gxfIO_hh
#include "typeOps.hh"
#include <queue>
#include <algorithm>
using namespace std;

/* it seems so stupid to need to keep writing one-off GFF/GTF parsers */

class GxfFeature;
class GxfRecord;
class AttrVals;
class FIOStream;

typedef enum {
    GXF_UNKNOWN_FORMAT,
    GFF3_FORMAT,
    GTF_FORMAT,
} GxfFormat;

/* Get format from file name, or error */
GxfFormat gxfFormatFromFileName(const string& fileName);

/*
 * Factory function use to obtain a feature
 */
typedef GxfFeature* (*GxfFeatureFactory)(const string& seqid, const string& source, const string& type,
                                         int start, int end, const string& score, const string& strand,
                                         const string& phase, const AttrVals& attrs);
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
    GxfFeatureFactory fGxfFeatureFactory;
    
    /* parse a feature */
    virtual GxfFeature* parseFeature(const StringVector& columns) = 0;
    
    /* constructor that opens file */
    GxfParser(const string& fileName,
              GxfFeatureFactory gxfFeatureFactory);

    public:
    /* destructor */
    virtual ~GxfParser();

    /* get the format being parser */
    virtual GxfFormat getFormat() const = 0;

    /* Factory to create a parser. file maybe compressed.  If gxfFormat is
     * unknown, guess from filename*/
    static GxfParser *factory(const string& fileName,
                              GxfFeatureFactory gxfFeatureFactory,
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
                              GxfFormat gxfFormat=GXF_UNKNOWN_FORMAT);

    /* copy a file to output, normally used for a header */
    void copyFile(const string& inFile);

    /* write one GxF record. */
    void write(const GxfRecord* gxfRecord);

    /* write one GxF line. */
    void write(const string& line);
    
};
#endif
