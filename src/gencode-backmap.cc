/*
 * program to map gencode datafiles to older assemblies.
 */
#include "jkcommon.hh"
#include "gxf.hh"
#include "typeOps.hh"
#include <fstream>
#include "transMap.hh"
#include "geneMapper.hh"

/* map to different assembly */
static void gencodeBackmap(const string& inGxfFile,
                           GxfFormat gxfFormat,
                           const string& mappingChains,
                           const string& outGxfFile) {
    TransMap transMap(mappingChains, false);
    GxfParser gxfParser(inGxfFile, gxfFormat);
    ofstream outGxf(outGxfFile);
    GeneMapper geneMapper(&transMap);
    geneMapper.mapGxf(&gxfParser, outGxf);
}

/* Entry point.  Parse arguments. */
int main(int argc, char *argv[]) {
    if (argc != 4) {
        errAbort(toCharStr("wrong # args: %s inGxf mappingChains outGxf"), argv[0]);
    }
    string inGxfFile = argv[1];
    GxfFormat gxfFormat = stringEndsWith(inGxfFile, ".gff3") ? GFF3_FORMAT : GTF_FORMAT;
    string mappingChains = argv[2];
    string outGxfFile = argv[3];

    gencodeBackmap(inGxfFile, gxfFormat, mappingChains, outGxfFile);

    return 0;
}
