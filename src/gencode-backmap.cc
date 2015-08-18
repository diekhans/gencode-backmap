/*
 * program to map gencode datafiles to older assemblies.
 */
#include "jkinclude.hh"
#include <getopt.h>
#include "gxf.hh"
#include "typeOps.hh"
#include <fstream>
#include <iostream>
#include "transMap.hh"
#include "geneMapper.hh"

/* map to different assembly */
static void gencodeBackmap(const string& inGxfFile,
                           GxfFormat gxfFormat,
                           const string& mappingChains,
                           bool swapMap,
                           const string& outGxfFile) {
    TransMap* genomeTransMap = TransMap::chainFileFactory(mappingChains, swapMap);
    GxfParser gxfParser(inGxfFile, gxfFormat);
    ofstream outGxf(outGxfFile);
    GeneMapper geneMapper(genomeTransMap);
    geneMapper.mapGxf(&gxfParser, outGxf);
    delete genomeTransMap;
}

/* Entry point.  Parse arguments. */
int main(int argc, char *argv[]) {
    const struct option long_options[] = {
        {"swapMap", 0, NULL, 's'},
        {NULL, 0, NULL, 0}
    };
    bool swapMap = false;
    opterr = 0;  // we print error message
    while (true) {
        int optc = getopt_long(argc, argv, "s", long_options, NULL);
        if (optc == -1) {
            break;
        } else if (optc == 's') {
            swapMap = true;
        } else {
            errAbort(toCharStr("invalid option %s"), argv[optind-1]);
        }
    }    

    if ((argc - optind) != 3) {
        errAbort(toCharStr("wrong # args: %s [-swapMap] inGxf mappingChains outGxf"), argv[0]);
    }
    string inGxfFile = argv[optind];
    string mappingChains = argv[optind+1];
    string outGxfFile = argv[optind+2];

    GxfFormat gxfFormat = (stringEndsWith(inGxfFile, ".gff3") or stringEndsWith(inGxfFile, ".gff3.gz"))
        ? GFF3_FORMAT : GTF_FORMAT;
    try {
        gencodeBackmap(inGxfFile, gxfFormat, mappingChains, swapMap, outGxfFile);
    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }
    return 0;
}
