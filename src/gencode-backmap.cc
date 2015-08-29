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
                           const string& mappingAligns,
                           bool chainAlignments,
                           bool swapMap,
                           const string& mappedGxfFile,
                           const string& unmappedGxfFile) {
    TransMap* genomeTransMap = chainAlignments
        ? TransMap::factoryFromChainFile(mappingAligns, swapMap)
        : TransMap::factoryFromPslFile(mappingAligns, swapMap);

    GxfParser gxfParser(inGxfFile, gxfFormat);
    ofstream mappedGxfFh(mappedGxfFile);
    ofstream unmappedGxfFh(unmappedGxfFile);
    GeneMapper geneMapper(genomeTransMap);
    geneMapper.mapGxf(&gxfParser, mappedGxfFh, unmappedGxfFh);
    delete genomeTransMap;
}

/* Entry point.  Parse arguments. */
int main(int argc, char *argv[]) {
    const string usage = "%s [options] inGxf mappingAligns mappedGxf unmappedGxf\n\n"
        "Map GENCODE annotations between assemblies projecting through genomic\n"
        "alignments. This operates on GENCODE GFF3 and GTF files and makes assumptions\n"
        "about their organization.\n\n"
        "Options:\n"
        "  --help - print this message and exit\n"
        "  --swapMap - swap the query and target sides of the mapping alignments\n"
        "Arguments:\n"
        "  inGxf - Input GENCODE GFF3 or GTF file. The format is recognize identified\n"
        "          a .gff3 or .gtf extension, it maybe compressed with gzip with an\n"
        "          additional .gz extensionn\n"
        "  mappingAligns - Alignments between the two genomes.  The \n"
        "  mappedGxf - \n"
        "  unmappedGxf - \n"

;

    const struct option long_options[] = {
        {"help", 0, NULL, 'h'},
        {"swapMap", 0, NULL, 's'},
        {NULL, 0, NULL, 0}
    };
    bool swapMap = false;
    bool help = false;
    opterr = 0;  // we print error message
    while (true) {
        int optc = getopt_long(argc, argv, "sch", long_options, NULL);
        if (optc == -1) {
            break;
        } else if (optc == 'h') {
            help = true;
            break;  // check no more
        } else if (optc == 's') {
            swapMap = true;
        } else {
            errAbort(toCharStr("invalid option %s"), argv[optind-1]);
        }
    }    
    if (help) {
        // don't check any other options with help
        cerr << usage << endl;
        return 1;
    }

    if ((argc - optind) != 4) {
        errAbort(toCharStr("wrong # args: " + usage), argv[0]);
    }
    string inGxfFile = argv[optind];
    string mappingAligns = argv[optind+1];
    string mappedGxfFile = argv[optind+2];
    string unmappedGxfFile = argv[optind+3];

    GxfFormat gxfFormat = UNKNOWN_FORMAT;
    if (stringEndsWith(inGxfFile, ".gff3") or stringEndsWith(inGxfFile, ".gff3.gz")) {
        gxfFormat = GFF3_FORMAT;
    } else if (stringEndsWith(inGxfFile, ".gtf") or stringEndsWith(inGxfFile, ".gtf.gz")) {
        gxfFormat = GTF_FORMAT;
    } else {
        errAbort(toCharStr("Error: expected input annotation with an extension of .gff3, .gff3.gz, .gtf, or .gtf.gz"));
    }
    
    bool chainAlignments = false;
    if (stringEndsWith(mappingAligns, ".chain") or stringEndsWith(mappingAligns, ".chain.gz")) {
        chainAlignments = true;
    } else if (stringEndsWith(mappingAligns, ".psl") or stringEndsWith(mappingAligns, ".psl.gz")) {
         chainAlignments = false;
    } else {
        errAbort(toCharStr("Error: expected mapping alignment file with an extension of .chain, .chain.gz, .psl, or .psl.gz"));
    }
    try {
        gencodeBackmap(inGxfFile, gxfFormat, mappingAligns, chainAlignments, swapMap, mappedGxfFile, unmappedGxfFile);
    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }
    return 0;
}
