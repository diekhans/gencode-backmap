/*
 * program to map gencode datafiles to older assemblies.
 */
#include "jkinclude.hh"
#include <getopt.h>
#include "gxf.hh"
#include "typeOps.hh"
#include "FIOStream.hh"
#include "transMap.hh"
#include "geneMapper.hh"
#include "targetAnnotations.hh"

/* map to different assembly */
static void gencodeBackmap(const string& inGxfFile,
                           const string& mappingAligns,
                           bool swapMap,
                           const string& substituteMissingTargetVersion,
                           const string& mappedGxfFile,
                           const string& unmappedGxfFile,
                           const string& mappingInfoTsv,
                           const string& targetGxf,
                           const string& transcriptPsls) {
    TransMap* genomeTransMap = TransMap::factoryFromFile(mappingAligns, swapMap);
    TargetAnnotations* targetAnnotations = (targetGxf.size() > 0) ? new TargetAnnotations(targetGxf) : NULL;
    GxfParser gxfParser(inGxfFile);
    FIOStream mappedGxfFh(mappedGxfFile, ios::out);
    FIOStream unmappedGxfFh(unmappedGxfFile, ios::out);
    FIOStream mappingInfoFh(mappingInfoTsv, ios::out);
    FIOStream *transcriptPslFh = (transcriptPsls.size() > 0) ? new FIOStream(transcriptPsls, ios::out) : NULL;
    GeneMapper geneMapper(genomeTransMap, targetAnnotations, substituteMissingTargetVersion);
    geneMapper.mapGxf(&gxfParser, mappedGxfFh, unmappedGxfFh, mappingInfoFh, transcriptPslFh);
    delete genomeTransMap;
    delete targetAnnotations;
    delete transcriptPslFh;
}

/* Entry point.  Parse arguments. */
int main(int argc, char *argv[]) {
    const string usage = "%s [options] target inGxf mappingAligns mappedGxf unmappedGxf mappingInfoTsv\n\n"
        "Map GENCODE annotations between assemblies projecting through genomic\n"
        "alignments. This operates on GENCODE GFF3 and GTF files and makes assumptions\n"
        "about their organization.\n\n"
        "Options:\n"
        "  --help - print this message and exit\n"
        "  --swapMap - swap the query and target sides of the mapping alignments\n"
        "  --targetGxf=gxfFile - GFF3 or GTF of gene annotations on target genome.\n"
        "    If specified, require mappings to location of previous version of\n"
        "    gene or transcript.\n"
        "  --transcriptPsls=pslFile - write all mapped transcript-level PSL to this file, including\n"
        "    multiple mappers.\n"
        "  --substituteMissingTargets=targetVersion - if target GxF is specified and no GENE maps to\n"
        "    the target locus, pass through the original target location.  Only a subset of the\n"
        "    biotypes are substituted. Argument is target GENCODE version that is stored as an attribute\n"
        "Arguments:\n"
        "  inGxf - Input GENCODE GFF3 or GTF file. The format is identified\n"
        "          by a .gff3 or .gtf extension, it maybe compressed with gzip with an\n"
        "          additional .gz extensionn\n"
        "  mappingAligns - Alignments between the two genomes.  The \n"
        "  mappedGxf - GxF file of mapped features on target genome\n"
        "  unmappedGxf - GxF file of unmapped features on source genome\n"
        "  mappingInfoTsv - TSV file with information about each gene and transcript mapping\n";

    const struct option long_options[] = {
        {"help", 0, NULL, 'h'},
        {"swapMap", 0, NULL, 's'},
        {"targetGxf", 1, NULL, 't'},
        {"transcriptPsls", 1, NULL, 'p'},
        {"substituteMissingTargets", 1, NULL, 'm'},
        {NULL, 0, NULL, 0}
    };
    const char* short_options = "hst:p:m:";
    
    bool swapMap = false;
    bool help = false;
    string targetGxf;
    string transcriptPsls;
    string substituteMissingTargetVersion;
    opterr = 0;  // we print error message
    while (true) {
        int optc = getopt_long(argc, argv, short_options, long_options, NULL);
        if (optc == -1) {
            break;
        } else if (optc == 'h') {
            help = true;
            break;  // check no more
        } else if (optc == 's') {
            swapMap = true;
        } else if (optc == 't') {
            targetGxf = string(optarg);
        } else if (optc == 'p') {
            transcriptPsls = string(optarg);
        } else if (optc == 'm') {
            substituteMissingTargetVersion = string(optarg);
        } else {
            errAbort(toCharStr("invalid option %s"), argv[optind-1]);
        }
    }    
    if (help) {
        // don't check any other options with help
        cerr << usage << endl;
        return 1;
    }

    if ((argc - optind) != 5) {
        errAbort(toCharStr("wrong # args: " + usage), argv[0]);
    }
    string inGxfFile = argv[optind];
    string mappingAligns = argv[optind+1];
    string mappedGxfFile = argv[optind+2];
    string unmappedGxfFile = argv[optind+3];
    string mappingInfoTsv = argv[optind+4];

    try {
        gencodeBackmap(inGxfFile, mappingAligns, swapMap, substituteMissingTargetVersion, mappedGxfFile, unmappedGxfFile, mappingInfoTsv, targetGxf, transcriptPsls);
    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }
    return 0;
}
