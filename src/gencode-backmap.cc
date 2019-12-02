
/*
 * program to map gencode datafiles to older assemblies.
 */
#include <getopt.h>
#include "gxf.hh"
#include "typeOps.hh"
#include "FIOStream.hh"
#include "transMap.hh"
#include "geneMapper.hh"
#include "annotationSet.hh"
#include "bedMap.hh"
#include "globals.hh"
#include "gxf.hh"
#include "./version.h"

/* verbose tracing enabled */
bool gVerbose = false;

/* check if a file is in the same format as a the input */
static bool checkGxfFormat(GxfFormat inFormat,
                           const string& gxfFile,
                           bool isOptional) {
    if ((gxfFile == "") and isOptional) {
        return true;
    }
    GxfFormat gxfFormat = gxfFormatFromFileName(gxfFile);
    if (not ((gxfFormat == inFormat) or (gxfFormat == DEV_NULL_FORMAT))) {
        cerr << "Error: all input and output formats must be consistently GFF3 or GTF" << endl;
        return false;
    } else {
        return true;
    }
}

/* check all files are in the same format */
static bool checkGxfFormats(const string& inGxfFile,
                            const string& mappedGxfFile,
                            const string& targetGxf,
                            const string& previousMappedGxf) {
    GxfFormat inFormat = gxfFormatFromFileName(inGxfFile);
    return checkGxfFormat(inFormat, mappedGxfFile, false)
        and checkGxfFormat(inFormat, targetGxf, true)
        and checkGxfFormat(inFormat, previousMappedGxf, true);
}

/* map to different assembly */
static void gencodeBackmap(const string& inGxfFile,
                           const string& mappingAligns,
                           bool swapMap,
                           const string& substituteMissingTargetVersion,
                           unsigned useTargetFlags,
                           bool onlyManualForTargetSubstituteOverlap,
                           ParIdHackMethod parIdHackMethod,
                           const string& headerFile,
                           const string& mappedGxfFile,
                           const string& mappingInfoTsv,
                           const string& targetGxf,
                           const string& targetPatchBed,
                           const string& previousMappedGxf,
                           const string& transcriptPsls) {
    TransMap* genomeTransMap = TransMap::factoryFromFile(mappingAligns, swapMap);
    AnnotationSet srcAnnotations(inGxfFile);
    AnnotationSet* targetAnnotations = (targetGxf.size() > 0)
        ? new AnnotationSet(targetGxf) : NULL;
    AnnotationSet* previousMappedAnnotations = (previousMappedGxf.size() > 0)
        ? new AnnotationSet(previousMappedGxf) : NULL;
    BedMap* targetPatchMap = (targetPatchBed.size() > 0)
        ? new BedMap(targetPatchBed) : NULL;
    GxfWriter* mappedGxfFh = GxfWriter::factory(mappedGxfFile, parIdHackMethod);
    if (headerFile.size() > 0) {
        mappedGxfFh->copyFile(headerFile);
    }
    FIOStream mappingInfoFh((mappingInfoTsv.size() > 0) ? mappingInfoTsv : "/dev/null" , ios::out);
    FIOStream* transcriptPslFh = (transcriptPsls.size() > 0) ? new FIOStream(transcriptPsls, ios::out) : NULL;
    GeneMapper geneMapper(&srcAnnotations, genomeTransMap, targetAnnotations, previousMappedAnnotations,
                          targetPatchMap, substituteMissingTargetVersion,
                          useTargetFlags, onlyManualForTargetSubstituteOverlap);
    geneMapper.mapGxf(*mappedGxfFh, mappingInfoFh, transcriptPslFh);
    delete mappedGxfFh;
    delete genomeTransMap;
    delete targetPatchMap;
    delete targetAnnotations;
    delete previousMappedAnnotations;
    delete transcriptPslFh;
}

const string usage = "%s [options] inGxf mappingAligns mappedGxf [mappingInfoTsv]\n\n"
    "Map GENCODE annotations between assemblies projecting through genomic\n"
    "alignments. This operates on GENCODE GFF3 and GTF files and makes assumptions\n"
    "about their organization.\n\n"
    "Options:\n"
    "  --help - print this message and exit\n"
    "  --verbose - verbose tracing to stderr\n"
    "  --swapMap - swap the query and target sides of the mapping alignments\n"
    "  --targetGxf=gxfFile - GFF3 or GTF of gene annotations on target genome.\n"
    "    If specified, require mappings to location of previous version of\n"
    "    gene or transcript.\n"
    "  --previousMappedGxf=gxfFile - GFF3 or GTF of gene annotations on previous mapping.\n"
    "    This is used to determine the mapped version number to append to the ids.\n"
    "  --headerFile=commentFile - copy contents of this file as comment header for GFF3/GTF output.\n"
    "    Doesn't include GFF3 file type meta comment.\n"
    "  --transcriptPsls=pslFile - write all mapped transcript-level PSL to this file, including\n"
    "    multiple mappers.\n"
    "  --substituteMissingTargets=targetVersion - if target GxF is specified and no GENE maps to\n"
    "    the target locus, pass through the original target location.  Only a subset of the\n"
    "    biotypes are substituted. Argument is target GENCODE version that is stored as an attribute\n"
    "  --targetPatches=targetPatchRegionBed Target genes with no corresponding mappings and that overlap\n"
    "    thee patched regions in the target genome will be passed through.\n"
    "  --useTargetForAutoSmallNonCoding - don't map automatic small non-coding transcripts, substituting\n"
    "    the target if requested.\n"
    "  --useTargetForAutoGenes - don't map automatic-only genes, substituting\n"
    "    the target automatic genes, even if they are not in the source \n"
    "    If a target set is not specified, they are skipped.\n"
    "    This does not work well in practice, left in for experimentation.\n"
    "  --useTargetForPseudoGenes - don't map pseudo genes, substituting\n"
    "    the target genes, even if they are not in the source.\n"
    "    If a target set is not specified, they are skipped.\n"
    "    This does not work well in practice, left in for experimentation.\n"
    "  --onlyManualForTargetSubstituteOverlap - when checking for overlap of\n"
    "    target with mapped before substituting a target gene, only consider\n"
    "    manual transcripts.\n"
    "  --oldStyleParIdHack - use ENSTR style PAR id unique on output rather than the\n"
    "    newer _PAR_Y.  Either form is recognized on input.\n"
    "Arguments:\n"
    "  inGxf - Input GENCODE GFF3 or GTF file. The format is identified\n"
    "          by a .gff3 or .gtf extension, it maybe compressed with gzip with an\n"
    "          additional .gz extensionn.  All GxF files types must be consistent;\n"
    "          either all GFF3 or all GTF.\n"
    "  mappingAligns - Alignments between the two genomes.  The \n"
    "  mappedGxf - GxF file of mapped features on target genome\n"
    "  mappingInfoTsv - TSV file with information about each gene and transcript mapping\n"
    "\n";

const struct option long_options[] = {
    {"help", 0, NULL, 'h'},
    {"verbose", 0, NULL, 'v'},
    {"swapMap", 0, NULL, 's'},
    {"targetGxf", 1, NULL, 't'}, 
    {"previousMappedGxf", 1, NULL, 'M'}, 
    {"targetPatches", 1, NULL, 'T'}, 
    {"headerFile", 1, NULL, 'H'},
    {"transcriptPsls", 1, NULL, 'p'},
    {"substituteMissingTargets", 1, NULL, 'm'},
    {"useTargetForAutoSmallNonCoding", 0, NULL, 'N'},
    {"useTargetForAutoGenes", 0, NULL, 'A'},
    {"useTargetForPseudoGenes", 0, NULL, 'P'},
    {"onlyManualForTargetSubstituteOverlap", 0, NULL, 'O'},
    {"oldStyleParIdHack", 0, NULL, 'Q'},
    {NULL, 0, NULL, 0}
};
const char* short_options = "hst:p:m:n";

/* print usage */
static void prUsage() {
    cerr << usage << "Version: " << VERSION << " (" << VERSION_HASH <<  ")" << endl;
}

/* Entry point.  Parse arguments. */
int main(int argc, char *argv[]) {
    bool swapMap = false;
    bool help = false;
    unsigned useTargetFlags = 0;
    string targetGxf;
    string targetPatchBed;
    string headerFile;
    string previousMappedGxf;
    string transcriptPsls;
    string substituteMissingTargetVersion;
    ParIdHackMethod parIdHackMethod = PAR_ID_HACK_NEW;
    bool onlyManualForTargetSubstituteOverlap = false;
    opterr = 0;  // we print error message
    while (true) {
        int optc = getopt_long(argc, argv, short_options, long_options, NULL);
        if (optc == -1) {
            break;
        } else if (optc == 'h') {
            help = true;
            break;  // check no more
        } else if (optc == 'v') {
            gVerbose = true;
        } else if (optc == 's') {
            swapMap = true;
        } else if (optc == 't') {
            targetGxf = string(optarg);
        } else if (optc == 'T') {
            targetPatchBed = string(optarg);
            useTargetFlags |= GeneMapper::useTargetForPatchRegions;
        } else if (optc == 'H') {
            headerFile = string(optarg);
        } else if (optc == 'M') {
            previousMappedGxf = string(optarg);
        } else if (optc == 'p') {
            transcriptPsls = string(optarg);
        } else if (optc == 'm') {
            substituteMissingTargetVersion = string(optarg);
        } else if (optc == 'O') {
            onlyManualForTargetSubstituteOverlap = true;
        } else if (optc == 'N') {
            useTargetFlags |= GeneMapper::useTargetForAutoNonCoding;
        } else if (optc == 'A') {
            useTargetFlags |= GeneMapper::useTargetForAutoGenes;
        } else if (optc == 'P') {
            useTargetFlags |= GeneMapper::useTargetForPseudoGenes;
        } else if (optc == 'Q') {
            parIdHackMethod = PAR_ID_HACK_OLD;
        } else {
            errAbort(toCharStr("invalid option %s"), argv[optind-1]);
        }
    }    
    if (help) {
        // don't check any other options with help
        prUsage();
        return 1;
    }

    int nposargs = (argc - optind);
    if ((nposargs < 3) or (nposargs > 4)) {
        cerr << "wrong # args: ";
        prUsage();
        return 1;
    }
    string inGxfFile = argv[optind];
    string mappingAligns = argv[optind+1];
    string mappedGxfFile = argv[optind+2];
    string mappingInfoTsv = (nposargs > 3) ? argv[optind+3] : "";

    if (not checkGxfFormats(inGxfFile, mappedGxfFile, targetGxf, previousMappedGxf)) {
        return 1;
    }
    
    try {
        gencodeBackmap(inGxfFile, mappingAligns, swapMap,
                       substituteMissingTargetVersion, useTargetFlags,
                       onlyManualForTargetSubstituteOverlap, parIdHackMethod,
                       headerFile, mappedGxfFile, mappingInfoTsv,
                       targetGxf, targetPatchBed, previousMappedGxf,
                       transcriptPsls);
    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }
    return 0;
}
