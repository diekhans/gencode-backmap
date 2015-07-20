/*
 * program to map gencode datafiles to older assemblies.
 */
#include "jkcommon.hh"
#include "gxf.hh"
#include "typeOps.hh"
#include <fstream>
#include "geneMapper.hh"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        errAbort(toCharStr("wrong # args: %s inGxf outGxf"), argv[0]);
    }
    string inGxfFile = argv[1];
    string outGxfFile = argv[2];
    
    GxfParser gxfParser(inGxfFile,
                        stringEndsWith(inGxfFile, ".gff3") ? GFF3_FORMAT : GTF_FORMAT);
    ofstream outGxf(outGxfFile);
    GeneMapper geneMapper(NULL);
    geneMapper.mapGxf(&gxfParser, outGxf);
    return 0;
}
