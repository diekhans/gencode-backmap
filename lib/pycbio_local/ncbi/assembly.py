# Copyright 2015-2015 Mark Diekhans
"""Parsing of NCBI assembly information files.
"""
import six
from pycbio_local.sys import PycbioException


def _noneIfNa(name):
    return None if name == "na" else name


def _naIfNone(name):
    return "na" if name is None else name


class AssemblyReport(object):
    """Parse assembly reports files, e.g.
      ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_report.txt
    Sequence ids of 'na' are converted to None.
    """

    class Record(object):
        def __init__(self, sequenceName, sequenceRole, assignedMolecule, locationType, genBankAccn, relationship, refSeqAccn, assemblyUnit, sequenceLength, ucscStyleName):
            self.sequenceName = sequenceName
            self.sequenceRole = sequenceRole
            self.assignedMolecule = assignedMolecule
            self.locationType = locationType
            self.genBankAccn = _noneIfNa(genBankAccn)
            self.relationship = relationship
            self.refSeqAccn = _noneIfNa(refSeqAccn)
            self.assemblyUnit = assemblyUnit
            self.sequenceLength = sequenceLength
            self.ucscStyleName = _noneIfNa(ucscStyleName)

        @property
        def gencodeName(self):
            """GENCODE uses UCSC name for primary chromosomes and mitochondia and GENBANK accessions
            for others"""
            if self.sequenceRole == "assembled-molecule":
                return self.ucscStyleName
            else:
                return self.genBankAccn

        @property
        def sequenceNameGenbankAccn(self):
            """Sequence name primary chromosomes and mitochondia and GENBANK accessions
            for others"""
            if self.sequenceRole == "assembled-molecule":
                return self.sequenceName
            else:
                return self.genBankAccn

        def __str__(self):
            return "\t".join([self.sequenceName, self.sequenceRole, self.assignedMolecule, self.locationType,
                              _naIfNone(self.genBankAccn), self.relationship, _naIfNone(self.refSeqAccn), self.assemblyUnit, str(self.sequenceLength), _naIfNone(self.ucscStyleName)])

    expectedHeader = "# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name"

    def __init__(self, asmReport):
        self.metaData = dict()  # main headers at started
        self.seqs = []
        self.bySequenceName = dict()
        self.byGenBankAccn = dict()
        self.byRefSeqAccn = dict()
        self.byUcscStyleName = dict()
        mode = "r" if six.PY3 else "rU"
        with open(asmReport, mode) as fh:
            self._parseMetaData(fh)
            self._skipToSeqTable(fh)
            self._parseRecords(fh)

    def _parseMetaData(self, fh):
        "parse metaData records at start"
        # Assembly Name:  GRCh38.p2
        for line in fh:
            line = line[0:-1]
            if line == "#":
                break
            self._parseMetaDataLine(line)

    def _parseMetaDataLine(self, line):
        colon = line.find(":")
        if colon < 0:
            raise PycbioException("invalid metaData line: " + line)
        self.metaData[line[2:colon]] = line[colon + 1:].strip()

    def _skipToSeqTable(self, fh):
        "skip past header line before sequence records"
        for line in fh:
            if line[0:-1] == self.expectedHeader:
                return
        raise PycbioException("expected assembly report header not found in " + fh.name)

    def _parseRecords(self, fh):
        for line in fh:
            self._parseRecord(fh, line[0:-1])

    def _parseRecord(self, fh, line):
        row = line.split('\t')
        if len(row) != 10:
            raise PycbioException("expected 10 columns in assemble report record, found " + str(len(row)) + " in " + fh.name)
        rec = self.Record(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], int(row[8]), row[9])
        self.seqs.append(rec)
        self.bySequenceName[rec.sequenceName] = rec
        if rec.genBankAccn is not None:
            self.byGenBankAccn[rec.genBankAccn] = rec
        if rec.refSeqAccn is not None:
            self.byRefSeqAccn[rec.refSeqAccn] = rec
        if rec.ucscStyleName is not None:
            self.byUcscStyleName[rec.ucscStyleName] = rec

    def ucscNameToSeqName(self, ucscName):
        try:
            return self.byUcscStyleName[ucscName].sequenceName
        except KeyError:
            raise ValueError("unknown UCSC chromosome name: '{}'".format(ucscName))

    def seqNameToUcscName(self, seqName):
        try:
            return self.bySequenceName[seqName].ucscStyleName
        except KeyError:
            raise ValueError("unknown chromosome name: '{}'".format(seqName))
