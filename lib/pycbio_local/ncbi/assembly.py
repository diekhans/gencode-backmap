# Copyright 2006-2022 Mark Diekhans
"""Parsing of NCBI assembly information files.
"""
from collections import namedtuple

from pycbio import PycbioException
from pycbio.sys import fileOps
from pycbio.sys.objDict import ObjDict


def _noneIfNa(name):
    return None if name == "na" else name


def _naIfNone(name):
    return "na" if name is None else name

class Record(namedtuple("Record",
                        ("sequenceName", "sequenceRole", "assignedMolecule", "locationType", "genBankAccn",
                         "relationship", "refSeqAccn", "assemblyUnit", "sequenceLength", "ucscStyleName"))):
    """Sequence record"""

    @property
    def gencodeName(self):
        """GENCODE uses UCSC name for primary chromosomes and mitochondia and GENBANK accessions
        for others. """
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
def _parseRecord(row):
    if len(row) != 10:
        raise PycbioException("expected 10 columns in assembly report record, found: " + str(len(row)))
    return Record(sequenceName=row[0],
                  sequenceRole=row[1],
                  assignedMolecule=row[2],
                  locationType=row[3],
                  genBankAccn=_noneIfNa(row[4]),
                  relationship=row[5],
                  refSeqAccn=_noneIfNa(row[6]),
                  assemblyUnit=row[7],
                  sequenceLength=int(row[8]),
                  ucscStyleName=_noneIfNa(row[9]))

class MetaData(ObjDict):
    "main headers at start, as fields, with ' ', converted to '_'"
    pass

def _parseMetaDataLine(metaData, line):
    colon = line.find(":")
    if colon < 0:
        raise PycbioException("invalid metaData line: " + line)
    name = line[2:colon].strip().lower().replace(' ', '_')
    value = line[colon + 1:].strip()
    metaData[name] = value

class AssemblyReport:
    """Parse assembly reports files, e.g.
      ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_report.txt
    Sequence ids of 'na' are converted to None.
    Metadata field names are convert lower case with the white-spaces between words are converted to underscores.
    This allows both the older and new capitalization to work. Thus both "Assembly Name" and  "Assembly Name" will
    become "assembly_name".
    """

    expectedHeader = "# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name"

    def __init__(self, asmReport):
        self.metaData = MetaData()
        self.seqs = []
        self.bySequenceName = dict()
        self.byGenBankAccn = dict()
        self.byRefSeqAccn = dict()
        self.byUcscStyleName = dict()
        try:
            self._parse(asmReport)
        except Exception as ex:
            raise PycbioException(f"parse of NCBI assembly report '{asmReport}' failed") from ex

    def _parse(self, asmReport):
        with fileOps.opengz(asmReport) as fh:
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
            _parseMetaDataLine(self.metaData, line)

    def _skipToSeqTable(self, fh):
        "skip past header line before sequence records"
        for line in fh:
            if line[0:-1] == self.expectedHeader:
                return
        raise PycbioException("expected assembly report header not found in " + fh.name)

    def _parseRecord(self, line):
        rec = _parseRecord(line.split('\t'))
        self.seqs.append(rec)
        self.bySequenceName[rec.sequenceName] = rec
        if rec.genBankAccn is not None:
            self.byGenBankAccn[rec.genBankAccn] = rec
        if rec.refSeqAccn is not None:
            self.byRefSeqAccn[rec.refSeqAccn] = rec
        if rec.ucscStyleName is not None:
            self.byUcscStyleName[rec.ucscStyleName] = rec

    def _parseRecords(self, fh):
        for line in fh:
            self._parseRecord(line[0:-1])

    @property
    def assemblyName(self):
        return self.metaData.assembly_name

    def getBySequenceName(self, seqName):
        try:
            return self.bySequenceName[seqName]
        except KeyError:
            raise PycbioException(f"unknown '{self.assemblyName}' sequence name: '{seqName}'")

    def getByGenBankAccn(self, genBankAccn):
        try:
            return self.byGenBankAccn[genBankAccn]
        except KeyError:
            raise PycbioException(f"unknown '{self.assemblyName}'  a GenBank accession: '{genBankAccn}'")

    def getByRefSeqAccn(self, refSeqAccn):
        try:
            return self.byRefSeqAccn[refSeqAccn]
        except KeyError:
            raise PycbioException(f"unknown '{self.assemblyName}' RefSeq accession: '{refSeqAccn}'")

    def getByUcscStyleName(self, ucscStyleName):
        try:
            return self.byUcscStyleName[ucscStyleName]
        except KeyError:
            raise PycbioException(f"unknown '{self.assemblyName}' UCSC-style name: '{ucscStyleName}'")

    def getByName(self, name):
        """try all of the various sequence names to find a record"""
        rec = self.bySequenceName.get(name)
        if rec is None:
            rec = self.byGenBankAccn.get(name)
        if rec is None:
            rec = self.byRefSeqAccn.get(name)
        if rec is None:
            rec = self.byUcscStyleName.get(name)
        if rec is None:
            raise PycbioException(f"unknown '{self.assemblyName}' sequence: '{name}'")
        return rec

    def ucscNameToSeqName(self, ucscName):
        return self.getByUcscStyleName[ucscName].sequenceName

    def seqNameToUcscName(self, seqName):
        return self.getBySequenceName[seqName].ucscStyleName
