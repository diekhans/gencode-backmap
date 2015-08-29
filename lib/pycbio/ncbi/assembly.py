"""Parsing of NCBI assembly information files.
"""
from pycbio.sys import PycbioException
from pycbio.sys import fileOps

class AssemblyReport(object):
    """Parse assembly reports files, e.g.
    ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.28.assembly.txt
    """

    class Record(object):
        def __init__(self, sequenceName, sequenceRole, assignedMolecule, locationType, genBankAccn, relationship, refSeqAccn, assemblyUnit, sequenceLength, ucscStyleName):
            self.sequenceName = sequenceName 
            self.sequenceRole = sequenceRole 
            self.assignedMolecule = assignedMolecule 
            self.locationType = locationType 
            self.genBankAccn = genBankAccn 
            self.relationship = relationship 
            self.refSeqAccn = refSeqAccn 
            self.assemblyUnit = assemblyUnit 
            self.sequenceLength = sequenceLength 
            self.ucscStyleName = ucscStyleName 

        def __str__(self):
            return "\t".join([self.sequenceName, self.sequenceRole, self.assignedMolecule, self.locationType, self.genBankAccn, self.relationship, self.refSeqAccn, self.assemblyUnit, str(self.sequenceLength), self.ucscStyleName])
            
    expectedHeader = "# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name"
    def __init__(self, asmReport):
        self.seqs = []
        self.bySequenceName = dict()
        self.byGenBankAccn = dict()
        self.byRefSeqAccn = dict()
        self.byUcscStyleName = dict()
        with open(asmReport, "rU") as fh:
            self.__skipToHeader(fh)
            self.__parseRecords(fh)

    def __skipToHeader(self, fh):
        "skip past header line before sequence records"
        for line in fh:
            if line[0:-1] == self.expectedHeader:
                return
        raise PycbioException("expected assembly report header not found in " + fh.name)

    def __parseRecords(self, fh):
        for line in fh:
            self.__parseRecord(fh, line[0:-1])

    def __parseRecord(self, fh, line):
        row = line.split('\t')
        if len(row) != 10:
            raise PycbioException("expected 10 columns in assemble report record, found " + str(len(row)) + " in " + fh.name)
        rec = self.Record(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], int(row[8]), row[9])
        self.seqs.append(rec)
        self.bySequenceName[rec.sequenceName] = rec
        self.byGenBankAccn[rec.genBankAccn] = rec
        self.byRefSeqAccn[rec.refSeqAccn] = rec
        self.byUcscStyleName[rec.ucscStyleName] = rec
