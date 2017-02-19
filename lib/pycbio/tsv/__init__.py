# Copyright 2006-2012 Mark Diekhans
"""" TSV (Tab Separated File) parsing"""

from pycbio.sys import PycbioException


class TsvError(PycbioException):
    "Error from reading or parsing a TSV file"
    def __init__(self, msg, reader=None, cause=None):
        if (reader is not None):
            msg = str(reader.fileName) + ":" + str(reader.lineNum) + ": " + msg
        PycbioException.__init__(self, msg, cause)


from pycbio.tsv.tsvRow import TsvRow
from pycbio.tsv.tsvReader import TsvReader, strOrNoneType, intOrNoneType
from pycbio.tsv.tsvTable import TsvTable
from pycbio.tsv.tabFile import TabFile
from pycbio.tsv.tabFile import TabFileReader

__all__ = (TsvError.__name__, TsvRow.__name__, TsvReader.__name__, "strOrNoneType", "intOrNoneType", TsvTable.__name__, TabFile.__name__, TabFileReader.__name__)
