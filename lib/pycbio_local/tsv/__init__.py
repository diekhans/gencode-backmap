# Copyright 2006-2012 Mark Diekhans
"""" TSV (Tab Separated File) parsing"""

from pycbio_local.sys import PycbioException


class TsvError(PycbioException):
    "Error from reading or parsing a TSV file"
    def __init__(self, msg, reader=None):
        if (reader is not None):
            msg = str(reader.fileName) + ":" + str(reader.lineNum) + ": " + msg
        super(TsvError, self).__init__(msg)


from pycbio_local.tsv.tsvRow import TsvRow, tsvRowToDict
from pycbio_local.tsv.tsvReader import TsvReader, strOrNoneType, intOrNoneType, printf_basic_dialect

__all__ = (TsvError.__name__, TsvRow.__name__, TsvReader.__name__,
           "strOrNoneType", "intOrNoneType",
           tsvRowToDict.__name__,
           printf_basic_dialect.__name__)
