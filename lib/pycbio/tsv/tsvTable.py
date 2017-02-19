# Copyright 2006-2012 Mark Diekhans
from pycbio.tsv.tsvReader import TsvReader
from pycbio.tsv import TsvError
from pycbio.sys.multiDict import MultiDict
import sys
import csv

# FIX: maybe make each index it's own class to handle uniq check, etc.


class TsvTable(list):
    """Class for reading and writing TSV files. Stores rows as a list of Row
    objects. Columns are indexed by name.

    - idx - index objects dict or MultiDict attribute per keyed columns.
    """
    class Indices (object):
        """object with attribute for each key column"""

        def __getitem__(self, key):
            return getattr(self, key)

    def __addIndex(self, keyCol, dictClass):
        if keyCol not in self.colMap:
            raise TsvError("key column \"" + keyCol + "\" is not defined"), None, sys.exc_info()[2]
        setattr(self.idx, keyCol, dictClass())

    def __createIndices(self, keyCols, dictClass):
        "keyCols maybe string or seq of strings"
        if type(keyCols) == str:
            self.__addIndex(keyCols, dictClass)
        else:
            for kc in keyCols:
                self.__addIndex(kc, dictClass)

    def __buildIndices(self, uniqKeyCols, multiKeyCols):
        self.idx = TsvTable.Indices()
        self.indices = self.idx  # FIXME: old name, delete
        if uniqKeyCols is not None:
            self.__createIndices(uniqKeyCols, dict)
        if multiKeyCols is not None:
            self.__createIndices(multiKeyCols, MultiDict)

    def __buildColDictTbl(self):
        """build an array, index by column number, of dict objects, or None if not
        indexed.  Used when loading rows. """
        if len(vars(self.idx)) == 0:
            return None
        tbl = []
        for iCol in xrange(len(self.columns)):
            tbl.append(getattr(self.idx, self.columns[iCol], None))
        return tbl

    def __indexCol(self, iCol, colDict, col, row):
        if (type(colDict) == dict) and colDict.get(col):
            raise Exception("column {} unique index value already entered: {} from {} ".format(self.columns[iCol], col, row))
        else:
            colDict[col] = row

    def __indexRow(self, colDictTbl, row):
        for i in xrange(len(row)):
            if colDictTbl[i] is not None:
                self.__indexCol(i, colDictTbl[i], row[i], row)

    # FIXME: need add row function, but colDict stuff conflicts, make member
    def __readBody(self, reader):
        colDictTbl = self.__buildColDictTbl()
        for row in reader:
            self.append(row)
            if colDictTbl is not None:
                self.__indexRow(colDictTbl, row)

    def __init__(self, fileName, uniqKeyCols=None, multiKeyCols=None, rowClass=None, typeMap=None,
                 defaultColType=None, columns=None, columnNameMapper=None, ignoreExtraCols=False, isRdb=False, inFh=None, allowEmpty=False, dialect=csv.excel_tab):
        """Read TSV file into the object

        fileName - name of file, opened unless inFh is specified
        uniqKeyCols - name or names of columns to index with uniq keys,
            can be string or sequence
        multiKeyCols - name or names of columns to index, allowing multiple keys.
            can be string or sequence
        rowClass - class or class factory function to use for a row. Must take
            TsvReader and list of string values of columns.
        typeMap - if specified, it maps column names to the type objects to
            use to convert the column.  Unspecified columns will not be
            converted. Key is the column name, value can be either a type
            or a tuple of (parseFunc, formatFunc).  If a type is use,
            str() is used to convert to a printable value.
        defaultColType - if specified, type of unspecified columns
        columns - if specified, the column names to use.  The header
            should not be in the file.
        ignoreExtraCols - should extra columns be ignored?
        isRdb - file is an RDB file, ignore second row (type map still needed).
        inFh - If not None, this is used as the open file, rather than
          opening it.  Closed when the end of file is reached.
        allowEmpty - an empty input results in an EOF rather than an error.
          Should specify this if reading from a database query.
        """
        reader = TsvReader(fileName, rowClass=rowClass, typeMap=typeMap, defaultColType=defaultColType, isRdb=isRdb, columns=columns, columnNameMapper=columnNameMapper,
                           ignoreExtraCols=ignoreExtraCols, inFh=inFh, allowEmpty=allowEmpty, dialect=dialect)
        try:
            self.columns = reader.columns
            self.colTypes = reader.colTypes
            self.colMap = reader.colMap
            self.__buildIndices(uniqKeyCols, multiKeyCols)
            self.__readBody(reader)
        except Exception as ex:
            raise TsvError("load failed", reader=reader, cause=ex), None, sys.exc_info()[2]

    def addColumn(self, colName, initValue=None, colType=None):
        "add a column to all rows in the table"
        if colName in self.colMap:
            raise TsvError("column \"" + colName + "\" is already defined"), None, sys.exc_info()[2]

        self.colMap[colName] = len(self.columns)
        if colType:
            assert(self.colTypes)
            self.colTypes.append(colType)
        elif self.colTypes:
            self.colTypes.append(None)
        self.columns.append(colName)

        # add column to each row
        for row in self:
            setattr(row, colName, initValue)

    def write(self, fh):
        fh.write(str.join("\t", self.columns))
        fh.write("\n")
        for row in self:
            row.write(fh)

# FIXME: duped in TabFile, also file ops


def tsvPrRow(fh, row):
    """Print a row (list or tupe) to a tab file.
    does string conversions on columns"""
    cnt = 0
    for col in row:
        if cnt > 0:
            fh.write("\t")
        fh.write(str(col))
        cnt += 1
    fh.write("\n")
