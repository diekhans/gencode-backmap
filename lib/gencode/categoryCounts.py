from collections import defaultdict
from pycbio.sys import fileOps
from pycbio.tsv import TsvReader

def fmtRate(count, total):
    if total == 0:
        return "0.0"
    else:
        return format(float(count)/total, "0.2f")


class CategoryCounts(object):
    """collect counts by a category expressed as a tuple.
    """

    defaultColumnKey = "count"

    def __init__(self, rowKeyHeader):
        self.rowKeyHeader = tuple(rowKeyHeader)
        # Each count is a row, keyed by the matrix column key.  When
        # this is not specified, there is a single column
        self.counts = defaultdict(lambda : defaultdict(int))
        self.rowTotals = defaultdict(int)
        self.columnTotals = defaultdict(int)
        self.allRowKeys = set()
        self.allColumnKeys = set()

    def count(self, rowKey, columnKey):
        self.counts[rowKey][columnKey] += 1
        self.rowTotals[rowKey] += 1
        self.columnTotals[columnKey] += 1
        self.allRowKeys.add(rowKey)
        self.allColumnKeys.add(columnKey)

    def countTsv(self, tsvFile, getRowKeys, typeMap=None, filterFunc=None, getColumnKey=None):
        "allows computed keys"
        for row in TsvReader(tsvFile, typeMap=typeMap):
            if (filterFunc == None) or filterFunc(row):
                self.count(getRowKeys(row),
                           (self.defaultColumnKey if getColumnKey == None else getColumnKey(row)))

    def writeCountsFreqs(self, fh, inclTotals=False):
        total = self.columnTotals[self.defaultColumnKey]
        fileOps.prRow(fh, self.rowKeyHeader + ("count", "freq"))
        for key in sorted(self.counts.iterkeys()):
            val = self.counts[key][self.defaultColumnKey]
            fileOps.prRow(fh, key + (val, fmtRate(val, total)))
        if inclTotals:
            fileOps.prRow(fh, ["all" for col in self.rowKeyHeader] + [total, fmtRate(total, total)])
                
    def __writeMatrixRow(self, fh, rowKey, rowData, columnKeys, rowFrequencies, columnFrequencies, inclTotals):
        if rowFrequencies:
            rowTotal = self.rowTotals[rowKey]
            data = [fmtRate(rowData[ck], rowTotal) for ck in columnKeys]
            if inclTotals:
                data.append(fmtRate(sum([rowData[ck] for ck in columnKeys]), rowTotal))
        elif columnFrequencies:
            data = [fmtRate(rowData[ck], self.columnTotals[ck]) for ck in columnKeys]
        else:
            data = [rowData[ck] for ck in columnKeys]
            if inclTotals:
                data.append(sum(data))
        fileOps.prRow(fh, list(rowKey) + data)

    def writeMatrix(self, fh, rowFrequencies=False, columnFrequencies=False, inclTotals=False):
        columnKeys = tuple(sorted(self.allColumnKeys))
        fileOps.prRow(fh, self.rowKeyHeader + columnKeys + (("total",) if inclTotals else ()))
        for rowKey in sorted(self.counts.iterkeys()):
            self.__writeMatrixRow(fh, rowKey, self.counts[rowKey], columnKeys, rowFrequencies, columnFrequencies, inclTotals)
        if inclTotals and not rowFrequencies:
            rowKey = ["all" for col in self.rowKeyHeader]
            self.__writeMatrixRow(fh, rowKey, self.columnTotals, columnKeys, rowFrequencies, columnFrequencies, inclTotals)
                
class CategoryStatsCombine(object):
    "combine statistic from multiple CategoryCounts TSVs"
    def __init__(self):
        self.inKeyColNames = None
        self.totalsKey = None
        self.allLabels = []  # column pairs
        self.allKeys = set()
        self.stats = defaultdict(defaultdict) # by key (row) then by label (column)

    def __addNewLabel(self, label):
        if label in self.allLabels:
            raise Exception("duplicate column label: " + label)
        self.allLabels.append(label)
            
    def __newSetInit(self, row):
        keyColNames = tuple(row._columns_[0:-2])
        if self.inKeyColNames == None:
            self.inKeyColNames = keyColNames
            self.totalsKey = tuple(["all" for col in self.inKeyColNames])
        elif keyColNames != self.inKeyColNames:
            raise Exception("expected columns of: (" + ", ".join(self.inKeyColNames) + ") got (" + ", ".join(keyColNames) + ")")

    def __addRow(self, label, row):
        key = tuple(row.getRow()[0:-2])
        self.allKeys.add(key)
        if label in self.stats[key]:
            raise Exception("duplicate row key for " + label + ": " + str(row))
        self.stats[key][label] = (row[-2], row[-1])
        
    def add(self, label, tsvFile):
        self.__addNewLabel(label)
        iRow = 0
        for row in TsvReader(tsvFile):
            if iRow == 0:
                self.__newSetInit(row)
            self.__addRow(label, row)
            iRow += 1
        
    def __writeHeader(self, fh):
        row = list(self.inKeyColNames)
        for label in self.allLabels:
            row.extend([label+"_count", label+"_freq"])
        fileOps.prRow(fh, row)

    __emptyCounts = ("0", "0.0")
    def __writeRow(self, fh, key):
        row = list(key)
        for label in self.allLabels:
            row.extend(self.stats[key].get(label, self.__emptyCounts))
        fileOps.prRow(fh, row)

    def write(self, fh):
        self.__writeHeader(fh)
        for key in sorted(self.allKeys):
            if key != self.totalsKey:
                self.__writeRow(fh, key)
        if self.totalsKey in self.stats:
            self.__writeRow(fh, self.totalsKey)
