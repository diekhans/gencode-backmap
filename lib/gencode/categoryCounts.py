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
    def __init__(self, keyColumnHeader):
        self.keyColumnHeader = tuple(keyColumnHeader)
        self.counts = defaultdict(int)
        self.total = 0

    def count(self, key):
        self.counts[key] += 1
        self.total += 1

    def __keyFromRow(self, keyCols, row):
        return tuple([row[k] for k in keyCols])
        
    def countTsv(self, tsvFile, keyCols=None, getKeys=None, typeMap=None, filterFunc=None):
        "allows computed keys"
        assert ((keyCols == None) and (getKeys != None)) or ((keyCols != None) and (getKeys == None))
        if keyCols != None:
            getKeys = lambda row: self.__keyFromRow(keyCols, row)
        for row in TsvReader(tsvFile, typeMap=typeMap):
            if (filterFunc == None) or filterFunc(row):
                self.count(getKeys(row))

    def write(self, fh, inclTotals=False):
        fileOps.prRow(fh, self.keyColumnHeader + ("count", "freq"))
        for key in sorted(self.counts.iterkeys()):
            val = self.counts[key]
            fileOps.prRow(fh, key + (val, fmtRate(val, self.total)))
        if inclTotals:
            fileOps.prRow(fh, ["all" for col in self.keyColumnHeader]
                          + [self.total, fmtRate(self.total, self.total)])
                
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
