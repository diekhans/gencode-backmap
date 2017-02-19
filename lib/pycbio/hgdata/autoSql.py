# Copyright 2006-2012 Mark Diekhans
"""support classes for parsing autoSql generated objects"""
import string


def strArraySplit(commaStr):
    "parser for comma-separated string list into a list"
    if len(commaStr) == 0:
        return []
    strs = commaStr.split(",")
    if commaStr.endswith(","):
        strs = strs[0:-1]
    return strs


def strArrayJoin(strs):
    "formatter for a list into a comma seperated string"
    if strs is not None:
        return string.join(strs, ",") + ","
    else:
        return ","

# TSV typeMap tuple for str arrays
strArrayType = (strArraySplit, strArrayJoin)


def intArraySplit(commaStr):
    "parser for comma-separated string list into a list of ints"
    ints = []
    for s in strArraySplit(commaStr):
        ints.append(int(s))
    return ints


def intArrayJoin(ints):
    "formatter for a list of ints into a comma seperated string"
    if ints is not None:
        strs = []
        for i in ints:
            strs.append(str(i))
        return string.join(strs, ",") + ","
    else:
        return ","

# TSV typeMap tuple for str arrays
intArrayType = (intArraySplit, intArrayJoin)
