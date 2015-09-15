# Copyright 2006-2012 Mark Diekhans
from pycbio.sys.typeOps import isListLike

class EnumValue(object):
    """A value of an enumeration.  The object id (address) is the
    unique value, with an associated display string and numeric value
    """
    __slots__ = ("enum", "name", "numValue", "strValue")
    def __init__(self, enum, name, numValue, strValue=None):
        self.enum = enum
        self.name = name
        self.numValue = numValue
        if strValue == None:
            self.strValue = name
        else:
            self.strValue = strValue

    def __getstate__(self):
        # optimize strValue if same as name
        return (self.enum, self.name, self.numValue, (None if self.strValue == self.name else self.strValue))

    def __setstate__(self, st):
        (self.enum,  self.name, self.numValue, self.strValue) = st
        if self.strValue == None:
            self.strValue = self.name

    def __rept__(self):
        return self.name

    def __str__(self):
        return self.strValue

    def __int__(self):
        return self.numValue

    def XXX__hash__(self):
        # FIXME: attempt to work around enum pickle problem in GenomeDefs.
        return hash(self.enum.name) + self.numValue

    def __cmp__(self, otherVal):
        # FIXME: attempt to work around enum pickle problem in GenomeDefs, compare names rather than
        # class objects.  Below should test be: not (isinstance(otherVal, EnumValue) and (otherVal.enum == self.enum)):
        if otherVal == None:
            return -1
        elif type(otherVal) == int:
            return cmp(self.numValue, otherVal)
        elif not isinstance(otherVal, EnumValue):
            raise TypeError("can't compare enumeration to type: " + str(type(otherVal)))
        elif otherVal.enum.name != self.enum.name:
            raise TypeError("can't compare enumerations of different types: "
                            + otherVal.enum.name + " and " + self.enum.name)
        else:
            return cmp(self.numValue, otherVal.numValue)

class Enumeration(object):
    """A class for creating enumeration objects.
    """

    def __init__(self, name, valueDefs, valueClass=EnumValue, bitSetValues=False):
        """Name is the name of the enumeration. ValueDefs is an ordered list of
        values.  If valueDefs contains a tuple, the first element is the
        value name, the second value is the __str__ value.  The third
        element is a list or tuple of string aliases that can be used to
        lookup the value under a different name, and the forth value
        is the numeric value.  Tuple elements maybe None or omitted.
        """
        self.name = name
        self.aliases = {}  # FIXME: not only aliases
        self.numValueMap = {}
        self.values = []
        self.maxNumValue = 0
        if bitSetValues:
            numValue = 1
        else:
            numValue = 0
        for valueDef in valueDefs:
            self.__defValue(valueClass, valueDef, numValue)
            if bitSetValues:
                numValue = numValue << 1
            else:
                numValue += 1
        self.values = tuple(self.values)

    def __len__(self):
        "return number of values"
        return len(self.values)

    def __createValue(self, valueClass, name, numValue, strValue):
        val = valueClass(self, name, numValue, strValue)
        setattr(self, name, val)
        self.aliases[name] = val
        if strValue != None:
            self.aliases[strValue] = val
        self.numValueMap[numValue] = val
        self.maxNumValue = max(self.maxNumValue, numValue)
        self.values.append(val)
        return val

    def __defValue(self, valueClass, valueDef, numValue):
        if isListLike(valueDef):
            return self.__defListValue(valueClass, valueDef, numValue)
        else:
            return self.__createValue(valueClass, valueDef, numValue, valueDef)

    def __defListValue(self, valueClass, valueDef, numValue):
        if (len(valueDef) > 3) and (valueDef[3] != None):
            numValue = valueDef[3]
            assert(isinstance(numValue, int))
        val = self.__createValue(valueClass, valueDef[0], numValue, valueDef[1])
        if (len(valueDef) > 2) and (valueDef[2] != None):
            if not isListLike(valueDef[2]):
                raise TypeError("valueDef[2] must be None, a list or tuple, found: " + str(valueDef[2]))
            for a in valueDef[2]:
                self.aliases[a] = val
        
    def lookup(self, name):
        """look up a value by name or aliases"""
        return self.aliases[name]

    def find(self, name):
        """find a value by name or aliases, or None if not found"""
        return self.aliases.get(name)

    def getValues(self, bitVals):
        "get a list of values associated with a bit set"
        vals = []
        for v in self.values:
            if v.numValue & bitVals:
                vals.append(v)
        return vals

    # FIXME: this doesn't belong (not an operation of self)
    def getValuesOr(self, vals):
        "get bit-wise or of the numeric values of a sequence of values"
        numVal = 0
        for v in vals:
            numVal |= int(v)
        return v

    def isValueOf(self, val):
        "is val an value of this enumeration?"
        for v in self.values:
            if val == v:
                return True
        return False

    def getByNumValue(self, numValue):
        "look up a value by numeric value"
        return self.numValueMap[numValue]

    # FIXME: emulates meta class new
    def __call__(self, name):
        return self.lookup(name)
