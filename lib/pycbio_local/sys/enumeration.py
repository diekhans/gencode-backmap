# Copyright 2006-2012 Mark Diekhans
from past.builtins import cmp
from builtins import object
import six
from pycbio_local.sys.typeOps import isListLike

# FIXME:
# - could have user add value object direcrly instead of complex value tuples
# - should be iterable
# - need to be able to pickle objects contains EnumValue and then
#   compare then with code constants after load
# - should really be a tuple instead of values field
# FIXME: should really use meta classes; some enumeration stuff:
# FIXME: should be iteratable
#
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/67107
# http://pytable.sourceforge.net/pydoc/basictypes.enumeration.html
# http://www.bgb.cc/garrett/originals/
# http://www.python.org/doc/essays/metaclasses/Enum.py
#   http://svn.python.org/projects/python/trunk/Demo/newmetaclasses/Enum.py
# http://www.python.org/cgi-bin/moinmoin/EnumerationProgramming
# http://techspot.zzzeek.org/2011/01/14/the-enum-recipe/
#
# This looks really good.
# http://stackoverflow.com/questions/100003/what-is-a-metaclass-in-python

# FIXME FIXME FIXME: immutable stuff and solts disabled due to problems with
# serialization.
# FIXME: should be able to specify a EnumValue; although need a way to add
# aliases in this case.


class EnumValue(object):
    """A value of an enumeration.  The object id (address) is the
    unique value, with an associated display string and numeric value
    """
    __slots__ = ("enum", "name", "numValue", "strValue")

    def __init__(self, enum, name, numValue, strValue=None):
        self.enum = enum
        self.name = name
        self.numValue = numValue
        if strValue is None:
            self.strValue = name
        else:
            self.strValue = strValue

    def __getstate__(self):
        # optimize strValue if same as name
        return (self.enum, self.name, self.numValue, (None if self.strValue == self.name else self.strValue))

    def __setstate__(self, st):
        (self.enum, self.name, self.numValue, self.strValue) = st
        if self.strValue is None:
            self.strValue = self.name

    def __rept__(self):
        return self.name

    def __str__(self):
        return self.strValue

    def __int__(self):
        return self.numValue

    def __hash__(self):
        return hash(self.numValue)

    if six.PY2:
        def __cmp__(self, otherVal):
            # FIXME: attempt to work around enum pickle problem in GenomeDefs, compare names rather than
            # class objects.  Below should test be: not (isinstance(otherVal, EnumValue) and (otherVal.enum == self.enum)):
            if otherVal is None:
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

    def __le__(self, other):
        if isinstance(other, EnumValue):
            return self.numValue <= other.numValue
        else:
            return self.numValue <= other

    def __lt__(self, other):
        if isinstance(other, EnumValue):
            return self.numValue < other.numValue
        else:
            return self.numValue < other

    def __ge__(self, other):
        if isinstance(other, EnumValue):
            return self.numValue >= other.numValue
        else:
            return self.numValue >= other

    def __gt__(self, other):
        if isinstance(other, EnumValue):
            return self.numValue > other.numValue
        else:
            return self.numValue > other

    def __eq__(self, other):
        if isinstance(other, EnumValue):
            return self.numValue == other.numValue
        else:
            return self.numValue == other

    def __ne__(self, other):
        if isinstance(other, EnumValue):
            return self.numValue != other.numValue
        else:
            return self.numValue != other


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
            self._defValue(valueClass, valueDef, numValue)
            if bitSetValues:
                numValue = numValue << 1
            else:
                numValue += 1
        self.values = tuple(self.values)

    def __len__(self):
        "return number of values"
        return len(self.values)

    def _createValue(self, valueClass, name, numValue, strValue):
        val = valueClass(self, name, numValue, strValue)
        setattr(self, name, val)
        self.aliases[name] = val
        if strValue is not None:
            self.aliases[strValue] = val
        self.numValueMap[numValue] = val
        self.maxNumValue = max(self.maxNumValue, numValue)
        self.values.append(val)
        return val

    def _defValue(self, valueClass, valueDef, numValue):
        if isListLike(valueDef):
            return self._defListValue(valueClass, valueDef, numValue)
        else:
            return self._createValue(valueClass, valueDef, numValue, valueDef)

    def _defListValue(self, valueClass, valueDef, numValue):
        if (len(valueDef) > 3) and (valueDef[3] is not None):
            numValue = valueDef[3]
            assert(isinstance(numValue, int))
        val = self._createValue(valueClass, valueDef[0], numValue, valueDef[1])
        if (len(valueDef) > 2) and (valueDef[2] is not None):
            if not isListLike(valueDef[2]):
                raise TypeError("valueDef[2] must be None, a list or tuple, found: " + str(valueDef[2]))
            for a in valueDef[2]:
                self.aliases[a] = val

    def X__getstate__(self):
        return (self.name, self.aliases, self.values, self.maxNumValue)

    def X__setstate__(self, st):
        (self.name, self.aliases, self.values, self.maxNumValue) = st
        for val in self.values:
            setattr(self, val.name, val)

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
