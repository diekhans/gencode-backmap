# Copyright 2006-2014 Mark Diekhans

# required enum34: https://pypi.python.org/pypi/enum34
from enum import Enum, EnumMeta

SymEnum = None  # class defined after use


class SymEnumValue(object):
    "Class used to define SymEnum member that have additional attributes."
    def __init__(self, value, externalName=None):
        self.value = value
        self.externalName = externalName


class _SysEnumExternalNameMap(object):
    "Mapping between internal and external member names"
    def __init__(self):
        self.intToExt = {}
        self.extToInt = {}

    def add(self, intName, extName):
        self.intToExt[intName] = extName
        self.extToInt[extName] = intName

    def toExtName(self, intName):
        "return name unchanged in no mapping"
        return self.intToExt.get(intName, intName)

    def toIntName(self, extName):
        "return name unchanged in no mapping"
        return self.extToInt.get(extName, extName)


class SymEnumMeta(EnumMeta):
    """metaclass for SysEnumMeta that implements looking up singleton members
    by string name."""

    @staticmethod
    def __symEnumValueUpdate(classdict, name, extNameMap):
        "record info about a member specified with SymEnum and update value in classdict"
        symValue = classdict[name]
        classdict[name] = symValue.value
        extNameMap.add(name, symValue.externalName)

    @staticmethod
    def __symEnumDerivedNew(metacls, cls, bases, classdict):
        "update class fields defined as SymEnumValue to register external names"
        extNameMap = classdict["__extNameMap__"] = _SysEnumExternalNameMap()
        for name in classdict.iterkeys():
            if isinstance(classdict[name], SymEnumValue):
                SymEnumMeta.__symEnumValueUpdate(classdict, name, extNameMap)
        return EnumMeta.__new__(metacls, cls, bases, classdict)

    def __new__(metacls, cls, bases, classdict):
        if SymEnum in bases:
            return SymEnumMeta.__symEnumDerivedNew(metacls, cls, bases, classdict)
        else:
            return EnumMeta.__new__(metacls, cls, bases, classdict)

    def __call__(cls, value, names=None, module=None, typ=None):
        "look up a value object, either by name of value,"
        if (names is None) and isinstance(value, str):
            # map string name to instance, check for external name
            value = cls.__extNameMap__.toIntName(value)
            member = cls._member_map_.get(value)
            if member is None:
                raise ValueError("'%s' is not a member or alias of %s" % (value, cls.__name__))
            else:
                return member
        else:
            return EnumMeta.__call__(cls, value, names, module, typ)


class SymEnum(Enum):
    """
    Metaclass for symbolic enumerations.  These are easily converted between
    string values and Enum objects.  This support construction from string
    values and str() returns value without class name.  Aliases can be
    added using the Enum approach of:
        val = 1
        valalias = val

    To handle string values that are not valid Python member names, an external
    name maybe associated with a field using a SymEnumValue object
        utr5 = SymEnumValue(1, "5'UTR")

    Either field name or external name maybe used to obtain a value.  The external
    name is returned with str().
    """
    __metaclass__ = SymEnumMeta

    def __str__(self):
        return self.__extNameMap__.toExtName(self.name)

    def __le__(self, other):
        if isinstance(other, SymEnum):
            return self.value <= other.value
        else:
            return self.value <= other

    def __lt__(self, other):
        if isinstance(other, SymEnum):
            return self.value < other.value
        else:
            return self.value < other

    def __ge__(self, other):
        if isinstance(other, SymEnum):
            return self.value >= other.value
        else:
            return self.value >= other

    def __gt__(self, other):
        if isinstance(other, SymEnum):
            return self.value > other.value
        else:
            return self.value > other

    def __eq__(self, other):
        if isinstance(other, SymEnum):
            return self.value == other.value
        else:
            return False

    def __ne__(self, other):
        if isinstance(other, SymEnum):
            return self.value != other.value
        else:
            return True
