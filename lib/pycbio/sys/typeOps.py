# Copyright 2006-2012 Mark Diekhans
"""Miscellaneous type operations"""
# FIXME: move to other modules or move set in here.
from builtins import range
from pycbio.sys import PycbioException


def isListLike(v):
    "is variable a list or tuple?"
    return isinstance(v, list) or isinstance(v, tuple)


def listInit(size, val):
    "create a list of length size, with each element containing val"
    lst = []
    for i in range(size):
        lst.append(val)
    return lst


def listAppend(lst, item):
    """if lst is None, create a new list with item, otherwise append item.
    Returns list"""
    if lst is None:
        return [item]
    else:
        lst.append(item)
    return lst


def listExtend(lst, items):
    """if lst is None, create a new list with items, otherwise extend with items.
    Returns list"""
    if lst is None:
        return list(items)
    else:
        lst.extend(items)
    return lst


def isIterable(v):
    "is variable a list, tuple, set, or hash? str doesn't count"
    # FIXME: bad name, as strings are iterable
    return isinstance(v, list) or isinstance(v, tuple) or isinstance(v, set) or isinstance(v, dict)


def mkiter(item):
    """create a iterator over item, if item is iterable, just return an iter,
    if item is not iterable or is a string, create an iterable to return just
    item, if item is none, return an empty iter"""
    # FIXME: don't really need to construct a list
    if item is None:
        return iter(())
    elif isIterable(item):
        return iter(item)
    else:
        return iter([item])


def mkset(item):
    """create a set from item.  If it's None, return an empty set, if it's
    iterable, convert to a set, if it's a single item, make a set of it,
    it it's already a set, just return as-is"""
    # FIXME: move to setOps
    if isinstance(item, set):
        return item
    elif item is None:
        return set()
    elif isIterable(item):
        return set(item)
    else:
        return set([item])


def noneOrZero(v):
    "test if a value is either None or len of zero"
    return (v is None) or (len(v) == 0)


def addUniq(d, k, v):
    "add to a dict, generating an error if the item already exists"
    if k in d:
        raise PycbioException("item \"{}\" already in dict".format(str(k)))
    d[k] = v


def dictObtain(d, key, mkFunc):
    "return entry d[key], creating with mkFunc if it doesn't exist"
    if key not in d:
        e = d[key] = mkFunc()
    else:
        e = d[key]
    return e


def _annonStr(self):
    "__str__ function for annon objects"
    return ", ".join(["{}={}".format(k, repr(getattr(self, k))) for k in sorted(dir(self)) if not k.startswith("__")])


def annon(**kwargs):
    """create an anonymous object with fields that are same as the keyword
    arguments.  This is different from a named tuple as it create an object without
    defining a class.
    """
    kwargs.update({"__str__": _annonStr})
    return type('', (), kwargs)()


def attrdict(obj):
    """Create a dictionary of all attributes in an object. This will work for
    classes with __slots__ or __dict__.    The returned object may or may not be the
    object and so should *not be modified*"""
    try:
        return vars(obj)
    except TypeError:
        return {a: getattr(obj, a) for a in dir(obj)}


__all__ = (isListLike.__name__, listInit.__name__, listAppend.__name__,
           isIterable.__name__, mkiter.__name__, mkset.__name__,
           noneOrZero.__name__, addUniq.__name__, dictObtain.__name__,
           annon.__name__, attrdict.__name__)
