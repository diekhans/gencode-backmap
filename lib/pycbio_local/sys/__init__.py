# Copyright 2006-2012 Mark Diekhans
from __future__ import print_function
import six
import sys
import inspect
import traceback

# FIXME: wanted to move this up a level, see doc/issues.org as why this didn't work


class PycbioException(Exception):
    """Base class for exceptions.  This implements exception chaining and
    stores a stack trace for python2 and the native exception chaining
    for python3.

    To chain an exception
       try:
          ...
       except Exception as ex:
          pycbioRaiseFrom(PycbioException("more stuff"(, ex)

    Optionally, traceback tb (sys.exc_info()[2] can be used to save get the stack
    trace.  Otherwise the current stack is used.

    Exceptions can be formatted using pycbioExFormat or
    pycbioExPrint.
    """

    def __init__(self, msg):
        super(PycbioException, self).__init__(msg)
        if six.PY2:
            self.__cause__ = None
            # Stack is saved as tuples to avoid traceback cycles
            self.__traceback_stack__ = _currentStack2(skip=0)


if six.PY2:  # noqa: C901
    def _frameToTuple2(frame):
        "convert a frame the traceback tuple format"
        finfo = inspect.getframeinfo(frame)
        if len(finfo.code_context) > 0:
            line = finfo.code_context[0].strip()
        else:
            line = ""
        return (finfo.filename, finfo.lineno, finfo.function, line)

    def _currentStack2(skip=0):
        frame = inspect.currentframe()
        stack = []
        while frame is not None:
            if skip <= 0:
                stack.insert(0, _frameToTuple2(frame))
            else:
                skip -= 1
            frame = frame.f_back
        return stack

    def _getTbStack2(ex):
        return getattr(ex, "__traceback_stack__", None)

    def _getTbStackSet2(ex):
        return frozenset(getattr(ex, "__traceback_stack__", ()))

    def _formatExInfo2(ex):
        return "{}: {}\n".format(type(ex).__name__, str(ex))

    def _formatExStack2(tbStack, higherTbEntries):
        results = ["Traceback (most recent call last):\n"]
        for tbEntry in tbStack:
            if tbEntry not in higherTbEntries:
                results.append(traceback.format_list([tbEntry])[0])
        return results

    def _formatOneEx2(ex, cause, higherTbEntries):
        results = []
        tbStack = _getTbStack2(ex)
        if tbStack is not None:
            results += _formatExStack2(tbStack, higherTbEntries)

        # if cause doe not have not a traceback, we need to include it now,
        if (cause is not None) and (_getTbStack2(cause) is None):
            results.append(_formatExInfo2(cause))
            results.append("Wrapped by:\n")
        results.append(_formatExInfo2(ex))
        return results

    def _dfsFormatStack2(ex, higherTbEntries=frozenset()):
        results = []
        cause = getattr(ex, "__cause__", None)
        if (cause is not None) and (_getTbStack2(cause) is not None):
            results += _dfsFormatStack2(cause, higherTbEntries | _getTbStackSet2(ex))
            results.append('\nThe above exception was the direct cause of the following exception:\n\n')
        results += _formatOneEx2(ex, cause, higherTbEntries)
        return results

    def _formatExceptChain2(ex):
        if _getTbStack2(ex) is not None:
            return _dfsFormatStack2(ex)
        else:
            # need some kind of stack, so use the current
            return _currentStack2(skip=0) + [_formatExInfo2(ex)]


def pycbioRaiseFrom(ex, cause=None):
    """Use this to raise with an cause or in cases where an exception maybe
    caught and added as a causes.  This will ensure capturing the stack."""
    if six.PY2:
        # save our stack information
        if isinstance(ex, PycbioException):
            ex.__cause__ = cause
            ex.__traceback_stack__ = _currentStack2(skip=0)
        else:
            # store traceback stack if we can, not all are writable
            try:
                setattr(ex, "__traceback_stack__", _currentStack2(skip=0))
            except (AttributeError, TypeError):
                pass
    six.raise_from(ex, cause)


def pycbioExFormat(ex):
    """Format any type of exception into list of new-line terminateds, handling
    PycbioException objects and saved stack traces"""
    if six.PY2:
        return _formatExceptChain2(ex)
    else:
        return traceback.format_exception(type(ex), ex, ex.__traceback__)


def pycbioExPrint(ex, file=None):
    """Format any type of exception into to a file, handling
    PycbioException objects and saved stack traces"""
    if six.PY2:
        print("".join(_formatExceptChain2(ex)), file=sys.stderr if file is None else file)
    else:
        return traceback.print_exception(type(ex), ex, ex.__traceback__, file=file)
