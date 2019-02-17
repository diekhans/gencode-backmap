# Copyright 2006-2012 Mark Diekhans
"""Miscellaneous file operations"""

import six
import os
import sys
import errno
import re
import socket
import tempfile
import pipettor
from threading import Lock
from pycbio_local.sys import PycbioException

# FIXME: normalize file line read routines to all take fh or name, remove redundant code.


def ensureDir(dir):
    """Ensure that a directory exists, creating it (and parents) if needed."""
    # catching exception rather than checking for existence prevents race condition
    try:
        os.makedirs(dir)
    except OSError as ex:
        if ex.errno != errno.EEXIST:
            raise ex


def ensureFileDir(fname):
    """Ensure that the directory for a file exists, creating it (and parents) if needed.
    Returns the directory path"""
    dir = os.path.dirname(fname)
    if len(dir) > 0:
        ensureDir(dir)
        return dir
    else:
        return "."


def rmFiles(*fileArgs):
    """Remove one or more files if they exist. Each file argument can be a
    single file name of a list of file names"""
    for files in fileArgs:
        if isinstance(files, six.string_types):
            files = [files]
        for f in files:
            if os.path.exists(f):
                os.unlink(f)


def rmTree(root):
    "remove a file hierarchy, root can be a file or a directory"
    if os.path.isdir(root):
        for dir, subdirs, files in os.walk(root, topdown=False):
            for f in files:
                os.unlink(dir + "/" + f)
            os.rmdir(dir)
    else:
        if os.path.lexists(root):
            os.unlink(root)


def isCompressed(path):
    "determine if a file appears to be compressed by extension"
    return path.endswith(".gz") or path.endswith(".bz2") or path.endswith(".Z")


def compressCmd(path, default="cat"):
    """return the command to compress the path, or default if not compressed, which defaults
    to the `cat' command, so that it just gets written through"""
    if path.endswith(".Z"):
        raise PycbioException("writing compress .Z files not supported")
    elif path.endswith(".gz"):
        return "gzip"
    elif path.endswith(".bz2"):
        return "bzip2"
    else:
        return default


def compressBaseName(path):
    """if a file is compressed, return the path without the compressed extension"""
    if isCompressed(path):
        return os.path.splitext(path)[0]
    else:
        return path


def decompressCmd(path, default="cat"):
    """"return the command to decompress the file to stdout, or default if not compressed, which defaults
    to the `cat' command, so that it just gets written through"""
    if path.endswith(".Z") or path.endswith(".gz"):
        return "zcat"
    elif path.endswith(".bz2"):
        return "bzcat"
    else:
        return default


def opengz(fileName, mode="r", buffering=-1, encoding=None, errors=None):
    """open a file, if it ends in an extension indicating compression, open
    with a compression or decompression pipe."""
    if isCompressed(fileName):
        if mode.startswith("r"):
            cmd = decompressCmd(fileName)
            return pipettor.Popen([cmd, fileName], mode=mode, buffering=buffering, encoding=encoding, errors=errors)
        elif mode.startswith("w"):
            cmd = compressCmd(fileName)
            return pipettor.Popen([cmd], mode=mode, stdout=fileName, buffering=buffering, encoding=encoding, errors=errors)
        else:
            raise PycbioException("mode {} not support with compression for {}".format(mode, fileName))
    elif six.PY3:
        return open(fileName, mode, buffering=buffering, encoding=encoding, errors=errors)
    else:
        return open(fileName, mode)

# FIXME: make these consistent and remove redundant code.  Maybe use
# keyword for flush. Do we even need them with print function?


def prLine(fh, *objs):
    "write each str(obj) followed by a newline"
    for o in objs:
        fh.write(str(o))
    fh.write("\n")


def prsLine(fh, *objs):
    "write each str(obj), seperated by a space followed by a newline"
    n = 0
    for o in objs:
        if n > 0:
            fh.write(' ')
        fh.write(str(o))
        n += 1
    fh.write("\n")


def prOut(*objs):
    "write each str(obj) to stdout followed by a newline"
    for o in objs:
        sys.stdout.write(str(o))
    sys.stdout.write("\n")


def prErr(*objs):
    "write each str(obj) to stderr followed by a newline"
    for o in objs:
        sys.stderr.write(str(o))
    sys.stderr.write("\n")


def prsOut(*objs):
    "write each str(obj) to stdout, separating with spaces and followed by a newline"
    n = 0
    for o in objs:
        if n > 0:
            sys.stdout.write(' ')
        sys.stdout.write(str(o))
        n += 1
    sys.stdout.write("\n")


def prsfErr(*objs):
    "write each str(obj) to stderr, separating with spaces and followed by a newline and a flush"
    n = 0
    for o in objs:
        if n > 0:
            sys.stderr.write(' ')
        sys.stderr.write(str(o))
        n += 1
    sys.stderr.write("\n")
    sys.stderr.flush()


def prfErr(*objs):
    "write each str(obj) to stderr followed by a newline and a flush"
    for o in objs:
        sys.stderr.write(str(o))
    sys.stderr.write("\n")
    sys.stderr.flush()


def prsErr(*objs):
    "write each str(obj) to stderr, separating with spaces and followed by a newline"
    n = 0
    for o in objs:
        if n > 0:
            sys.stderr.write(' ')
        sys.stderr.write(str(o))
        n += 1
    sys.stderr.write("\n")


def prStrs(fh, *objs):
    "write each str(obj), with no newline"
    for o in objs:
        fh.write(str(o))


def prRow(fh, row):
    """Print a row (list or tupe) to a tab file.
    Does string conversion on each columns"""
    first = True
    for col in row:
        if not first:
            fh.write("\t")
        fh.write(str(col))
        first = False
    fh.write("\n")


def prRowv(fh, *objs):
    """Print a row from each argument to a tab file.
    Does string conversion on each columns"""
    first = True
    for col in objs:
        if not first:
            fh.write("\t")
        fh.write(str(col))
        first = False
    fh.write("\n")


class FileAccessor(object):
    """Context manager that opens a file (possibly compressed) if specified as
    a string, otherwise assume it is file-like and don't open/close"""
    def __init__(self, fspec, mode="r"):
        self.fspec = fspec
        self.mode = mode
        self.fh = None

    def __enter__(self):
        self.fh = opengz(self.fspec, self.mode) if isinstance(self.fspec, six.string_types) else self.fspec
        return self.fh

    def __exit__(self, typ, value, traceback):
        if isinstance(self.fspec, six.string_types):
            self.fh.close()


def iterLines(fspec):
    """generator over lines in file, dropping newlines.  If fspec is a string,
    open the file and close at end. Otherwise it is file-like object and will
    not be closed."""
    with FileAccessor(fspec) as fh:
        for line in fh:
            yield line[:-1]


def iterRows(fspec):
    """generator over rows in a tab-separated file.  Each line of the file is
    parsed, split into columns and returned.  If fspec is a string, open the
    file and close at end. Otherwise it is file-like object and will not be
    closed."""
    with FileAccessor(fspec) as fh:
        for line in fh:
            yield line[0:-1].split("\t")


def readFileLines(fspec):
    "read lines from a open file or a file name into a list, removing the newlines"
    return [l for l in iterLines(fspec)]


def readNonCommentLines(fspec):
    """read lines from an open file on file by name into a list, removing the
    newlines, striping leading and training white space, and skipping blank
    lines and those with the first non-space character is '#'."""
    lines = []
    for line in iterLines(fspec):
        line = line.strip()
        if (len(line) > 0) and (line[0] != '#'):
            lines.append(line)
    return lines


def readLine(fh):
    "read a line from a file, dropping a newline; None on eof"
    # FIXME: delete?
    line = fh.readline()
    if len(line) == 0:
        return None
    if line[-1:] == "\n":
        line = line[:-1]
    return line


def findTmpDir(tmpDir=None):
    """find the temporary directory to use, if tmpDir is not None, it is use"""
    if tmpDir is not None:
        return tmpDir
    tmpDir = os.getenv("TMPDIR")
    if tmpDir is not None:
        return tmpDir
    # UCSC special checks
    for tmpDir in ("/data/tmp", "/scratch/tmp", "/var/tmp", "/tmp"):
        if os.path.exists(tmpDir):
            return tmpDir
    raise PycbioException("can't find a tmp directory")


def setTmpEnv(tmpDir=None):
    """Setup TMPDIR env. If tmpDir arg is not None, set TMPDIR to this value.
    If tmpDir arg is not, keep TMPDIR if set, otherwise set to value returned
    by findTmpDir"""
    os.environ["TMPDIR"] = findTmpDir(tmpDir)


def tmpFileGet(prefix=None, suffix=".tmp", tmpDir=None):
    """Obtain a tmp file with a unique name in a secure way. File
    will only be accessible to user."""
    fh = tempfile.NamedTemporaryFile(prefix=prefix, suffix=suffix,
                                     dir=findTmpDir(tmpDir), delete=False)
    fh.close()
    return fh.name


def tmpDirGet(prefix=None, suffix=".tmp", tmpDir=None):
    """Obtain a tmp directory with a unique name in a secure way.  Directory
    will only be accessible to user."""
    return tempfile.mkdtemp(prefix=prefix, suffix=suffix, dir=findTmpDir(tmpDir))


_hostName = None  # get first time when needed
_atomicNextNum = 0  # number to include in atomicTmpFile, just in case same process tries creates multiple
_atomicTmpFileMutex = Lock()


def _tryForNewAtomicTmpFile(finalDir, finalBasename, finalExt):
    "Attempt to get a tmp file name, return None in the rare change that it exists"
    global _atomicNextNum
    baseName = "{}.{}.{}.{}.tmp{}".format(finalBasename, _hostName, os.getpid(), _atomicNextNum, finalExt)
    _atomicNextNum += 1
    tmpPath = os.path.join(finalDir, baseName)
    if os.path.exists(tmpPath):
        return None
    else:
        return tmpPath


def atomicTmpFile(finalPath):
    """Return a tmp file name to use with atomicInstall.  This will be in the
    same directory as finalPath. The temporary file will have the same extension
    as finalPath.  In final path is in /dev (/dev/null, /dev/stdout), it is
    returned unchanged and atomicTmpInstall will do nothing..  Thread-safe."""
    # FIXME: this would make a good object and maybe contact manager
    # note: this can't use tmpFileGet, since file should not be created or be private
    finalDir = os.path.dirname(os.path.normpath(finalPath))
    if finalDir == '/dev':
        return finalPath
    if finalDir == "":
        finalDir = '.'
    finalBasename = os.path.basename(finalPath)
    finalExt = os.path.splitext(finalPath)[1]
    global _hostName
    if _hostName is None:
        _hostName = socket.gethostname()
    with _atomicTmpFileMutex:
        tmpPath = None
        while tmpPath is None:
            tmpPath = _tryForNewAtomicTmpFile(finalDir, finalBasename, finalExt)
    return tmpPath


def atomicInstall(tmpPath, finalPath):
    "atomic install of tmpPath as finalPath"
    if os.path.dirname(os.path.normpath(finalPath)) != '/dev':
        os.rename(tmpPath, finalPath)


def uncompressedBase(path):
    "return the file path, removing a compression extension if it exists"
    if path.endswith(".gz") or path.endswith(".bz2") or path.endswith(".Z"):
        return os.path.splitext(path)[0]
    else:
        return path


_devNullFh = None


def getDevNull():
    "get a file object open to /dev/null, caching only one instance"
    global _devNullFh
    if _devNullFh is None:
        _devNullFh = open("/dev/null", "r+")
    return _devNullFh


def _parseMd5(line):
    "parse output of openssl md5"
    m = re.match("^MD5\\((.+)\\)= ([a-f0-9]+)\\n?", line)
    return m.group(1), m.group(2)


def md5sum(filePath):
    "compute md5 on a file"
    return _parseMd5(pipettor.runout(["openssl", "md5", filePath]))[1]


_sc_arg_max = None


def getArgMax():
    global _sc_arg_max
    if _sc_arg_max is None:
        _sc_arg_max = os.sysconf("SC_ARG_MAX")
    return _sc_arg_max


def _mkMd5SumCmd(filePaths, i, maxCmdLen):
    cmd = ["openssl", "md5"]
    cmdLen = 0
    while (cmdLen < maxCmdLen) and (i < len(filePaths)):
        cmd.append(filePaths[i])
        cmdLen += len(filePaths[i])
        i += 1
    return i, cmd


def _runMd5SumCmd(cmd):
    return [_parseMd5(line) for line in pipettor.runout(cmd)[0:-1].split('\n')]


def md5sums(filePaths):
    "compute md5 on a list of files, returning list of list of (path, sum)"
    maxCmdLen = getArgMax() - 1024  # a little padding
    i = 0
    results = []
    while i < len(filePaths):
        i, cmd = _mkMd5SumCmd(filePaths, i, maxCmdLen)
        results.extend(_runMd5SumCmd(cmd))
    return results
