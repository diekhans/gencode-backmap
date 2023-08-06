#include "FIOStream.hh"
#include "gzstream.hh"
#include <unistd.h>

// FIXME: drop file name of "-" convention

/* get file name based on "-" or "" being stdio */
static const string getRealFileName(const string& fileName,
                                    int ioMode) {
    if (fileName.empty() || (fileName == "-")) {
        return (ioMode == ios::in) ? "/dev/stdin" : "/dev/stdout";
    } else {
        return fileName;
    }
}

/**
 * Construct and open the stream.
 * @param fileName Name of file to open.  If `-' or empty, then
 *  stdin or stdout will be used.
 * @param ioMode ios::in or ios::out
 */
FIOStream::FIOStream(const string& fileName,
                     ios_base::openmode ioMode):
    iostream(openStreamBuf(getRealFileName(fileName, ioMode), isGZipFile(fileName), ioMode)),
    fFileName(getRealFileName(fileName, ioMode)) {
    // must call getRealFileName twice to deal with initialization order
}

/**
 * Open the file and allocate the streambuf.  Even though zlib can read
 * uncompressed files, we don't use it due to seek being slow.
 */
streambuf* FIOStream::openStreamBuf(const string& fileName,
                                    bool compressed,
                                    ios_base::openmode ioMode) {
    fGZFileBuf= NULL;
    fFileBuf = NULL;

    streambuf* strBuf = NULL;
    if (compressed) {
        // zlib based read or write
        fGZFileBuf = new gzstreambuf();
        strBuf = fGZFileBuf->open(fileName.c_str(), ioMode);
    } else {
        // FIXME: find a more portable way
        fFileBuf = new basic_filebuf<char>();
        strBuf = fFileBuf->open(fileName.c_str(), ioMode);
    }
    if (strBuf == NULL) {
        throw ios_base::failure(string("can't open \"") + fileName + "\" for " + ((ioMode & ios::out) ? "write" : "read") + " access");
    }
    return strBuf;
}

/**
 * Close the underlying file.
 */
void FIOStream::close() {
    if (fGZFileBuf != NULL) {
        fGZFileBuf->close();
    } else if (fFileBuf != NULL) {
        fFileBuf->close();
    }
    if (bad()) {
        throw ios_base::failure("I/O error on " + getFileName());
    }
}

/**
 * Destructor.
 */
FIOStream::~FIOStream() {
    close();
    if (fGZFileBuf != NULL) {
        delete fGZFileBuf;
    } else if (fFileBuf != NULL) {
        delete fFileBuf;
    }
}


/*
 * Local Variables:
 * mode: c++
 * End:
 */
