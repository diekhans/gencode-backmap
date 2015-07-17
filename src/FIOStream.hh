#ifndef SHARED_FIOSTREAM_H
#define SHARED_FIOSTREAM_H

#include <typeinfo>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
class gzstreambuf;
class filebuf;

/**
 * File iostream that supports gzipped files.  Automatically detects
 * gzipped input stream and decompresses.  On output, if the file
 * ends in .gz, the file is compressed.  Otherwise it is not compressed.
 */
class FIOStream: public iostream {
 private:
    /** Saved file name */
    string fFileName;

    /** zlib streambuf, if using zlib */
    gzstreambuf* fGZFileBuf;
    
    /** filebuf, if not using zlib */
    basic_filebuf<char>* fFileBuf;

    /** Internal methods */
    streambuf* openStreamBuf(const string& fileName,
                             bool compressed,
                             ios_base::openmode ioMode = ios::in);

    /**
     * Determine if a file name ends with .gz.
     */
    static bool isGZipFile(const string& fileName) {
        int len = fileName.size();
        return ((len > 3)
                && (fileName[len-3] == '.')
                && (fileName[len-2] == 'g')
                && (fileName[len-1] == 'z'));
    }

 public:    
    /**
     * Construct and open the stream.
     * @param fileName Name of file to open.  If `-' or empty, then
     *  stdin or stdout will be used.
     * @param ioMode ios::in or ios::out
     */
    FIOStream(const string& fileName,
              ios_base::openmode ioMode = ios::in);

    /** Close the underlying file. */
    void close();

    /** Destructor.*/
    virtual ~FIOStream();

    /** Get the file name. */
    const string& getFileName() const {
        return fFileName;
    }

    /** Determined if the file is compressed. */
    bool isCompressed() const {
        return (fGZFileBuf != NULL);
    }

    /** read a line, return false if on EOF */
    bool readLine(string& line) {
        std::getline(*this, line);
        if (eof()) {
            return false;  //EOF
        }
        if (fail()) {
            throw ios_base::failure("I/O error on " + getFileName());
        }
        return true;
    }
};

#endif

/*
 * Local Variables:
 * mode: c++
 * End:
 */
