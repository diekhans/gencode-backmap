#include "typeOps.hh"
#include <stdexcept>
#include <errno.h>
#include <stdlib.h>

/*
 * Convert a string to an int.
 */
int stringToInt(const string& str,
                bool* isOk,
                int base) {
    const char* cstr = str.c_str();
    char *endPtr;
    errno = 0;
    long lnum = strtol(cstr, &endPtr, base);
    if ((endPtr == cstr) || (*endPtr != '\0')) {
        if (isOk != NULL) {
            *isOk = false;
            return 0;
        } else {
            throw invalid_argument("Invalid integer \"" + string(str) + "\"");
        }
    }
     
    int num = (int)lnum;
    if ((errno != 0) || ((long)num != lnum)) {
        if (isOk != NULL) {
            *isOk = false;
            return 0;
        } else {
            throw invalid_argument("Integer out of range \"" + string(str) + "\"");
        }
    }
    if (isOk != NULL) {
        *isOk = true;
    }
    return num;
}

/*
 * Split a string into a vector of string given a separator character.
 */
StringVector stringSplit(const string& str,
                         char separator) {
    StringVector strs;
    
    int prevIdx = 0;
    int sepIdx;
    while ((sepIdx = (int)str.find_first_of(separator, prevIdx)) >= 0) {
        strs.push_back(str.substr(prevIdx, sepIdx-prevIdx));
        prevIdx = sepIdx+1;
    }
    strs.push_back(str.substr(prevIdx));

    return strs;
}

