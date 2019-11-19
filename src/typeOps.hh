/*
 * Some basic type operations
 */
#ifndef typeOps_hh
#define typeOps_hh
#include <string>
#include <vector>
#include <set>
#include "jkinclude.hh"
using namespace std;

#include <typeinfo>

/**
 * @name instanceOf
 *
 * Check at run time if an object is an instance of a class.  This differs
 * from comparing the result of type_id, as it takes inheritance into account.
 *
 * @param objPtr A pointer to the object to check.  A pointer must be used, not
 *  an object.  If necessary, generate an address.
 * @param The class name to test against.
 * @return true if the object is an instance of the class.
 */
#define instanceOf(objPtr, className) \
    ((bool)(dynamic_cast<const className*>(objPtr) != NULL))


/* vector of integers */
typedef vector<int> IntVector;

/* a pair of strings */
typedef pair<string, string> StringPair;

/* set of strings */
typedef set<string> StringSet;

/* whitespace characters */
extern const string whitespace;

/* convert a string to a non-const char* for calling kent library */
inline char *toCharStr(const string& s) {
    return const_cast<char*>(s.c_str());
}

/* does a string start with a prefix? */
inline bool stringStartsWith(const string& s, const string& prefix) {
    return (s.size() >= prefix.size()) && std::equal(prefix.begin(), prefix.end(), s.begin());
}

/* does a string end with a suffix? */
inline bool stringEndsWith(const string& s, const string& suffix) {
    return (s.size() >= suffix.size()) && (s.rfind(suffix) == (s.size()-suffix.size()));
}

/*
 * Convert a string to an int.
 */
int stringToInt(const string& str,
                bool* isOk = NULL,
                int base = 10);

/* an empty string constant */
extern const string emptyString;

/* string vector */
typedef vector<string> StringVector;

/*
 * Split a string into a vector of string given a separator character.
 */
StringVector stringSplit(const string& str,
                         char separator);

/*
 * Join a string into a vector into a string
 */
string stringJoin(const StringVector& strvec,
                  char separator);

/* is a string empty (only whitespace) */
inline bool stringEmpty(const string& s) {
    return s.find_first_not_of(whitespace) == string::npos;
}

/* trim from end of string */
inline string stringRtrim(const string& s, const string& t = whitespace) {
    string ls = s;
    ls.erase(ls.find_last_not_of(t) + 1);
    return ls;
}

/* trim from beginning of string */
inline string stringLtrim(const string& s, const string& t = whitespace) {
    string ls = s;
    ls.erase(0, ls.find_first_not_of(t));
    return ls;
}

/* trim from both ends of string */
inline string stringTrim(const string& s, const string& t = whitespace) {
    return stringLtrim(stringRtrim(s, t), t);
}

/** Convert an integer to a string. */
string toString(int num);

/** Convert an character to a string. */
inline string charToString(char ch) {
    return string(1, ch);
}

#endif
