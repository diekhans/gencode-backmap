/*
 * Some basic type operations
 */
#ifndef typeOps_hh
#define typeOps_hh
#include <string>
#include <vector>
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

/*
 * Convert a string to an int.
 */
int stringToInt(const string& str,
                bool* isOk = NULL,
                int base = 10);

/* string vector */
typedef vector<const string> StringVector;

/*
 * Split a string into a vector of string given a separator character.
 */
StringVector stringSplit(const string& str,
                         char separator);
#endif
