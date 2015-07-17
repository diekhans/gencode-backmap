/*
 * Some basic type operations
 */
#ifndef typeOps_hh
#define typeOps_hh
#include <string>
#include <vector>
using namespace std;

/*
 * Convert a string to an int.
 */
int stringToInt(const string& str,
                bool* isOk = NULL,
                int base = 10);

/*
 * Split a string into a vector of string given a separator character.
 */
vector<const string> stringSplit(const string& str,
                                 char separator);
#endif
