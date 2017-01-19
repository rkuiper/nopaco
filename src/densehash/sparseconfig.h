/* src/densehash/sparseconfig.h.  Generated from sparseconfig.h.in by configure.  */
/* src/sparsehash/internal/sparseconfig.h.in.  Generated from configure.ac by autoheader.  */

/* Namespace for Google classes */
#define GOOGLE_NAMESPACE ::nopacoHASH

/* the location of the header defining hash functions */
#define HASH_FUN_H <functional>

/* the location of <unordered_map> or <hash_map> */
#define HASH_MAP_H <unordered_map>

/* the namespace of the hash<> function */
#define HASH_NAMESPACE std

/* the location of <unordered_set> or <hash_set> */
#define HASH_SET_H <unordered_set>

/* define if the compiler has hash_map */
#define HAVE_HASH_MAP 1

/* define if the compiler has hash_set */
#define HAVE_HASH_SET 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if the system has the type `long long'. */
#define HAVE_LONG_LONG 1

/* Define to 1 if you have the `memcpy' function. */
#define HAVE_MEMCPY 1

/* Define to 1 if you have the `memmove' function. */
#define HAVE_MEMMOVE 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if the system has the type `uint16_t'. */
#define HAVE_UINT16_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* define if the compiler supports unordered_{map,set} */
#define HAVE_UNORDERED_MAP 1

/* Define to 1 if the system has the type `u_int16_t'. */
#define HAVE_U_INT16_T 1

/* Define to 1 if the system has the type `__uint16'. */
/* #undef HAVE___UINT16 */

/* Name of package */
/* #undef PACKAGE */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "R.Kuiper"

/* Define to the full name of this package. */
#define PACKAGE_NAME "nopaco"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "nopaco 0.99.5"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "nopaco"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.99.5"

/* The system-provided hash function including the namespace. */
#define SPARSEHASH_HASH HASH_NAMESPACE::hash

/* The system-provided hash function, in namespace HASH_NAMESPACE. */
#define SPARSEHASH_HASH_NO_NAMESPACE hash

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
/* #undef VERSION */

/* Stops putting the code inside the Google namespace */
#define _END_GOOGLE_NAMESPACE_ }

/* Puts following code inside the Google namespace */
#define _START_GOOGLE_NAMESPACE_ namespace nopacoHASH {
