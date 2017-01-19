## Process this file with autoconf to produce configure.
## In general, the safest way to proceed is to run ./autogen.sh

# make sure we're interpreted by some minimal autoconf
AC_PREREQ(2.57)

AC_INIT(nopaco, 0.99.5, R.Kuiper)

AC_CONFIG_HEADERS(src/densehash/sparseconfig.h)

AC_CONFIG_MACRO_DIRS([src/densehash/m4])

AX_CXX_COMPILE_STDCXX([11],,[mandatory])

# Check whether some low-level functions/files are available
AC_HEADER_STDC
AC_CHECK_FUNCS(memcpy memmove)
AC_CHECK_TYPES([uint16_t])     # defined in C99 systems
AC_CHECK_TYPES([u_int16_t])    # defined in BSD-derived systems, and gnu
AC_CHECK_TYPES([__uint16])     # defined in some windows systems (vc7)
AC_CHECK_TYPES([long long])    # probably defined everywhere, but...

# Figure out where hash_map lives and also hash_fun.h (or stl_hash_fun.h).
# This also tells us what namespace hash code lives in.
AC_CXX_STL_HASH
AC_CXX_STL_HASH_FUN

# Find out what namespace the user wants our classes to be defined in.
# TODO(csilvers): change this to default to sparsehash instead.
AC_DEFINE_GOOGLE_NAMESPACE(nopacoHASH)

# In unix-based systems, hash is always defined as hash<> (in namespace.
# HASH_NAMESPACE.)  So we can use a simple AC_DEFINE here.  On
# windows, and possibly on future unix STL implementations, this
# macro will evaluate to something different.)
AC_DEFINE(SPARSEHASH_HASH_NO_NAMESPACE, hash,
          [The system-provided hash function, in namespace HASH_NAMESPACE.])

# Do *not* define this in terms of SPARSEHASH_HASH_NO_NAMESPACE, because
# SPARSEHASH_HASH is exported to sparseconfig.h, but S_H_NO_NAMESPACE isn't.
AC_DEFINE(SPARSEHASH_HASH, HASH_NAMESPACE::hash,
          [The system-provided hash function including the namespace.])


# Write generated configuration file
AC_OUTPUT
