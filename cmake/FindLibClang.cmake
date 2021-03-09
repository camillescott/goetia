#.rst:
# FindLibClang
# ------------
#
# Find LibClang
#
# Find LibClang headers and library
#
# ::
#
#   LibClang_FOUND             - True if libclang is found.
#   LibClang_LIBRARY           - Clang library to link against.
#   LibClang_VERSION           - Version number as a string (e.g. "3.9").

#
# Python support for clang might not be available for Python3. We need to
# find what we have.
#
find_library(LibClang_LIBRARY libclang.so 
             PATHS $ENV{CONDA_PREFIX}/lib/ /usr/lib/llvm-6.0/lib/ /usr/lib/llvm-7/lib/ x86_64-linux-gnu )
message(STATUS ${LibClang_LIBRARY})

if(LibClang_LIBRARY)
    set(LibClang_LIBRARY ${LibClang_LIBRARY})
    string(REGEX REPLACE ".*clang-\([0-9]+.[0-9]+\).*" "\\1" LibClang_VERSION_TMP "${LibClang_LIBRARY}")
    set(LibClang_VERSION ${LibClang_VERSION_TMP} CACHE STRING "LibClang version" FORCE)
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibClang REQUIRED_VARS  LibClang_LIBRARY 
                                           VERSION_VAR    LibClang_VERSION)

find_program(CLANG_EXE clang++)
EXECUTE_PROCESS( COMMAND ${CLANG_EXE} --version OUTPUT_VARIABLE clang_full_version_string )
string (REGEX REPLACE ".*clang version ([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" CLANG_VERSION_STRING ${clang_full_version_string})

set(CLANG_VERSION_STRING ${CLANG_VERSION_STRING} PARENT_SCOPE)

mark_as_advanced(LibClang_VERSION)
unset(_filename)
unset(_find_libclang_filename)
