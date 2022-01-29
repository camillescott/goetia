#.rst:
# FindCppyy
# -------
#
# Find Cppyy
#
# This module finds an installed Cppyy.  It sets the following variables:
#
# ::
#
#   Cppyy_FOUND - set to true if Cppyy is found
#   Cppyy_DIR - the directory where Cppyy is installed
#   Cppyy_EXECUTABLE - the path to the cppyy-generator executable
#   Cppyy_INCLUDE_DIRS - Where to find the ROOT header files.
#   Cppyy_VERSION - the version number of the Cppyy backend.
#
#
# The module also defines the following functions:
#
#   cppyy_add_bindings - Generate a set of bindings from a set of header files.
#
# The minimum required version of Cppyy can be specified using the
# standard syntax, e.g.  find_package(Cppyy 4.19)
#
#

execute_process(COMMAND cling-config --cmake OUTPUT_VARIABLE CPPYY_MODULE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)

if(CPPYY_MODULE_PATH)
    #
    # Cppyy_DIR: one level above the installed cppyy cmake module path
    #
    set(Cppyy_DIR ${CPPYY_MODULE_PATH}/../)
    #
    # Cppyy_INCLUDE_DIRS: Directory with cppyy headers
    #
    set(Cppyy_INCLUDE_DIRS ${Cppyy_DIR}include)
    #
    # Cppyy_VERSION.
    #
    find_package(ROOT QUIET REQUIRED PATHS ${CPPYY_MODULE_PATH})
    if(ROOT_FOUND)
        set(Cppyy_VERSION ${ROOT_VERSION})
    endif()
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cppyy
                                  REQUIRED_VARS ROOT_genreflex_CMD Cppyy_DIR Cppyy_INCLUDE_DIRS CPPYY_MODULE_PATH
                                  VERSION_VAR Cppyy_VERSION
)
mark_as_advanced(Cppyy_VERSION)

# Get the cppyy libCling library. Not sure if necessary?
find_library(LibCling_LIBRARY libCling.so PATHS ${Cppyy_DIR}/lib)


#
# Generate a set of bindings from a set of header files. Somewhat like CMake's
# add_library(), the output is a compiler target. In addition ancilliary files
# are also generated to allow a complete set of bindings to be compiled,
# packaged and installed.
#
#   cppyy_add_bindings(
#       pkg
#       [URL url]
#       [LICENSE license]
#       [LANGUAGE_STANDARD std]
#       [GENERATE_OPTIONS option...]
#       [COMPILE_OPTIONS option...]
#       [INCLUDE_DIRS dir...]
#       [LINK_LIBRARIES library...]
#
# The bindings are based on https://cppyy.readthedocs.io/en/latest/, and can be
# used as per the documentation provided via the cppyy.gbl namespace. First add
# the directory of the <pkg>.rootmap file to the LD_LIBRARY_PATH environment
# variable, then "import cppyy; from cppyy.gbl import <some-C++-entity>".
#
# Alternatively, use "import <pkg>". This convenience wrapper supports
# "discovery" of the available C++ entities using, for example Python 3's command
# line completion support.
#
#
# Arguments and options:
#
#   pkg                 The name of the package to generate. This can be either
#                       of the form "simplename" (e.g. "Akonadi"), or of the
#                       form "namespace.simplename" (e.g. "KF5.Akonadi").
#
#
#
#   LANGUAGE_STANDARD std
#                       The version of C++ in use, "14" by default.
#
#   INTERFACE_FILE      Header to be passed to genreflex. Should contain template
#                       specialization declarations if required.
#
#   INTERFACE_HEADERS   List of headers to run cppyy-generator on.
#
#   HEADERS             Library headers from which to generate the map. Should match up with
#                       interface file includes.
#
#   SELECTION_XML       selection XML file passed to genreflex.
#
#
#   
#
#   GENERATE_OPTIONS option
#                       Options which will be passed to the rootcling invocation
#                       in the cppyy-generate utility. cppyy-generate is used to
#                       create the bindings map.
# 
#  
#
#
#   COMPILE_OPTIONS option
#                       Options which are to be passed into the compile/link
#                       command.
#
#   INCLUDE_DIRS dir    Include directories.
#
#   LINK_LIBRARIES library
#                       Libraries to link against.
#

#
#
# Returns via PARENT_SCOPE variables:
#
#   CPPYY_LIB_TARGET    The target cppyy bindings shared library.
#
#
function(cppyy_add_bindings pkg)
    set(simple_args LANGUAGE_STANDARD INTERFACE_FILE SELECTION_XML )
    set(list_args HEADERS INTERFACE_HEADERS COMPILE_OPTIONS INCLUDE_DIRS LINK_LIBRARIES 
        GENERATE_OPTIONS NAMESPACES)

    cmake_parse_arguments(
        ARG
        ""
        "${simple_args}"
        "${list_args}"
        ${ARGN}
    )
    if(NOT "${ARG_UNPARSED_ARGUMENTS}" STREQUAL "")
        message(SEND_ERROR "Unexpected arguments specified '${ARG_UNPARSED_ARGUMENTS}'")
    endif()
    string(REGEX MATCH "[^\.]+$" pkg_simplename ${pkg})
    string(REGEX REPLACE "\.?${pkg_simplename}" "" pkg_namespace ${pkg})
    set(lib_name "${pkg_namespace}${pkg_simplename}Cppyy")
    set(lib_file ${CMAKE_SHARED_LIBRARY_PREFIX}${lib_name}${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(cpp_file ${CMAKE_CURRENT_BINARY_DIR}/${pkg_simplename}.cpp)
    set(pcm_file ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${lib_name}_rdict.pcm)
    set(rootmap_file ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${lib_name}.rootmap)
    set(extra_map_file ${CMAKE_CURRENT_BINARY_DIR}/${pkg_simplename}.map)

    #
    # Language standard.
    #
    if("${ARG_LANGUAGE_STANDARD}" STREQUAL "")
        set(ARG_LANGUAGE_STANDARD "14")
    endif()

    #
    # Includes
    #
    foreach(dir ${ARG_INCLUDE_DIRS})
        list(APPEND includes "-I${dir}")
    endforeach()
    
    find_package(LibClang REQUIRED)

    #
    # Set up genreflex args.
    #
    set(genreflex_cxxflags "--extra-cxxflags")
    list(APPEND genreflex_cxxflags "-std=c++${ARG_LANGUAGE_STANDARD} -march=native")

    set(genreflex_args)
    if("${ARG_INTERFACE_FILE}" STREQUAL "")
        message(SEND_ERROR "No Interface specified")
    endif()
    list(APPEND genreflex_args "${ARG_INTERFACE_FILE}")
    if(NOT "${ARG_SELECTION_XML}" STREQUAL "")
        list(APPEND genreflex_args "--selection=${ARG_SELECTION_XML}")
    endif()

    list(APPEND genreflex_args "-o" "${cpp_file}")
    list(APPEND genreflex_args "--rootmap=${rootmap_file}")
    list(APPEND genreflex_args "--rootmap-lib=${lib_file}")
    list(APPEND genreflex_args "-l" "${lib_file}")
    foreach(dir ${includes})
        list(APPEND genreflex_args "${dir}")
    endforeach(dir)

    #
    # run genreflex
    #
    add_custom_command(
        OUTPUT ${cpp_file} ${rootmap_file} ${pcm_file}
        COMMAND ${CMAKE_SOURCE_DIR}/build-utils/genreflex-wrapper.py ${genreflex_cxxflags} ${genreflex_args}
        DEPENDS ${ARG_INTERFACE_FILE} ${ARG_HEADERS} ${ARG_SELECTION_XML}
    )

    #
    # Set up cppyy-generator args.
    #
    list(APPEND ARG_GENERATE_OPTIONS "-std=c++${ARG_LANGUAGE_STANDARD}")
    list(APPEND ARG_GENERATE_OPTIONS "-march=native")
    if(${CONDA_ACTIVE})
        # get the libcxx headers installed by conda
        execute_process(COMMAND bash -c "cat `find $CONDA_PREFIX -iname 'libcxx-*.json' -type f` | python -c \"import sys, json; print(json.load(sys.stdin)['extracted_package_dir'], end='')\""
                        OUTPUT_VARIABLE conda_libcxx_header_base)
        message(STATUS "adding conda libcxx includes to cppyy-generator options (${conda_libcxx_header_base})")
        list(APPEND ARG_GENERATE_OPTIONS "-I${conda_libcxx_header_base}/include/c++/v1")

        # now the headers from libclang...
        set(CLANGDEV_INCLUDE $ENV{CONDA_PREFIX}/lib/clang/${CLANG_VERSION_STRING}/include)
        message(STATUS "adding conda clangdev includes to cppyy-generator options (${CLANGDEV_INCLUDE})")
        list(APPEND ARG_GENERATE_OPTIONS "-I${CLANGDEV_INCLUDE}")

    endif()
    foreach(dir ${includes})
        list(APPEND ARG_GENERATE_OPTIONS "${dir}")
    endforeach(dir)
    #
    # Run cppyy-generator. First check dependencies. TODO: temporary hack: rather
    # than an external dependency, enable libclang in the local build.
    #
    get_filename_component(Cppyygen_EXECUTABLE ${ROOT_genreflex_CMD} DIRECTORY)
    set(Cppyygen_EXECUTABLE ${Cppyygen_EXECUTABLE}/cppyy-generator)

    set(generator_args)
    foreach(arg IN LISTS ARG_GENERATE_OPTIONS)
        string(REGEX REPLACE "^-" "\\\\-" arg ${arg})
        list(APPEND generator_args ${arg})
    endforeach()

    foreach(interface_header ${ARG_INTERFACE_HEADERS})
        string(REPLACE "${GOETIA_INCLUDE_ROOT}" "" header_basename "${interface_header}")
        set(header_map "${CMAKE_CURRENT_BINARY_DIR}/maps/${header_basename}.map")
        list(APPEND header_maps ${header_map})
        add_custom_command(OUTPUT ${header_map}
                           COMMAND python ${Cppyygen_EXECUTABLE} 
                                   --flags "\"${generator_args}\""
                                   ${header_map} ${interface_header}
                           DEPENDS ${interface_header} 
        )
    endforeach()
    add_custom_command(OUTPUT ${extra_map_file}
                       COMMAND python ${CMAKE_SOURCE_DIR}/build-utils/merge-cppyy-maps.py ${header_maps} -o ${extra_map_file}
                       DEPENDS ${header_maps}
    )
    #
    # Compile/link.
    #
    add_library(${lib_name} SHARED ${cpp_file} ${pcm_file} ${rootmap_file} ${extra_map_file})
    set_target_properties(${lib_name} PROPERTIES LINKER_LANGUAGE CXX)
    set_property(TARGET ${lib_name} PROPERTY VERSION ${version})
    set_property(TARGET ${lib_name} PROPERTY CXX_STANDARD ${ARG_LANGUAGE_STANDARD})
    set_property(TARGET ${lib_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${pkg_dir})
    set_property(TARGET ${lib_name} PROPERTY LINK_WHAT_YOU_USE TRUE)
    set_property(TARGET ${lib_name} PROPERTY VISIBILITY_INLINES_HIDDEN 0)
    # set_property(TARGET ${lib_name} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    target_include_directories(${lib_name} PRIVATE ${ARG_INCLUDE_DIRS} ${Cppyy_INCLUDE_DIRS})
    if (NOT "${ARG_COMPILE_OPTIONS}" STREQUAL "")
        target_compile_options(${lib_name} PRIVATE ${ARG_COMPILE_OPTIONS})
    endif()
    target_link_libraries(${lib_name} PUBLIC ${LibCling_LIBRARY} 
                          -Wl,--whole-archive ${ARG_LINK_LIBRARIES} -Wl,--no-whole-archive)

    set(LIBCLING         ${LibCling_LIBRARY} PARENT_SCOPE)
    set(CPPYY_LIB_TARGET ${lib_name} PARENT_SCOPE)
    set(CPPYY_ROOTMAP    ${rootmap_file} PARENT_SCOPE)
    set(CPPYY_PCM        ${pcm_file} PARENT_SCOPE)
    set(CPPYY_EXTRA_MAP  ${extra_map_file} PARENT_SCOPE)

endfunction(cppyy_add_bindings)



