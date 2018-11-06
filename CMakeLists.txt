# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(noiseAna)
set(libname "noiseAna")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  ANITA_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)
#find_package(FFTW REQUIRED)


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR} )
include_directories($ENV{BOOST_INC})
include_directories($ENV{CANVAS_INC})
include_directories($ENV{CANVAS_ROOT_IO_INC})
include_directories($ENV{CETLIB_EXCEPT_INC})
include_directories($ENV{CETLIB_INC})
include_directories($ENV{FHICLCPP_INC})
include_directories($ENV{GALLERY_INC})
include_directories($ENV{LARCOREOBJ_INC})
include_directories($ENV{LARDATAOBJ_INC})
include_directories($ENV{NUSIMDATA_INC})
include_directories($ENV{ROOT_INC})

add_definitions(${ROOT_CXX_FLAGS})
link_directories(${UTIL_LIB_DIR})
link_directories($ENV{BOOST_LIB})
link_directories($ENV{CANVAS_LIB})
link_directories($ENV{CANVAS_ROOT_IO_LIB})
link_directories($ENV{CETLIB_EXCEPT_LIB})
link_directories($ENV{CETLIB_LIB})
link_directories($ENV{FHICLCPP_LIB})
link_directories($ENV{GALLERY_LIB})
link_directories($ENV{LARCOREOBJ_LIB})
link_directories($ENV{LARDATAOBJ_LIB})
link_directories($ENV{NUSIMDATA_LIB})

file(GLOB HEADER_FILES "AveragePowerSpectrum.h" )

file(GLOB_RECURSE SOURCE_FILES AveragePowerSpectrum.cxx)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(DICTNAME G__${libname})

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ROOT_LIBRARIES})


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()


message("DICTNAME is set to ${DICTNAME}" )
message("HEADER_FILES is set to ${HEADER_FILES}" )
message("SOURCE_FILES is set to ${SOURCE_FILES}" )

#---Create  a main program using the library
add_executable(makeSampleRayleigh makeSampleRayleigh.cxx)
target_link_libraries(makeSampleRayleigh RootFftwWrapper gallery ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} noiseAna)

add_executable(makeRayleighsRun makeRayleighsRun.cxx)
target_link_libraries(makeRayleighsRun RootFftwWrapper gallery ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} noiseAna)

add_executable(makeRayleighsAllChannels makeRayleighsAllChannels.cxx)
target_link_libraries(makeRayleighsAllChannels RootFftwWrapper gallery ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} noiseAna)

add_executable(makeRayleighsOneChannel makeRayleighsOneChannel.cxx)
target_link_libraries(makeRayleighsOneChannel RootFftwWrapper gallery ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} noiseAna)

