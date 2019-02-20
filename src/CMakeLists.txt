file(GLOB SOURCES "*.cpp")
add_library(nfm SHARED ${SOURCES})
add_library(nfm_static STATIC ${SOURCES})
