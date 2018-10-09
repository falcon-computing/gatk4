ExternalProject_Add(pairhmm-client-download
    PREFIX "pairhmm-client"
    URL https://s3.amazonaws.com/fcs-build-public/pairhmm-client.tgz
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${CMAKE_BINARY_DIR}/pairhmm-client/install"
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(PMMClient)
add_dependencies(PMMClient pairhmm-client-download)

set(PMMClient_DIR "${CMAKE_BINARY_DIR}/pairhmm-client/install")
set(PMMClient_INCLUDE_DIRS "${PMMClient_DIR}/include")
set(PMMClient_LIBRARY_DIRS "${PMMClient_DIR}/lib")
set(PMMClient_LIBRARIES pmm_client pmm_interface)
