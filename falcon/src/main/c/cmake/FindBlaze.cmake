ExternalProject_Add(blaze-download
    PREFIX "blaze"
    URL https://s3.amazonaws.com/fcs-build-public/blaze-v0.6.0.tgz
    URL_MD5 216d0a263c2a0c4788e813d6f4012b04
    SOURCE_DIR "${CMAKE_BINARY_DIR}/blaze/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(Blaze)
add_dependencies(Blaze blaze-download)

set(Blaze_DIR "${CMAKE_BINARY_DIR}/blaze/install")
set(Blaze_INCLUDE_DIRS "${Blaze_DIR}/include")
set(Blaze_LIBRARY_DIRS "${Blaze_DIR}/lib")
set(Blaze_LIBRARIES "blaze" "blaze_message")
