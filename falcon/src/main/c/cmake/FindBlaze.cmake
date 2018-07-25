ExternalProject_Add(blaze-download
    PREFIX "blaze"
    URL https://s3.amazonaws.com/fcs-build-public/blaze-v0.5.0.tar.gz
    URL_MD5 995366bbf937d66f8e8568e551c4aea8
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
