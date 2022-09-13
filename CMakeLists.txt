cmake_minimum_required(VERSION 3.7...3.23)

if (${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(SimplePath)

set(CMAKE_CXX_STANDARD 20)

add_executable(SimplePath main.cpp base/Util.h base/FileParser.cpp base/FileParser.h)

if (UNIX)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(SimplePath PRIVATE Threads::Threads atomic)
endif(UNIX)
