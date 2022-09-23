cmake_minimum_required(VERSION 3.7...3.23)

if (${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(SimplePath)

set(CMAKE_CXX_STANDARD 23)

if (MSVC)
    # We need the latest to get things like ranges::views
    set(CMAKE_CXX_STANDARD 23)

    # Windows doesn't define these: it only sets the expected definitions for __AVX__ and up.
    # https://docs.microsoft.com/en-us/cpp/preprocessor/predefined-macros
    add_definitions(/D__SSE4_1__)
    add_definitions(/D__SSE4_2__)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
else()
    set(CMAKE_CXX_STANDARD 20)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma") # Only set this for axv2 or above
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
endif(MSVC)

add_executable(SimplePath main.cpp base/Util.h base/FileParser.cpp base/FileParser.h base/Scene.h math/Vector3.h base/Constants.h math/Math.h math/Matrix4x4.h math/BBox.h math/ONB.h math/LinearSpace3x3.h math/Quaternion.h math/AffineSpace.h shapes/Shape.h shapes/Intersection.h math/Ray.h shapes/Sphere.h shapes/ListAccelerator.h base/not_null.h shapes/Aggregate.h shapes/Primitive.h materials/Material.h base/SmartToRawPointerIterator.h base/Array2D.h math/RGB.h base/Endian.h Image/Image.h math/Sampler.h math/VectorType.h base/Tile.h base/TileScheduler.h)

if (UNIX)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(SimplePath PRIVATE Threads::Threads atomic)
endif(UNIX)
