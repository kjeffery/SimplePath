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
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
else()
    set(CMAKE_CXX_STANDARD 20)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma") # Only set this for axv2 or above
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
endif(MSVC)

set(SOURCES
        Cameras/Camera.h
        Image/Image.cpp
        Image/Image.h
        Integrators/Integrator.h
        Lights/Light.h
        base/Array2D.h
        base/Constants.h
        base/Endian.h
        base/FileParser.cpp
        base/FileParser.h
        base/Logger.cpp
        base/Logger.h
        base/PlyReader.cpp
        base/PlyReader.h
        base/ProgressBar.h
        base/RunningStats.h
        base/Scene.h
        base/SmartToRawPointerIterator.h
        base/Tile.h
        base/TileScheduler.h
        base/Util.h
        base/not_null.h
        main.cpp
        materials/Material.cpp
        materials/Material.h
        math/AffineSpace.h
        math/Angles.h
        math/BBox.h
        math/LinearSpace3x3.h
        math/Math.h
        math/Matrix4x4.h
        math/ONB.h
        math/Quaternion.h
        math/RGB.h
        math/Ray.h
        math/Sampler.h
        math/Sampling.cpp
        math/Sampling.h
        math/Vector3.h
        math/VectorType.h
        shapes/Aggregate.h
        shapes/Intersection.h
        shapes/ListAccelerator.h
        shapes/Plane.h
        shapes/Primitive.h
        shapes/Shape.h
        shapes/Sphere.h
        shapes/Triangle.cpp
        shapes/Triangle.h
)

add_executable(SimplePath ${SOURCES})

if (UNIX)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(SimplePath PRIVATE Threads::Threads atomic)
endif(UNIX)
