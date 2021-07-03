cmake_minimum_required(VERSION 3.19)
project(2imv15_project_2)

set(CMAKE_CXX_STANDARD 17)

include_directories(include/common)
include_directories(include/linux)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(PNG REQUIRED)
include_directories(library ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

set(HEADER_FILES
        include/linux/gfx/config.h
        include/linux/gfx/geom3d.h
        include/linux/gfx/geom4d.h
        include/linux/gfx/gfx.h
        include/linux/gfx/mat2.h
        include/linux/gfx/mat3.h
        include/linux/gfx/mat4.h
        include/linux/gfx/quat.h
        include/linux/gfx/symmat2.h
        include/linux/gfx/symmat3.h
        include/linux/gfx/symmat4.h
        include/linux/gfx/vec2.h
        include/linux/gfx/vec3.h
        include/linux/gfx/vec4.h)

set(SOURCE_FILES
        src/demo.cpp
        src/solver.cpp
        src/Square.cpp
        src/SquareRigid.cpp
        src/Boundary.cpp
        src/BaseObject.cpp)

#set(LIBRARY_FILES
#        lib/x64/libfreeglut.a
#        lib/x64/libfreeglut_static.a
#        lib/libfreeglut.a
#        lib/libfreeglut_static.a
#        lib/libpng.dll.a
#        lib/libpng.la
#        lib/libpng12.a
#        lib/libpng12.dll.a
#        lib/libpng12.la)

add_executable(2imv15_project_2 ${SOURCE_FILES} ${HEADER_FILES} ${LIBRARY_FILES})
target_link_libraries(2imv15_project_2 ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_link_libraries(2imv15_project_2 PNG::PNG)
