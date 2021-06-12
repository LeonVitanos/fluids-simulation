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
        include/common/AngularSpringForce.h
        include/common/CircularWireConstraint.h
        include/common/imageio.h
        include/common/linearSolver.h
        include/common/Particle.h
        include/common/RodConstraint.h
        include/common/SlidingConstraint.h
        include/common/SpringForce.h
        include/common/Wall.h
        include/common/ConstraintSolver.h
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
        src/AngularSpringForce.cpp
        src/CircularWireConstraint.cpp
        src/imageio.cpp
        src/linearSolver.cpp
        src/Particle.cpp
        src/RodConstraint.cpp
        src/SpringForce.cpp
        src/Gravity.cpp
        src/Force.cpp
        src/Scene.cpp
        src/SlidingConstraint.cpp
        src/FixedConstraint.cpp
        src/ConstraintSolver.cpp
        src/Wall.cpp
        src/demo.cpp
        src/solver.cpp)

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
