cmake_minimum_required(VERSION 3.1)
project(relax-field)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
message("build with libigl")
option(LIBIGL_GLFW                "Build target igl::glfw"                ON)
option(LIBIGL_PNG                 "Build target igl::png"                 OFF)
option(LIBIGL_IMGUI               "Build target igl::imgui"               ON)
option(LIBIGL_OPENGL              "Build target igl::opengl"              OFF)
option(LIBIGL_PREDICATES          "Build target igl::predicates"          OFF)
option(LIBIGL_COPYLEFT_COMISO     "Build target igl_copyleft::comiso"     ON)
option(LIBIGL_COPYLEFT_CORE       "Build target igl_copyleft::core"       ON)
option(LIBIGL_COPYLEFT_TETGEN     "Build target igl_copyleft::tetgen"     OFF)
option(LIBIGL_COPYLEFT_CGAL       "Build target igl_copyleft::cgal"       OFF)
option(LIBIGL_RESTRICTED_TRIANGLE "Build target igl_restricted::triangle" OFF)
include(libigl)


# find_package (SuiteSparse 5.12 NO_MODULE)
message("build with suite-sparse")
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::glfw igl::imgui igl_copyleft::comiso ${SUITESPARSE_LIBRARIES})
