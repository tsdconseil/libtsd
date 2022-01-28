# D'après https://cmake.org/cmake/help/latest/guide/tutorial/index.html
# et http://sirien.metz.supelec.fr/depot/SIR/TutorielCMake/index.html


# Faire ça sous Windows ?
# cmake -G "MinGW Makefiles" ..
# mingw32-make.exe

cmake_minimum_required(VERSION 3.8)


set(CMAKE_BUILD_TYPE Release)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(tsd VERSION 1.0)

option(BUILD_TESTS "Build auto-tests" OFF)


##############################
# Dépendances obligatoires
##############################
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (fmt REQUIRED)

##############################
# Dépendances optionnelles (..._FOUND)
##############################
find_package(Freetype)
find_package(PkgConfig)

if(PkgConfig_FOUND)
pkg_check_modules(GTKMM gtkmm-3.0)
else()
set(GTKMM_found TRUE)
endif()

find_package(PNG)

message("Result of optional dependencies lookup:")
message(" - FreeType : " ${Freetype_FOUND})
message(" - GTKMM : " ${GTKMM_FOUND})
message(" - libpng : " ${PNG_FOUND})


if(Freetype_FOUND AND PNG_FOUND)

add_compile_definitions(LIBTSD_USE_FREETYPE=1)
add_compile_definitions(LIBTSD_USE_PNG=1)

if(GTKMM_FOUND)

message(" => Libtsd figures are fully enabled (including display with GTKMM).")

add_compile_definitions(LIBTSD_USE_GTKMM=1)
else()


message(" => Libtsd figures are partially enabled (only PNG saving, no real-time diplay).")

add_compile_definitions(LIBTSD_USE_GTKMM=0)
endif()

else()

message(" => Libtsd figures are disabled.")

add_compile_definitions(LIBTSD_USE_FREETYPE=0)
add_compile_definitions(LIBTSD_USE_PNG=0)
add_compile_definitions(LIBTSD_USE_GTKMM=0)
endif()


install(DIRECTORY include/ DESTINATION include/${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE liste_tests CONFIGURE_DEPENDS "tests/*.cc")

file(GLOB_RECURSE liste_sources CONFIGURE_DEPENDS "src/*.cc")

add_library(tsd ${liste_sources})
#add_library(tsd SHARED ${liste_sources})
#add_library(tsd-static STATIC ${liste_sources})

add_executable(tsd-tests ${liste_tests})

target_include_directories(tsd PUBLIC include PRIVATE ${FREETYPE_INCLUDE_DIRS})  
#target_include_directories(tsd-static PUBLIC include PRIVATE ${FREETYPE_INCLUDE_DIRS})  
target_include_directories(tsd-tests PUBLIC include PRIVATE ${FREETYPE_INCLUDE_DIRS})  

link_directories (${GTKMM_LIBRARY_DIRS})

include_directories(${GTKMM_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS})

target_link_libraries (tsd Eigen3::Eigen ${FREETYPE_LIBRARIES}  ${GTKMM_LIBRARIES} fmt ${PNG_LIBRARY})
#target_link_libraries (tsd-static Eigen3::Eigen ${FREETYPE_LIBRARIES}  ${GTKMM_LIBRARIES} fmt ${PNG_LIBRARY})
target_link_libraries (tsd-tests tsd Eigen3::Eigen ${FREETYPE_LIBRARIES}  ${GTKMM_LIBRARIES} fmt ${PNG_LIBRARY})

install (TARGETS tsd DESTINATION lib)

#file(
#  GLOB 
#  usage_examples
#  *.cc
#)

#foreach(f ${usage_examples})
#  get_filename_component(exampleName ${f} NAME_WE) 
#  add_executable (${exampleName} ${f}) 
#  target_link_libraries(${exampleName} toto)
#  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${exampleName}
#    DESTINATION bin
#    RENAME ${CMAKE_PROJECT_NAME}-${exampleName})
#endforeach(f)
                           



# add the executable
# add_executable(Tutorial tutorial.cxx)

