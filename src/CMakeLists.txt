file(GLOB header *.h)
file(GLOB source *.cc)

enable_language(Fortran)

set(MAXIMA_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/usr/share/maxima}")

add_custom_command(OUTPUT stvk_tet.f90
   COMMAND "maxima" ARGS --userdir=${MAXIMA_DIR} -b stvk_tet.mac
   DEPENDS stvk_tet.mac)

add_library(stvk SHARED ${header} ${source} stvk_tet.f90)
