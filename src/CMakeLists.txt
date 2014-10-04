file(GLOB header *.h)
file(GLOB source *.cc)

enable_language(Fortran)

set(MAXIMA_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/usr/share/maxima}")

add_custom_command(OUTPUT constitutive.f90
   COMMAND "maxima" ARGS --userdir=${MAXIMA_DIR} -b constitutive.mac
   DEPENDS constitutive.mac)

add_library(elastic SHARED ${header} ${source} constitutive.f90)

# install the libraries
