file(GLOB header *.h)
file(GLOB source *.cc)

enable_language(Fortran)

set(MAXIMA_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/usr/share/maxima}")

add_custom_command(OUTPUT constitutive.f90
   COMMAND "maxima" ARGS --userdir=${MAXIMA_DIR} -b constitutive.mac
   DEPENDS constitutive.mac)

add_custom_command(OUTPUT AX_B.f90
   COMMAND "maxima" ARGS --userdir=${MAXIMA_DIR} -b AX_B.mac
   DEPENDS AX_B.mac)

# add_custom_command(OUTPUT reduced2rs.f90
#    COMMAND "maxima" ARGS --userdir=${MAXIMA_DIR} -b reduced2rs.mac
#    DEPENDS reduced2rs.mac)

add_library(elastic SHARED ${header} ${source} constitutive.f90 AX_B.f90 reduced2rs.f90)
target_link_libraries(elastic
${ARPACK_LIBRARIES}
)
# install the libraries
