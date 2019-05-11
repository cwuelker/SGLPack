cmake_minimum_required( VERSION 2.6 )

set( CMAKE_BUILD_TYPE Release )
set( CMAKE_CXX_STANDARD 11 )

# Set your path preferences here
# ----------------------------------------------------------------------
set( PATH_TO_PROJECT /home/wuelker/www/NFSGLFT )
set( PATH_TO_INCLUDE /usr/local/include )
set( GSLPATH /home/wuelker/www/NFSGLFT/include/GSL/gsl-master )
# ----------------------------------------------------------------------

project( NFSGLFT )

add_executable( NFSGLFT example_NFSGLFT.cpp ${GSLPATH}/specfunc/legendre_poly.c ${GSLPATH}/specfunc/hyperg_1F1.c ${GSLPATH}/specfunc/bessel_J0.c ${GSLPATH}/specfunc/bessel_Jn.c ${GSLPATH}/specfunc/bessel_Jnu.c ${GSLPATH}/specfunc/bessel_olver.c ${GSLPATH}/specfunc/airy.c ${GSLPATH}/specfunc/airy_der.c ${GSLPATH}/specfunc/gamma.c ${GSLPATH}/err/error.c ${GSLPATH}/specfunc/hyperg.c ${GSLPATH}/specfunc/exp.c ${GSLPATH}/specfunc/laguerre.c ${GSLPATH}/specfunc/beta.c ${GSLPATH}/specfunc/psi.c ${GSLPATH}/complex/math.c ${GSLPATH}/specfunc/bessel_J1.c ${GSLPATH}/specfunc/bessel_I1.c ${GSLPATH}/specfunc/bessel_Inu.c ${GSLPATH}/specfunc/bessel_In.c ${GSLPATH}/specfunc/zeta.c ${GSLPATH}/specfunc/bessel_amp_phase.c ${GSLPATH}/specfunc/log.c ${GSLPATH}/specfunc/hyperg_U.c ${GSLPATH}/specfunc/poch.c ${GSLPATH}/specfunc/pow_int.c ${GSLPATH}/specfunc/bessel_temme.c ${GSLPATH}/specfunc/bessel_I0.c ${GSLPATH}/specfunc/bessel.c ${GSLPATH}/specfunc/result.c ${GSLPATH}/specfunc/trig.c ${GSLPATH}/err/message.c ${GSLPATH}/sys/infnan.c ${GSLPATH}/err/stream.c ${GSLPATH}/specfunc/elementary.c ${GSLPATH}/sys/fdiv.c ${GSLPATH}/specfunc/bessel_Knu.c ${GSLPATH}/specfunc/bessel_K0.c ${GSLPATH}/sys/coerce.c ./include/SpharmonicKit27/flt_classic.c ./include/SpharmonicKit27/cospmls.c ./include/SpharmonicKit27/primitive.c ./include/SpharmonicKit27/newFCT.c ./include/SpharmonicKit27/OURperms.c ./include/SpharmonicKit27/OURmods.c ./include/SpharmonicKit27/FST_semi_memo.c ./include/SpharmonicKit27/fft_grids.c ./include/SpharmonicKit27/naive_synthesis.c ./include/SpharmonicKit27/primitive_FST.c ./include/SpharmonicKit27/seminaive.c ./include/SpharmonicKit27/csecond.c ./include/SpharmonicKit27/oddweights.c ./include/SpharmonicKit27/weights.c ./include/SpharmonicKit27/FFTcode.c ./include/SpharmonicKit27/indextables.c ./include/SpharmonicKit27/permroots.c ./include/SpharmonicKit27/precomp_flt_classic.c ./include/clenshaw.cpp ./include/fftw.cpp ./include/nfft_use.cpp ./include/flt.cpp ./include/nfsglft.cpp )

link_directories( /usr/local/lib64 )

include_directories( ${PATH_TO_PROJECT}/include )
include_directories( ${PATH_TO_INCLUDE} )
include_directories( ${GSLPATH} )

target_link_libraries( NFSGLFT fftw3 nfft3 m stdc++ )
