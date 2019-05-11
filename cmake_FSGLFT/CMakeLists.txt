cmake_minimum_required( VERSION 2.6 )

set( CMAKE_BUILD_TYPE Release )
set( CMAKE_CXX_STANDARD 11 )

# Set your path preferences here
# ----------------------------------------------------------------------
set( PATH_TO_PROJECT /home/wuelker/www/SGLPack )
set( PATH_TO_SPHARMONIC_KIT ${PATH_TO_PROJECT}/include/SpharmonicKit27 )
# ----------------------------------------------------------------------

project( SGLPack )

add_executable( SGLPack example_FSGLFT.cpp ${PATH_TO_PROJECT}/include/fsft.cpp ${PATH_TO_PROJECT}/include/fsglft.cpp ${PATH_TO_SPHARMONIC_KIT}/cospmls.c ${PATH_TO_SPHARMONIC_KIT}/primitive.c ${PATH_TO_SPHARMONIC_KIT}/newFCT.c ${PATH_TO_SPHARMONIC_KIT}/OURperms.c ${PATH_TO_SPHARMONIC_KIT}/OURmods.c ${PATH_TO_SPHARMONIC_KIT}/FST_semi_memo.c ${PATH_TO_SPHARMONIC_KIT}/fft_grids.c ${PATH_TO_SPHARMONIC_KIT}/naive_synthesis.c ${PATH_TO_SPHARMONIC_KIT}/primitive_FST.c ${PATH_TO_SPHARMONIC_KIT}/seminaive.c ${PATH_TO_SPHARMONIC_KIT}/csecond.c ${PATH_TO_SPHARMONIC_KIT}/oddweights.c ${PATH_TO_SPHARMONIC_KIT}/weights.c ${PATH_TO_SPHARMONIC_KIT}/FFTcode.c ${PATH_TO_SPHARMONIC_KIT}/indextables.c ${PATH_TO_SPHARMONIC_KIT}/permroots.c )

include_directories( ${PATH_TO_PROJECT}/include )
include_directories( ${PATH_TO_SPHARMONIC_KIT} )
