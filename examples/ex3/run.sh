#!/bin/bash

source ../../config.sh

OS_NAME=$(uname)

\rm -f exe
\rm -f *.o

# project root directory
ROOT_FOLDER=$(dirname $(dirname $(pwd)))

#runtime dynamic library path
RPATH="${ROOT_FOLDER}:${MCI_FOLDER}:${NFM_FOLDER}:${FFNN_FOLDER}"

# Build the debugging main executable
echo "$CC $FLAGS $OPTFLAGS -pthread -Wall $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp"
$CC $FLAGS $OPTFLAGS -pthread -Wall $IMCI $INFM $IFFNN -I${ROOT_FOLDER}/src/ -I/usr/local/include -c *.cpp

# For Mac OS, the install name is wrong and must be corrected
case ${OS_NAME} in
   "Darwin")
      echo "$CC $FLAGS $OPTFLAGS -pthread -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBFFNN $LIBGSL"
      $CC $FLAGS $OPTFLAGS -pthread -L${ROOT_FOLDER} $LMCI $LNFM $LFFNN $LGSL -o exe *.o -l$LIBNAME $LIBMCI $LIBNFM $LIBFFNN $LIBGSL
      
      echo "install_name_tool -change lib${LIBNAME}.so ${ROOT_FOLDER}/lib${LIBNAME}.so exe"
      install_name_tool -change lib${LIBNAME}.so ${ROOT_FOLDER}/lib${LIBNAME}.so exe
      ;;
   "Linux")
      echo "$CC $FLAGS $OPTFLAGS -pthread $LMCI $LNFM $LFFNN -I${ROOT_FOLDER}/src -L$${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}" $LIBMCI $LIBNFM $LIBFFNN
      $CC $FLAGS $OPTFLAGS -pthread $LMCI $LNFM $LFFNN -I${ROOT_FOLDER}/src/ -L${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} $LIBMCI $LIBNFM $LIBFFNN
      ;;
esac

echo ""
# Run the debugging executable
echo "Ready to run!"
echo ""
echo "--------------------------------------------------------------------------"
echo ""
echo ""
echo ""
./exe
