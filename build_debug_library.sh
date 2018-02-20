#!/bin/bash

source config.sh

\rm -f *.so
cd src/
\rm -f *.o *.so
echo "$CC $FLAGS $DEBUGFLAGS -fpic -c *.cpp"
$CC $FLAGS $DEBUGFLAGS -fpic -c *.cpp

case ${OS_NAME} in
    "Darwin")
        ROOT_FOLDER=$(dirname $(pwd))
        echo "$CC $FLAGS $DEBUGFLAGS -shared -install_name ${ROOT_FOLDER}/lib${LIBNAME}.so  -o lib${LIBNAME}.so *.o"
        $CC $FLAGS $DEBUGFLAGS -shared -install_name ${ROOT_FOLDER}/lib${LIBNAME}.so -o lib${LIBNAME}.so *.o
        ;;
    "Linux")
        echo "$CC $FLAGS $DEBUGFLAGS -shared -o lib${LIBNAME}.so *.o"
        $CC $FLAGS $DEBUGFLAGS -shared -o lib${LIBNAME}.so *.o
        ;;
esac

mv lib${LIBNAME}.so ../
cd ..

echo
echo "Library ready!"
echo
echo "Help, how can I use it?"
case ${OS_NAME} in
    "Linux")
        echo "1)   $CC $FLAGS -I$(pwd)/src/ -c example.cpp"
        echo "     $CC $FLAGS -L$(pwd) example.o -l${LIBNAME}"
        echo "2)   $CC $FLAGS -I$(pwd)/src/ -L$(pwd) example.cpp -l${LIBNAME}"
        ;;
    "Darwin")
        echo "1)   $CC $FLAGS -I$(pwd)/src/ -c example.cpp"
        echo "     $CC $FLAGS -I$(pwd)/src/ -L$(pwd) example.o -l${LIBNAME}"
        echo "2)   $CC $FLAGS -I$(pwd)/src/ -L$(pwd) example.cpp -l${LIBNAME}"
        ;;
esac
