#!/bin/sh
if [ ! ${0%%/*} ]; then
    my_dir=$(dirname "$0")
else
    my_dir=$(dirname "$PWD/$0")
fi
install_dir=$(dirname "$my_dir")
echo $install_dir
OS=$(uname)

#setup env
if [ $OS = "Linux" ]; then
    export LD_LIBRARY_PATH=$install_dir/lib
fi
export FONTCONFIG_FILE=$install_dir/etc/fonts/fonts.conf
export XDG_DATA_HOME=$install_dir/share

# run gnuplot
if [ "$#" = "0" ]; then
    $install_dir/libexec/gnuplot
else
    $install_dir/libexec/gnuplot $1 > /dev/null 2>&1
fi
