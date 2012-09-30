#!/bin/sh

# This file downloads the chemcpp sources,

# 1) It checks wheter the directory "src/chemcpp" exists
# if NOT:
# 2) It downloads and extracts them into the chemcpp_1.0.2 folder
# 3) Patches them (a bug)
# 4) It moves the sources to the "src" folder

#So: if you want to try out a new/different version of chemcpp, just put it unter src/chemcpp
#If you are on windows, put the directory there yourself (and patch it)...

if [ ! -d ../src/chemcpp ]; then

rm -rf chemcpp_1.0.2
wget http://downloads.sourceforge.net/project/chemcpp/chemcpp/1.0.2/chemcpp_1.0.2.tar.gz
tar xfvz chemcpp_1.0.2.tar.gz
rm chemcpp_1.0.2.tar.gz
patch -p0 < chemcpp_1.0.2_uintbug.diff

cp ./Makefile ./chemcpp_1.0.2/src/
cp ./constant.cpp ./chemcpp_1.0.2/src/

grep -v "#define CHEMCPPPATH" ./chemcpp_1.0.2/src/constant.h > ./chemcpp_1.0.2/src/constant.h.tmp
mv ./chemcpp_1.0.2/src/constant.h.tmp ./chemcpp_1.0.2/src/constant.h
cat ./constant.h.append >> ./chemcpp_1.0.2/src/constant.h

mv ./chemcpp_1.0.2 ../src/chemcpp


#data-directory
#(no symbolic links on windows...)
rm -f ../inst/chemcpp/data/*
cp ../src/chemcpp/data/*.csv ../inst/chemcpp/data/


fi

#you might want to strip chemcpp afterwards...




