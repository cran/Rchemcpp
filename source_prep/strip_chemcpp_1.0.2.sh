#!/bin/sh

# Removes everything from chemcpp which is not necessary for compiling the R-package

if [ -d ../src/chemcpp ]; then

rm ../src/chemcpp/* 2> /dev/null
rm -rf ../src/chemcpp/bin
rm -rf ../src/chemcpp/data
rm -rf ../src/chemcpp/debug
rm -rf ../src/chemcpp/doc
rm -rf ../src/chemcpp/examples
rm -rf ../src/chemcpp/html
rm -rf ../src/chemcpp/templates
rm -rf ../src/chemcpp/tools
rm -rf ../src/chemcpp/src/Makefile.am
rm -rf ../src/chemcpp/src/Makefile.in
rm -rf ../src/chemcpp/src/libchemcpp.so

echo "stripping done"

fi





