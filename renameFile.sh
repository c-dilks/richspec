#!/bin/bash
# append datadir name to a pdf file

if [ $# -ne 1 ]; then exit; fi

dataname=$(readlink -f datadir | sed 's/\/$//' | sed 's/^.*\///g')
fname=$(echo $1 | sed 's/\.pdf$//g')
mv -v ${fname}{,.${dataname}}.pdf
