#!/bin/bash

if [ $# -ne 1 ]; then echo "usage: $0 [logFile]"; exit; fi
f=$1

function grab { 
  echo $1 $(
    sed -n "/$1 :/,/}/p" $f |\
    grep -wE $2 |\
    sed "s/^.*$1/$1/g" |\
    sed 's/;$//g' 
  )
}

grab run id
grab laser w
grab laser x
grab laser y
grab mapmt hv
grab mapmt id
grab maroc id
grab log events
