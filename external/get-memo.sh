#!/bin/bash

# memo

rm -rf memo

git clone https://github.com/jimporter/memo.git

pushd memo > /dev/null

  if [ $? != "0" ]
  then
    exit 1
  fi

  mkdir -p ../include/memo
  cp -rp include/memo.hpp  ../include/memo

popd > /dev/null

