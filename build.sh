#!/bin/bash

flags="-Prelease"

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
  -p|--platform)
    flags="$flags -Pcloud=$2"
    shift
    ;;
  -f|--profiling)
    flags="$flags -Pprofiler"
    ;;
  *)
    # unknown option
    echo "Failed to recongize argument '$1'"
    exit 1
    ;;
  esac
  shift # past argument or value
done

rm -f build.log

# build falcon_genomics
cd ./falcon

echo "building falcon_genomics..."
./gradlew clean install $flags &>> ../build.log
version=$(cat ../build.log | grep "Version:" | awk '{print $2}')

echo "falcon genomics version: $version"
cd -

echo "building gatk..."
./gradlew bundle -Dfalcon.version="$version" &>> build.log

echo "installing gatk"
mkdir -p ./export/
cp build/libs/gatk.jar ./export/
