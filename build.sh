#!/bin/bash

rm -f build.log

# build falcon_genomics
cd ./falcon

echo "building falcon_genomics..."
if [ $# -gt 0 ]; then
	  cloud=$1
	    ./gradlew clean install -Prelease -Pcloud=$cloud &>> ../build.log
	else
		  ./gradlew clean install -Prelease &>> ../build.log
	  fi
	  version=$(cat ../build.log | grep "Version:" | awk '{print $2}')

	  echo "falcon genomics version: $version"
	  cd -

	  echo "building gatk..."
	  #mvn clean package -Ddisable.queue -Dfalcon-genomics.version="$version" &>> build.log
      ./gradlew bundle -Dfalcon.version="$version" &>> build.log
	  echo "installing gatk to ~/.falcon-genome/gatk/$version"
	  mkdir -p ~/.falcon-genome/gatk/$version
	  cp build/libs/gatk.jar  ~/.falcon-genome/gatk/$version/
