#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
s3_dir=s3://fcs-data/gatk3-build/testdata
aws s3 sync $s3_dir/ $DIR/
