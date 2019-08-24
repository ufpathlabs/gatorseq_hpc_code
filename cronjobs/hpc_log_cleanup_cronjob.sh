#!/bin/bash

echo $(date)

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR
rm *.log.txt

echo "" 

