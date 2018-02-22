#!/bin/bash

cd /home/mlibydt3/galaxy/tools/

folder=$(echo $(readlink -f $1) | rev | cut -c 5- | rev)
folder_path=$folder$(echo _files)
cp $folder_path/* $2/