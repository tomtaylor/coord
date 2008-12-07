#! /bin/bash

PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/home/matt/bin
HOME=/home/matt
LANGUAGE=en_GB:en_US:en_GB:en
LANG=en_GB.UTF-8


OP="$1"
ARGS="$OP --inline --diagram --op $OP/doc"

if [ $1 == ""]; then
		echo 'No dir given.' > commit-docgen
else
		echo "Post commit: $1" > commit-docgen
		/usr/bin/rdoc $ARGS >> commit-docgen 2>&1
fi
