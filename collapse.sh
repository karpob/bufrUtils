#!/bin/bash
# collapse input directory and move all bufr files to current directory.
find $1 -iname '*.bufr' -exec mv \{\} . \;
