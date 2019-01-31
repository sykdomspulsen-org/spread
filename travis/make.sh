#!/bin/bash
DATE=`date +%Y.%m.%d`
sed -i "s/^Version: .*$/Version: $DATE/" DESCRIPTION
head DESCRIPTION

sed -i "s/packageStartupMessage(\"Version.*$/packageStartupMessage(\"Version $DATE\")/" R/onAttach.R
## Other options:
## Only add if the commit is tagged: so something like:
#if [ $TRAVIS_TAG ] ; then
#   addToDrat
#fi
##but will need to edit .travis.yml since $TRAVIS_BRANCH will now equal $TRAVIS_TAG
