#!/bin/bash
#TZ=":CET"
DATE=`date +%Y.%-m.%-d`
DATETIME=`date +%Y.%m.%d\ at\ %H:%M`
sed -i "s/^Version: .*$/Version: $DATE/" DESCRIPTION
head DESCRIPTION

sed -i "s/packageStartupMessage(\"Version.*$/packageStartupMessage(\"Version: $DATETIME\")/" R/3_onAttach.R
