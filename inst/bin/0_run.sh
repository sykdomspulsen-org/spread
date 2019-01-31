#!/bin/bash

COMPUTER=$(cat /tmp/computer)

(
  flock -n 200 || exit 1

  source /etc/environment

  echo
  echo
  echo
  echo
  echo "****START****noispiah****"

  /usr/local/bin/Rscript /r/noispiah/src/RunProcess.R

  echo "****END****noispiah****"
  
) 200>/var/lock/.noispiah.exclusivelock
