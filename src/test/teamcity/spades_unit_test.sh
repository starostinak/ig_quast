#!/bin/bash

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


echo "### PREPROCESSING ###"

./prepare_cfg
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "prepare_cfg finished with exit code $errlvl"
    exit $errlvl
fi

make rdt
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "make rdt finished with exit code $errlvl"
    exit $errlvl
fi

echo "### RUNNING ###"

set -e
./run rdt

errlvl=$?

echo "### TEAMCITY INVOKATION COMPLETE, EXIT CODE = $errlvl ###"

exit $errlvl

