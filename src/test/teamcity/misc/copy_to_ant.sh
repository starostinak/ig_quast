#!/bin/bash

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

cd /tmp/
mkdir -p data/output
cp -rv /Johnny/data/ANT_BACKUP/ant$1/input ./data/
chmod -R 777 data/

