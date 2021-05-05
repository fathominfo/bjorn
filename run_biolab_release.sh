#!/bin/bash
outdir=$1
metafn=$2
date=$3
# TEST RUN
# python3 src/biolabs_release.py
# DEV RUN
python3 src/biolabs_release.py --out-dir $outdir --date $date --metadata $metafn
# REAL RUN /Users/al/Documents/scripps/analysis/jordan/biolabs_metadata_template.xls
# python3 src/biolabs_release.py
