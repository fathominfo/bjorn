#!/bin/bash
metafn=$1
# TEST RUN
# python3 src/biolabs_release.py
# DEV RUN
python3 src/biolabs_release_v0.py --out-dir /Users/al/Documents/scripps/analysis/jordan/2021-02_run2_release3 --date 2021-02-17 --metadata $metafn
# REAL RUN /Users/al/Documents/scripps/analysis/jordan/biolabs_metadata_template.xls
# python3 src/biolabs_release.py
