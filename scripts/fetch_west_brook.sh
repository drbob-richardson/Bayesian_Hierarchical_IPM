#!/usr/bin/env bash
# Fetch the West Brook brook trout PIT-tag dataset from USGS ScienceBase.
# Public domain (CC0). DOI: 10.5066/P14PDHXM. ~30 MB; not stored in git.
set -euo pipefail
dest="data/west_brook"
item="66cc8fe5d34e98e8a9243552"
base="https://www.sciencebase.gov/catalog/file/get/${item}"
mkdir -p "$dest"
curl -fSL -o "$dest/cdWB_electro_DR.csv"            "${base}?name=cdWB_electro_DR.csv"
curl -fSL -o "$dest/WestBrookTroutData_bhl_Final.xml" "${base}?name=WestBrookTroutData_bhl_Final.xml"
echo "Downloaded West Brook data to $dest/"
