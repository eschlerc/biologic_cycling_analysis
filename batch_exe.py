# Author: Chris Eschler
# This script batch executes the navani_analysis scripts to quickly analyze all data in a folder.
# Paste the navani_analysis directory into the top level of the data directory.

import os
import shutil

for rootpath, filepaths, filenames in os.walk('/home/eschlerc/Dropbox (MIT)/MIT/_Grad/Thesis/Data/GCPL/PROPEL-1K'):
    print(filepaths)