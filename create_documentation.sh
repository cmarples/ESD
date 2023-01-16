#!/bin/bash
# Generate LEOD documentation using Doxygen and open the index.html file.
# This script assumes that the current directory is the LEOD home directory.
cd docs
doxygen Doxyfile
cd html
explorer.exe index.html