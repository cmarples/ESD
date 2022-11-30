# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 11:10:40 2022
"""

import os.path as path
import matplotlib as mpl
print("Your style sheets are located at: {}".format(path.join(mpl.__path__[0], 'mpl-data', 'stylelib')))