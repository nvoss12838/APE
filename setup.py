# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 14:59:18 2015

@author: radar
"""

import sys
from distutils.core import setup, Extension

compile_args = [-03]
setup(name='APE',
      version = '0.1.0',
      description = 'Automated PEST Enviornment',
      author = 'Nick Voss',
      author_email = 'nvoss@mail.usf.edu',
      licsence = 'MIT',
      py_modules = ['pcf'],
      zip_safe = False,
      )