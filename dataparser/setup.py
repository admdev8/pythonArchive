# ---------------------------------------------------------------------------------
# 	dataparser -> setup.py
#	Copyright (C) 2019 Michael Winters
#	github: https://github.com/mesoic
#	email:  mesoic@protonmail.com
# ---------------------------------------------------------------------------------
#	
#	Permission is hereby granted, free of charge, to any person obtaining a copy
#	of this software and associated documentation files (the "Software"), to deal
#	in the Software without restriction, including without limitation the rights
#	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#	copies of the Software, and to permit persons to whom the Software is
#	furnished to do so, subject to the following conditions:
#	
#	The above copyright notice and this permission notice shall be included in all
#	copies or substantial portions of the Software.
#	
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#	SOFTWARE.
#

#!/usr/bin/env python 
# -*- coding: utf-8 -*-
import sys

# dataparser contains a series of unmaintained utilities for working with 
# data from electrical measurements. The scripts are archived as example 
# material in the event that such utilites are required in the future. 
try:
	from setuptools import setup

except ImportError:
	print('Please install or upgrade setuptools or pip to continue')
	sys.exit(1)

setup(name='minispice',
		description='Data parser utilities (archived)',
		version='0.0',
		author='mesoic',
		author_email="mesoic@protonmail.com",
		keywords='dataparser',
		license='MIT License',
		python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*',
		install_requires=[],
		packages=['dataparser', 'dataparser.parsers'],
		platforms="Linux, Windows, Mac",
		use_2to3=False,
		zip_safe=False,
)