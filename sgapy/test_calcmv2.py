# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:42:07 2023

@author: bmelosh
"""

import unittest
import numpy as np

from calcmv import *


def test_calcmv(self):

		self.assertEqual(calcmv([1,2,3],[1,2,3]), [357.48567847616096, 52.59761786149341])
        self.assertEqual(calcmv([2],[3]), [112, 3])