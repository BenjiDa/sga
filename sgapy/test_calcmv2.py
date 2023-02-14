# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:42:07 2023

@author: bmelosh
"""

import unittest
#import numpy as np

from calcmv import calcmv

class CalcMvTestCase(unittest.TestCase):
    
    def test_calcmv(self):
        self.assertEqual(calcmv([1,2,3],[1,2,3]), (6.239302118134834, 0.9180016103888388, 0.7938543729768796, 2.155973186831644, 0.22838940273739197, 0.3424509207746073))
        self.assertEqual(calcmv([1,2, 0.1, 0.23],[1,3, 0.2, 0.1]), (6.237570311145323, 0.45021086631275453, 0.7359978441152887, 2.1306644186861927, 0.18273351363673201, 0.3133153233988649))
        
        


if __name__ == '__main__':
    unittest.main()