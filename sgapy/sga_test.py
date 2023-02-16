# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 14:43:52 2023

@author: bmelosh
"""

import unittest
import numpy as np

from calcmv import calcmv
from sphtocart import sphtocart
from carttosph import carttosph
from zerotwopi import zerotwopi
from stcoordline import stcoordline
from pole import pole
from dircosaxes import dircosaxes
from cauchy import cauchy
from principalstress import principalstress


class SgaTestCase(unittest.TestCase):
    
    def test_calcmv(self):
        self.assertEqual(calcmv([1,2,3],[1,2,3]), (6.239302118134834, 0.9180016103888388, 0.7938543729768796, 2.155973186831644, 0.22838940273739197, 0.3424509207746073))
        self.assertEqual(calcmv([1,2, 0.1, 0.23],[1,3, 0.2, 0.1]), (6.237570311145323, 0.45021086631275453, 0.7359978441152887, 2.1306644186861927, 0.18273351363673201, 0.3133153233988649))
    
    def test_sphtocart(self):
        self.assertEqual(sphtocart(1, 1, 0), [0.2919265817264289, 0.4546487134128409, 0.8414709848078965])
        self.assertEqual(sphtocart(0,0,0), [1.0,0,0])
        self.assertEqual(sphtocart(0,0,1), [0,0,1.0])
    
    def test_carttosph(self):
        self.assertEqual(carttosph(1,1,1), [0.7853981633974483, 1.5707963267948966])
        
    def test_zerotwopi(self):
        self.assertEqual(zerotwopi(1.4), (1.4))
        self.assertEqual(zerotwopi(-7), -0.7168146928204138)
        self.assertEqual(zerotwopi(-2), 4.283185307179586)
        self.assertEqual(zerotwopi(7), 0.7168146928204138)
        
    def test_stcoordline(self):
        self.assertEqual(stcoordline(1,1,0), ([0.24689431284211627, 0.15852901519210347]))
        self.assertEqual(stcoordline(1,1,1), ([0.3350375824927938, 0.21512515777911212]))
        
    def test_pole(self):
        self.assertEqual(pole(2,2, 0), [3.5707963267948966, -0.42920367320510344])
        self.assertEqual(pole(2,2,1), [0.4292036732051034, -0.4292036732051034])
        
    def test_dircosaxes(self):
        outlist = [[-1.00000000e+00, -0.00000000e+00,  1.22464680e-16], [-9.18338009e-49,  1.00000000e+00, -7.49879891e-33], [-6.12323400e-17, -7.49879891e-33, -1.00000000e+00]]
        a = np.array(outlist)
        self.assertEqual(dircosaxes(0,np.pi,np.pi*-1).all(), a.all())
    
    def test_cauchy(self):
        stress = np.array(([[1,1,1],[1,2,3],[4,3,2]])).reshape(3,3)
        tx1 = 1
        px1 = 2
        tx3 =1
        strike = 3
        dip = 1
        T = np.array(([1,2,1]))
        pT = np.array(([1,1,1]))
        self.assertEqual(cauchy(stress, tx1, px1, tx3, strike, dip)[0].all(), T.all())
        self.assertEqual(cauchy(stress, tx1, px1, tx3, strike, dip)[1].all(), pT.all())
        
    def test_principalstress(self):
        stress = np.array(([[1,1,1],[1,2,1],[1,2,3]]))
        tx1 = 1
        px1 = 2
        tx3 = 1
        pstress = np.array(([[1, 1.83286821, 0.29868141],[2,2,2],[2,2,2]]))
        dCp = np.array(([[1,1,1],[1,1,1],[1,1,1]]))
        print(principalstress(stress, tx1, px1, tx3)[0][0], pstress[0])
        self.assertEqual(principalstress(stress, tx1, px1, tx3)[0][0], pstress[0])
        self.assertEqual(principalstress(stress, tx1, px1, tx3)[1], dCp

    

if __name__ == '__main__':
    unittest.main()