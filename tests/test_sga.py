import unittest
import numpy as np

from sga import *


class test_sga():


	def test_zerotwopi(self):

		self.assertEqual(zerotwopi(5), 5)
		self.assertEqual(zerotwopi(-1), 5.283185307179586)

	def test_sphtocart(self, trd,plg,k):

		self.assertEqual(sphtocart(1, 1, 1), [0.7080734182735712, -0.4546487134128409, 0.5403023058681398])
		self.assertEqual(sphtocart(1,1,0), [0.2919265817264289, 0.4546487134128409, 0.8414709848078965])

	def test_carttosph(self, cn,ce,cd):

		self.assertEqual(carttosph(0.2919265817264289, 0.4546487134128409, 0.8414709848078965), [1,1])
		self.assertEqual(carttosph(0.5, 0.7, 0.1), [0.9505468408120751, 0.1001674211615598])
		self.assertEqual(carttosph(3, 7, 1), [1.1659045405098132, 1.5707963267948966])

	def test_stcoordline(self, trd, plg, sttype):

		self.assertEqual(stcoord(1,1,0), [0.24689431284211627, 0.15852901519210347])
		self.assertEqual(stcoord(1,1,1), [0.3350375824927938, 0.21512515777911212])

	def test_pole(self, trd, plg, k):

		self.assertEqual(pole(1,1,0), [2.5707963267948966, 0.5707963267948966])
		self.assertEqual(pole(1,1,1), [5.71238898038469, 0.5707963267948967])
		self.assertEqual(pole(1,2,0), [2.5707963267948966, -0.42920367320510344])

	def test_dircosaxes(tX1,pX1,tX3):
		pass
		#self.assertEqual(dircosaxes(0, 1.570796, 1.570796), np.array([[ 3.26794897e-07,  0.00000000e+00,  1.00000000e+00], [ 1.00000000e+00, -3.26794897e-07, -3.26794897e-07],[ 3.26794897e-07,  1.00000000e+00, -1.06794904e-13]])
		#self.assertEqual(dircosaxes(0, 1.570796, -1.570796), np.array([[ 3.26794897e-07,  0.00000000e+00,  1.00000000e+00], [-1.00000000e+00, -3.26794897e-07,  3.26794897e-07],[ 3.26794897e-07, -1.00000000e+00, -1.06794904e-13]])



if __name__ == '__main__':
    unittest.main()