from ..Pole import Pole
from ..SphToCart import SphToCart
from ..CartToSph import CartToSph
from ..ZeroTwoPi import ZeroTwoPi, calculate_trend_less_than_zero, calculate_trend_greater_than_twopi
import math
import numpy


class TestSphToCart:
	def test_sph_to_cart_and_cart_to_sphere_k_equals_1(self):
		cn,ce,cd =  SphToCart(trd=4.7123889,plg=1.57,k=1)
		assert cn == -0.99999968293183139
		assert cd == 0.0007963267107332633

		trd, plg = CartToSph(cn=cn,ce=ce,cd=cd)

		assert round(trd,2) == round(math.pi,2)
		assert round(plg,2) == 0 


class TestAzimuthConstrains:
	def test_calculate_trend_less_than_zero(self):
		twopi=numpy.dot(2.0,math.pi)
		trend = -2
		expected_trend = twopi + trend
		trend = calculate_trend_less_than_zero(trend, twopi)

		assert trend == expected_trend

	def test_calculate_trend_greater_than_twopi(self):
		twopi=numpy.dot(2.0,math.pi)
		trend = 7.2
		expected_trend = trend - twopi 
		trend = calculate_trend_greater_than_twopi(trend, twopi)
		assert trend == expected_trend


class TestZeroTwoPi:
	def setup(self):
		self.trend_less_than_zero = -1.7
		self.trend = 0 
		self.trend_over_two_pi = 7.2

	def test_trend_less_than_zero(self):
		trend = ZeroTwoPi(self.trend_less_than_zero)
		assert trend == 4.5831853071795861

	def test_trend_greater_than_or_equal_twopi(self):
		trend = ZeroTwoPi(self.trend_over_two_pi)
		assert trend == 0.91681469282041395

	def test_trend_is_within_0_to_twopi(self):
		trend = ZeroTwoPi(self.trend)
		assert round(trend,4) == 0  


class TestTrendAndPlunge:
	def test_k_equals_1(self):
		trend, plunge = Pole(trd=4.7123889,plg=1.57,k=1)
		assert round(trend,2) == round(math.pi,2)
		assert round(plunge,2) == 0 
 
	def test_k_equals_0(self):
		pass 

