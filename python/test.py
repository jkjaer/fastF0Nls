"""
Unittest of shared library provided through module single_pitch.py

Can be executed as 

$ python test.py

or

$ nosetests

or

$ python -m unittest

Requires that the library ../cpp/lib/single_pitch.so is build

"""

import unittest
import single_pitch
import numpy as np
import time

def model(N, ell, omega, sigma):
    data = sigma*np.random.randn(N)
    for l in range(1, ell+1):
        data += np.sin(l*omega*np.arange(0, N))
    
    return data

class TestSinglePitch(unittest.TestCase):

    def setUp(self):
        self.N = 500
        self.omega = 0.1
        self.sigma = 0.01
        self.L = 15
        self.pitch_bounds = np.array([0.01, 0.45])
        self.sp = single_pitch.single_pitch(self.N, self.L, self.pitch_bounds)
        
    def test_order_one_standard(self):

        ell = 1
        rpt = 3
    
        for k in range(rpt):
            data = model(self.N, ell, self.omega, self.sigma)
            
            omega_h = self.sp.est(data, method=1)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, rtol=1e-3, atol=1e-3))
            self.assertEqual(ell, ell_h)
    
    def test_order_two_standard(self):

        ell = 2
        rpt = 3
        for k in range(rpt):
            data = model(self.N, ell, self.omega, self.sigma)
            
            omega_h = self.sp.est(data, method=1)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, rtol=2e-4, atol=2e-4))
            self.assertEqual(ell, ell_h)

    def test_order_one_fast(self):

        ell = 1
        rpt = 3
        for k in range(rpt):
            data = model(self.N, ell, self.omega, self.sigma)
            
            omega_h = self.sp.est(data)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, rtol=1e-3, atol=1e-3))
            self.assertEqual(ell, ell_h)

    def test_order_two_fast(self):

        ell = 2
        rpt = 3
        for k in range(rpt):
            data = model(self.N, ell, self.omega, self.sigma)
            
            omega_h = self.sp.est(data)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, rtol=1e-3, atol=1e-3))
            self.assertEqual(ell, ell_h)

    def test_order_two_fast_high_accuracy(self):

        ell = 2
        rpt = 3
        for k in range(rpt):
            data = model(self.N, ell, self.omega, 0.0)
            
            omega_h = self.sp.est(data, eps=1e-7)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, atol=1e-7))
            self.assertEqual(ell, ell_h)

    def test_order_two_fast_lnBfZero(self):

        ell = 2
        rpt = 3
        lnBFZeroOrder = 10.0
        for k in range(rpt):
            data = model(self.N, ell, self.omega, 0.0)
            
            omega_h = self.sp.est(data, lnBFZeroOrder=lnBFZeroOrder)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, atol=1e-3, rtol=1e-3))
            self.assertEqual(ell, ell_h)

    def test_order_ten_fast_lnBfZero(self):

        ell = 10
        rpt = 3
        for k in range(rpt):
            data = model(self.N, ell, self.omega, 0.0)
            
            omega_h = self.sp.est(data)
            ell_h = self.sp.modelOrder()
            
            self.assertTrue(np.allclose(self.omega, omega_h, atol=1e-3, rtol=1e-3))
            self.assertEqual(ell, ell_h)

    def tearDown(self):
        del self.sp


if __name__ == '__main__':
    unittest.main()
