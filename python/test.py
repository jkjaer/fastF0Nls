"""
Unittest of shared library provided through module single_pitch.py

"""
import unittest
import single_pitch
import numpy as np
import time

class TestModelOrderOne(unittest.TestCase):

    def test_standard(self):

        N = 500
        omega = 0.1
        sigma = 0.01
        L = 15
        pitch_bounds = np.array([0.01, 0.045])
        F = 5*L*N
        sp = single_pitch.single_pitch(L, F, N, pitch_bounds)
        lnBFZeroOrder = 10.0
        rpt = 3
    
        for k in range(rpt):
            data = np.sin(omega*np.arange(0, N)) + sigma*np.random.randn(N)    
            
            omega_h = sp.est(data, lnBFZeroOrder, 1e-6, method=1)
            ell_h = sp.modelOrder()
            
            self.assertTrue(np.allclose(omega, omega_h, rtol=2e-4, atol=2e-4))
            self.assertEqual(1, ell_h)
            
        del sp
        
    def test_fast(self):

        N = 500
        omega = 0.1
        sigma = 0.01
        L = 15
        pitch_bounds = np.array([0.01, 0.045])
        F = 5*L*N
        sp = single_pitch.single_pitch(L, F, N, pitch_bounds)
        lnBFZeroOrder = 10.0
        rpt = 3
    
        for k in range(rpt):
            data = np.sin(omega*np.arange(0, N)) + sigma*np.random.randn(N)    
            
            omega_h = sp.est(data, lnBFZeroOrder, 1e-3)
            ell_h = sp.modelOrder()
            
            self.assertTrue(np.allclose(omega, omega_h, rtol=1e-3, atol=1e-3))
            self.assertEqual(1, ell_h)
            
        del sp

class TestModelOrderTwo(unittest.TestCase):

    def test_standard(self):

        N = 500
        omega = 0.1
        sigma = 0.01
        L = 15
        pitch_bounds = np.array([0.01, 0.045])
        F = 5*L*N
        sp = single_pitch.single_pitch(L, F, N, pitch_bounds)
        lnBFZeroOrder = 10.0
        rpt = 3
    
        for k in range(rpt):
            data = np.sin(omega*np.arange(0, N)) + 0.9*np.sin(2*omega*np.arange(0, N) + 1.0) + sigma*np.random.randn(N)    
            
            omega_h = sp.est(data, lnBFZeroOrder, 1e-3, method=1)
            ell_h = sp.modelOrder()
            
            self.assertTrue(np.allclose(omega, omega_h, rtol=1e-3, atol=1e-3))
            self.assertEqual(2, ell_h)
            
        del sp

    def test_fast(self):

        N = 500
        omega = 0.1
        sigma = 0.01
        L = 15
        pitch_bounds = np.array([0.01, 0.045])
        F = 5*L*N
        sp = single_pitch.single_pitch(L, F, N, pitch_bounds)
        lnBFZeroOrder = 10.0
        rpt = 3
    
        for k in range(rpt):
            data = np.sin(omega*np.arange(0, N)) + 0.9*np.sin(2*omega*np.arange(0, N) + 1.0) + sigma*np.random.randn(N)    
            
            omega_h = sp.est(data, lnBFZeroOrder, 1e-3, method=0)
            ell_h = sp.modelOrder()
            
            self.assertTrue(np.allclose(omega, omega_h, rtol=1e-3, atol=1e-3))
            self.assertEqual(2, ell_h)
            
        del sp

if __name__ == '__main__':
    unittest.main()
