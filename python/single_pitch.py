import os
import ctypes
from ctypes import c_void_p, c_float, c_double, c_int
libpath = os.path.dirname(os.path.realpath(__file__))+"/../lib/single_pitch.so"
lib = ctypes.cdll.LoadLibrary(libpath)


lib.single_pitch_new.argtypes = [c_int, c_int, c_int, c_void_p];
lib.single_pitch_new.restype = c_void_p

lib.single_pitch_est.argtypes = [c_void_p, c_void_p, c_double];
lib.single_pitch_est.restype = c_double

lib.single_pitch_del.argtypes = [c_void_p];
lib.single_pitch_del.restype = None;

lib.single_pitch_model_order.argtypes = [c_void_p];
lib.single_pitch_model_order.restype = int;

class single_pitch(object):
    def __init__(self, maxModelOrder, nFftGrid, nData, pitchBounds):

        self.obj = lib.single_pitch_new(maxModelOrder, nFftGrid, nData,
                                        pitchBounds.ctypes.data)

    def est(self, data, lnBFZeroOrder=0.0):    
        return lib.single_pitch_est(self.obj, data.ctypes.data, 
                                    lnBFZeroOrder);

    def modelOrder(self):    
        return lib.single_pitch_model_order(self.obj);
    
    def __del__(self):
        lib.single_pitch_del(self.obj)


if __name__ == '__main__':
    import numpy as np
    import time

    N = 500
    omega = 0.1
    sigma = 0.01
    L = 15
    pitch_bounds = np.array([0.01, 0.045])
    F = 5*L*N
    sp = single_pitch(L, F, N, pitch_bounds)
    lnBFZeroOrder = 100.0
    rpt = 100

    t = time.time()
    
    for k in range(rpt):
        data = np.sin(omega*np.arange(0, N)) + sigma*np.random.randn(N)    
    
        omega_h = sp.est(data, lnBFZeroOrder)
        ell_h = sp.modelOrder()
        
        print("omega = {:8.4}, omega_h = {:8.4}, ell={:}".format(omega, omega_h, ell_h))

    t0 = time.time()-t

    print("Average time per solve = {:.6} [s]".format(t0/rpt))
    del sp
            
