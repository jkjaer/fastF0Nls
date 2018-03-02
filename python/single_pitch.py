import os
import ctypes
from ctypes import c_void_p, c_float, c_double, c_int
libpath = os.path.dirname(os.path.realpath(__file__))+"/../cpp/lib/single_pitch.so"
lib = ctypes.cdll.LoadLibrary(libpath)


lib.single_pitch_new.argtypes = [c_int, c_int, c_int, c_void_p];
lib.single_pitch_new.restype = c_void_p

lib.single_pitch_est.argtypes = [c_void_p, c_void_p, c_double, c_double];
lib.single_pitch_est.restype = c_double

lib.single_pitch_est_fast.argtypes = [c_void_p, c_void_p, c_double, c_double];
lib.single_pitch_est_fast.restype = c_double

lib.single_pitch_del.argtypes = [c_void_p];
lib.single_pitch_del.restype = None;

lib.single_pitch_model_order.argtypes = [c_void_p];
lib.single_pitch_model_order.restype = int;

class single_pitch(object):
    def __init__(self, nData, maxModelOrder, pitchBounds, nFftGrid=None):
        """
        Initialize the object
        
        Example:
        N = 500
        L = 15
        pitch_bounds = np.array([0.01, 0.45])
        sp = single_pitch(L, N, pitch_bounds)

        Uses default nFftGrid = 5NL. Otherwise set via parameter nFftGrid

        Pitch bound is in normalized frequency (1.0 = Nyquist frequency)

        """

        if nFftGrid == None:
            nFftGrid = 5*nData*maxModelOrder
            
        self.obj = lib.single_pitch_new(maxModelOrder, nFftGrid, nData,
                                        pitchBounds.ctypes.data)

    def est(self, data, lnBFZeroOrder=0.0, eps=1e-3, method=0):
        """
        Estimates based on double vector (data)

        The function returns the estimated frequency in radians per sample

        Addtional parameters
        lnBFZeroOrder: log Baysian factor for order zero. (default 0.0)
               May be increased if non desirable low power signals interfere 
               and should considered noise.
        eps:  Requested accuracy in radians per sample. (Defaut 1e-3)
              Higher accuracy (smaller eps) may require more solve time
              and vice-versa

        If method = 0 (default), then the algorithm uses the following steps

         1. Calculate the objective function for all candidate model
           order and on the Fourier grid
         2. Perform model order selection
         3. Refine the for the selected model order

        If method != 0, then the algorithm uses the more computational
        demanding steps (but possible more accurate)

         1. Calculate the objective function for all candidate model
           order and on the Fourier grid
         2. Refine the best estimates for each model order
         3. Perform model order selection

        """

        if method == 0:
            return lib.single_pitch_est_fast(self.obj, data.ctypes.data, 
                                        lnBFZeroOrder, eps)
        else:
            return lib.single_pitch_est(self.obj, data.ctypes.data, 
                                        lnBFZeroOrder, eps)

    def modelOrder(self):
        """
        Retuns the estimated model order for the lastest solve.

        """
        
        return lib.single_pitch_model_order(self.obj);
    
    def __del__(self):
        lib.single_pitch_del(self.obj)
            
