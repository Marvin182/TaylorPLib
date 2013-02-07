# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_TaylorPLib', [dirname(__file__)])
        except ImportError:
            import _TaylorPLib
            return _TaylorPLib
        if fp is not None:
            try:
                _mod = imp.load_module('_TaylorPLib', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _TaylorPLib = swig_import_helper()
    del swig_import_helper
else:
    import _TaylorPLib
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


MAX_MESSAGE_SIZE = _TaylorPLib.MAX_MESSAGE_SIZE
class Polynomial(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Polynomial, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Polynomial, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _TaylorPLib.new_Polynomial(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _TaylorPLib.delete_Polynomial
    __del__ = lambda self : None;
    def order(self): return _TaylorPLib.Polynomial_order(self)
    def ncoeff(self): return _TaylorPLib.Polynomial_ncoeff(self)
    def __eq__(self, *args): return _TaylorPLib.Polynomial___eq__(self, *args)
    def __ne__(self, *args): return _TaylorPLib.Polynomial___ne__(self, *args)
    def __lt__(self, *args): return _TaylorPLib.Polynomial___lt__(self, *args)
    def __le__(self, *args): return _TaylorPLib.Polynomial___le__(self, *args)
    def __gt__(self, *args): return _TaylorPLib.Polynomial___gt__(self, *args)
    def __ge__(self, *args): return _TaylorPLib.Polynomial___ge__(self, *args)
    def __add__(self, *args): return _TaylorPLib.Polynomial___add__(self, *args)
    def __iadd__(self, *args): return _TaylorPLib.Polynomial___iadd__(self, *args)
    def __neg__(self): return _TaylorPLib.Polynomial___neg__(self)
    def __sub__(self, *args): return _TaylorPLib.Polynomial___sub__(self, *args)
    def __isub__(self, *args): return _TaylorPLib.Polynomial___isub__(self, *args)
    def __mul__(self, *args): return _TaylorPLib.Polynomial___mul__(self, *args)
    def __imul__(self, *args): return _TaylorPLib.Polynomial___imul__(self, *args)
    def __div__(self, *args): return _TaylorPLib.Polynomial___div__(self, *args)
    def __idiv__(self, *args): return _TaylorPLib.Polynomial___idiv__(self, *args)
    def sqr(self): return _TaylorPLib.Polynomial_sqr(self)
    def setSqr(self): return _TaylorPLib.Polynomial_setSqr(self)
    def sqrt(self): return _TaylorPLib.Polynomial_sqrt(self)
    def setSqrt(self): return _TaylorPLib.Polynomial_setSqrt(self)
    def _print(self, *args): return _TaylorPLib.Polynomial__print(self, *args)
    def eval(self, *args): return _TaylorPLib.Polynomial_eval(self, *args)
    def feval(self): return _TaylorPLib.Polynomial_feval(self)
    def shift(self): return _TaylorPLib.Polynomial_shift(self)
    def isConst(self, *args): return _TaylorPLib.Polynomial_isConst(self, *args)
    def isId(self, *args): return _TaylorPLib.Polynomial_isId(self, *args)
    def isZero(self, *args): return _TaylorPLib.Polynomial_isZero(self, *args)
    def set2Const(self, *args): return _TaylorPLib.Polynomial_set2Const(self, *args)
    def set2Id(self): return _TaylorPLib.Polynomial_set2Id(self)
    def set2Zero(self, *args): return _TaylorPLib.Polynomial_set2Zero(self, *args)
    def setCoeffs(self, *args): return _TaylorPLib.Polynomial_setCoeffs(self, *args)
Polynomial_swigregister = _TaylorPLib.Polynomial_swigregister
Polynomial_swigregister(Polynomial)

class Matrix(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Matrix, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Matrix, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _TaylorPLib.new_Matrix(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _TaylorPLib.delete_Matrix
    __del__ = lambda self : None;
    def nrows(self): return _TaylorPLib.Matrix_nrows(self)
    def ncols(self): return _TaylorPLib.Matrix_ncols(self)
    def dimT(self): return _TaylorPLib.Matrix_dimT(self)
    def get(self, *args): return _TaylorPLib.Matrix_get(self, *args)
    def __call__(self, *args): return _TaylorPLib.Matrix___call__(self, *args)
    def __eq__(self, *args): return _TaylorPLib.Matrix___eq__(self, *args)
    def __ne__(self, *args): return _TaylorPLib.Matrix___ne__(self, *args)
    def __add__(self, *args): return _TaylorPLib.Matrix___add__(self, *args)
    def __iadd__(self, *args): return _TaylorPLib.Matrix___iadd__(self, *args)
    def __sub__(self, *args): return _TaylorPLib.Matrix___sub__(self, *args)
    def __isub__(self, *args): return _TaylorPLib.Matrix___isub__(self, *args)
    def __neg__(self): return _TaylorPLib.Matrix___neg__(self)
    def __imul__(self, *args): return _TaylorPLib.Matrix___imul__(self, *args)
    def __mul__(self, *args): return _TaylorPLib.Matrix___mul__(self, *args)
    def mmCaABbC(self, *args): return _TaylorPLib.Matrix_mmCaABbC(self, *args)
    def bmmCaABbC(self, *args): return _TaylorPLib.Matrix_bmmCaABbC(self, *args)
    def mmCasABbC(self, *args): return _TaylorPLib.Matrix_mmCasABbC(self, *args)
    def mmCaAsBbC(self, *args): return _TaylorPLib.Matrix_mmCaAsBbC(self, *args)
    def mmCaAUTBPbC(self, *args): return _TaylorPLib.Matrix_mmCaAUTBPbC(self, *args)
    def mmCaAATbC(self, *args): return _TaylorPLib.Matrix_mmCaAATbC(self, *args)
    def mmCaATAbC(self, *args): return _TaylorPLib.Matrix_mmCaATAbC(self, *args)
    def mmCaATBbC(self, *args): return _TaylorPLib.Matrix_mmCaATBbC(self, *args)
    def mmCaATBPbC(self, *args): return _TaylorPLib.Matrix_mmCaATBPbC(self, *args)
    def mmCaABTbC(self, *args): return _TaylorPLib.Matrix_mmCaABTbC(self, *args)
    def bmmCaABTbC(self, *args): return _TaylorPLib.Matrix_bmmCaABTbC(self, *args)
    def mmCaIBbC(self, *args): return _TaylorPLib.Matrix_mmCaIBbC(self, *args)
    def mmCaAIbC(self, *args): return _TaylorPLib.Matrix_mmCaAIbC(self, *args)
    def cpermutem(self, *args): return _TaylorPLib.Matrix_cpermutem(self, *args)
    def rpermutem(self, *args): return _TaylorPLib.Matrix_rpermutem(self, *args)
    def transpose(self): return _TaylorPLib.Matrix_transpose(self)
    def asTranspose(self): return _TaylorPLib.Matrix_asTranspose(self)
    def shift(self): return _TaylorPLib.Matrix_shift(self)
    def isId(self): return _TaylorPLib.Matrix_isId(self)
    def isZero(self): return _TaylorPLib.Matrix_isZero(self)
    def set2Id(self, *args): return _TaylorPLib.Matrix_set2Id(self, *args)
    def set2IdFromIndices(self, *args): return _TaylorPLib.Matrix_set2IdFromIndices(self, *args)
    def set2Zero(self, *args): return _TaylorPLib.Matrix_set2Zero(self, *args)
    def set2ZeroFromIndices(self, *args): return _TaylorPLib.Matrix_set2ZeroFromIndices(self, *args)
    def set2Val(self, *args): return _TaylorPLib.Matrix_set2Val(self, *args)
    def set2ValFromIndices(self, *args): return _TaylorPLib.Matrix_set2ValFromIndices(self, *args)
    def _print(self, *args): return _TaylorPLib.Matrix__print(self, *args)
Matrix_swigregister = _TaylorPLib.Matrix_swigregister
Matrix_swigregister(Matrix)

class MathException(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MathException, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MathException, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _TaylorPLib.new_MathException(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _TaylorPLib.delete_MathException
    __del__ = lambda self : None;
    def what(self): return _TaylorPLib.MathException_what(self)
MathException_swigregister = _TaylorPLib.MathException_swigregister
MathException_swigregister(MathException)

# This file is compatible with both classic and new-style classes.


