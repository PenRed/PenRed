import sys

_pyPenred = None
_checked = False

def get_pyPenred():
    global _pyPenred, _checked
    if _pyPenred is not None:
        return _pyPenred
    if _checked:
        return None
    _checked = True
    try:
        import pyPenred
        _pyPenred = pyPenred
        return _pyPenred
    except Exception:
        return None

def is_available():
    return get_pyPenred() is not None
