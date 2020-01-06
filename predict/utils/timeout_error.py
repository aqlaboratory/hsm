import os
import errno, signal
from functools import wraps

"""
Code for TimeoutError/timeout decorator taken from:
    https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish

The use of this code is to cut off computations that take too long.
An example of the use here is:

'''
@timeout(s)
def function(x,y,z):
    # do stuff

try:
    function(x,y,z)
except TimeoutError as te:
    # do something if it doesn't complete in time

'''

The function will be cutoff after `s` seconds and a TimeoutError will be thrown. In 
the above example, we catch that error and do something else in this case.
"""
class TimeoutError(Exception):
    pass

def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator

