
import numpy as np

from numba import jit, njit, jitclass, float64
from numba.types import Array


@njit(float64(Array(float64, 1, 'C'), Array(float64, 1, 'C')))
def _free_func(x0, x1):
    return np.dot(x0, x1)


# In practice, MyClass is very complicated and I don't want to decorate it with @jitclass
# since this requires me to provide a full numba spec for each class member.
# @jitclass([('x_ref', Array(float64, 1, 'C'))])
class MyClass():

    def __init__(self, x_ref):
        self.x_ref = x_ref
        # self._unused = ['a']*5

    def call_member(self, x_other):
        return self._member_func(self.x_ref, x_other)

    def call_free(self, x_other):
        return _free_func(self.x_ref, x_other)

    @staticmethod
    @njit
    def _member_func(x0, x1):
        return np.dot(x0, x1)


def main():
    test_x0 = np.random.random(10)
    test_x1 = np.random.random(10)
    test_instance = MyClass(test_x0)
    result1 = test_instance.call_free(test_x1)
    print(result1)
    result2 = test_instance.call_member(test_x1)
    print(result2)


if __name__ == '__main__':
    main()
