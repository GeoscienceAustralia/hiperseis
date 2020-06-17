#!/usr/bin/env python
# coding: utf-8
"""
Decorator to count number of times a function is called.
"""


def call_counter(func):
    """
    Decorator to count calls to a function. The number of calls can be queryied from func.counter.

    :param func: Function whose calls to count
    :return: func wrapper
    """
    def _inner(*args):
        _inner.counter += 1
        return func(*args)
    # end func
    _inner.counter = 0
    return _inner
# end func
