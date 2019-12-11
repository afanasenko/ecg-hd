# coding: utf-8

import time

class IntervalTimer:

    def __init__(self):
        self.reset()

    def reset(self):
        self.t = time.clock()

    def interval(self):
        ts = time.clock()
        interval = ts - self.t
        self.t = ts
        return interval