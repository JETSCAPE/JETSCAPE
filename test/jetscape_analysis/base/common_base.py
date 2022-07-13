#!/usr/bin/env python3

"""
  Analysis base class.

  Author: Mateusz Ploskon
"""

################################################################
class CommonBase(object):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            self.__setattr__(key, value)

    # ---------------------------------------------------------------
    # Add an arbitrary attribute to the class
    # ---------------------------------------------------------------
    def set_attribute(self, **kwargs):
        for key, value in kwargs.items():
            self.__setattr__(key, value)

    # ---------------------------------------------------------------
    # Return formatted string of class members
    # ---------------------------------------------------------------
    def __str__(self):
        s = []
        variables = self.__dict__.keys()
        for v in variables:
            s.append("{} = {}".format(v, self.__dict__[v]))
        return "[i] {} with \n .  {}".format(self.__class__.__name__, "\n .  ".join(s))
