#!/usr/bin/env python3

"""
  Reader base class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
from jetscape_analysis.base import common_base

################################################################
class ReaderBase(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_file="", **kwargs):
        super(ReaderBase, self).__init__(**kwargs)

    # ---------------------------------------------------------------
    # Generator (in pythonic sense) to loop over all events
    # ---------------------------------------------------------------
    def __call__(self, n_events):

        for i in range(0, n_events):
            yield self.next_event()
