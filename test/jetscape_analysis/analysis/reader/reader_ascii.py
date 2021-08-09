#!/usr/bin/env python3

"""
  Ascii reader class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

import os
import sys
import numpy as np
import pandas as pd

# Event class
from event import event_ascii

# Base class
from reader import reader_base

################################################################
class ReaderAscii(reader_base.ReaderBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_file_hadrons='', **kwargs):
        super(ReaderAscii, self).__init__(**kwargs)
        
        self.event_list_hadrons = self.parse_event(input_file_hadrons)
        self.current_event = 0
        self.n_events = len(self.event_list_hadrons)
  
    # ---------------------------------------------------------------
    # Parse the file into a list of events, each consisting of a list of lines
    # ---------------------------------------------------------------
    def parse_event(self, input_file):

        event_list = []
        event = None
        with open(input_file, 'r') as f:
            for line in f.readlines():
                            
                # If a new event, write the previous event and then clear it
                if line.startswith('#'):
                    if 'JETSCAPE_FINAL_STATE' in line:
                        continue
                        
                    if event:
                        event_list.append(event)
                    event = []
                    
                else:
                    event.append(np.array(line.rstrip('\n').split(), dtype=float))
                                    
            # Write the last event
            event_list.append(event)
        
        return event_list
        
    # ---------------------------------------------------------------
    # Get next event
    # Return event if successful, False if unsuccessful
    # ---------------------------------------------------------------
    def next_event(self):

        if self.current_event < self.n_events:
            self.current_event += 1
            event_hadrons = self.event_list_hadrons[self.current_event-1]
            return event_ascii.EventAscii(event_hadrons)
        else:
            sys.exit('Current event {} greater than total n_events {}'.format(self.current_event,
                                                                              self.n_events))
