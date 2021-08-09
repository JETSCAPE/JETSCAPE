#!/usr/bin/env python3

"""
  Event class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

import pyhepmc_ng

# Base class
from event import event_base

################################################################
class EventAscii(event_base.EventBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, event_hadrons=None, **kwargs):
        super(EventAscii, self).__init__(**kwargs)
        
        self.event_hadrons = event_hadrons

    # ---------------------------------------------------------------
    # Get list of hadrons.
    # ---------------------------------------------------------------
    def hadrons(self, min_track_pt=0.):
    
        return self.particles(self.event_hadrons, min_track_pt=min_track_pt)

    # ---------------------------------------------------------------
    # Construct list of HepMC particles from event's list of particles
    # ---------------------------------------------------------------
    def particles(self, event, min_track_pt=0.):
    
        particles = []
        for particle in event:
            pid = int(particle[1])
            status = int(particle[2])
            e = particle[3]
            px = particle[4]
            py = particle[5]
            pz = particle[6]

            four_vector = pyhepmc_ng.FourVector(px, py, pz, e)
            particle = pyhepmc_ng.GenParticle(four_vector, pid, status)
        
            pt = particle.momentum.pt()
            if pid != 12 and pid != 14 and pid != 16: # Remove neutrinos
                if pt > min_track_pt:
                    particles.append(particle)
                
        return particles
