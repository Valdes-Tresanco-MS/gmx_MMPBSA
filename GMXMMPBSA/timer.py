"""
This is a module that contains the timer facility for gmx_MMPBSA. It's
useful for profiling the performance of the script.                   
"""

# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA                  #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

from time import time


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Timer(object):
    """ Timer class. It adds new timers and keeps track of how much time has been
        spent. """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self):
        """ Declares the global timer """
        self.timers = {'global': -time()}
        self.descriptions = {'global': 'Total time taken:'}
        self.active_timers = ['global']
        self.timer_names = ['global']
        self.units = 'sec.'

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def add_timer(self, timer_name, description):
        """ Add a new timer """
        # Check to make sure it doesn't already exist. Then add it if it doesn't
        if timer_name in self.timer_names:
            return

        self.timers[timer_name] = 0
        self.timer_names.append(timer_name)
        self.descriptions[timer_name] = description

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def start_timer(self, timer_name):
        """ Start the specified timer """
        # Make sure the timer is not on already
        if timer_name in self.active_timers:
            return

        # Check to make sure we've added the timer or not. If not, then add it
        if not timer_name in self.timer_names:
            self.add_timer(timer_name, '%s timer' % timer_name)

        self.timers[timer_name] -= time()

        # This timer is now on
        self.active_timers.append(timer_name)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def stop_timer(self, timer_name):
        """ End the specified timer """
        # First make sure that it is a timer to begin with
        if not timer_name in self.timer_names:
            return

        # Next make sure that the timer is on in the first place
        if not timer_name in self.active_timers:
            return

        # Now if it's on, end it and remove it from the list of active timers
        self.timers[timer_name] += time()
        self.active_timers.pop(self.active_timers.index(timer_name))

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def end_all(self):
        """ End all of the timers """
        while len(self.active_timers) > 0:
            self.stop_timer(self.active_timers[0])

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def done(self):
        """ Tell Timer we are done, and we can convert units to best
            human-readable option
        """

        # Make sure all timers are ended
        self.end_all()

        tfactor = 1
        # Now test the magnitude of the global timer so we can decide what the
        # reported units should be
        if self.timers['global'] > 60 * 60 * 48:
            self.units = 'days'
            tfactor = 60 * 60 * 24
        elif self.timers['global'] > 60 * 60:
            self.units = 'hr.'
            tfactor = 60 * 60
        elif self.timers['global'] > 60:
            self.units = 'min.'
            tfactor = 60

        for timer in self.timer_names:
            self.timers[timer] /= tfactor

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def print_(self, timer, outfile, newline=True):
        """ Prints the value of the timer """
        outfile.write("%-40s %8.3f %s" % (self.descriptions[timer],
                                          self.timers[timer], self.units))
        if newline:
            from os import linesep as ls
            outfile.write(ls)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
