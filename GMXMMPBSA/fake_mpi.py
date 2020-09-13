"""
This is a module that provides the same attributes as mpi4py so that
we can reduce the number of #ifdef MPIs that are used. Ideally, the
only ones that will be needed after using this module are the ones to
make sure that this is loaded #ifndef MPI and related ones.

                           GPL LICENSE INFO

Copyright (C) 2009 - 2011 Dwight McGee, Billy Miller III, and Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MPI(object):
   """ Fake MPI class :) """

   class Communicator(object):

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      def Get_rank(self):
         """ The rank of a serial process is always 0 """
         return 0

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      def Get_size(self):
         """ The size of a serial version is always 1 """
         return 1

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      def bcast(self, inp, root=0):
         """ Mimics an MPI_Bcast """
         return inp

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      def Barrier(self):
         """ Mimics an MPI_Barrier """
         pass

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      def Abort(self, status=1):
         """ Just exits """
         from sys import exit
         exit(status)

   @staticmethod
   def Finalize():
       return

   COMM_WORLD = Communicator()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
