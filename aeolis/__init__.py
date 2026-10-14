'''This file is part of AeoLiS.
   
AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.
   
AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''


from __future__ import absolute_import

import os
import sys
import logging

# Threads for the compiled parallel kernels. By default numba starts one
# OpenMP worker per logical processor and the workers spin after every
# parallel region, competing with the single-threaded parts of a model step
# (measured: 30 % of the step on a hybrid 10-core laptop). Memory-bound
# kernels gain nothing from hyperthreads, so use half the logical processors
# unless AEOLIS_THREADS or NUMBA_NUM_THREADS says otherwise, and let idle
# workers sleep. All settings are only defaults: the environment wins.
os.environ.setdefault('OMP_WAIT_POLICY', 'PASSIVE')
os.environ.setdefault('KMP_BLOCKTIME', '0')
_threads = int(os.environ.get('AEOLIS_THREADS') or os.environ.get('NUMBA_NUM_THREADS')
               or max(1, (os.cpu_count() or 2) // 2))
if 'numba' not in sys.modules:
    # numba reads this at import and refuses a changed value afterwards
    os.environ.setdefault('NUMBA_NUM_THREADS', str(_threads))

class Logger(logging.getLoggerClass()):
    '''Custom logger class that supports multiline logging and raising
       exceptions that are logged with the full traceback'''

    
    def _log(self, lvl, msgs, *args, **kwargs):
        if not isinstance(msgs, list):
            msgs = [msgs]
        for msg in msgs:
            super(Logger, self)._log(lvl, msg, *args, **kwargs)

            
    def log_and_raise(self, msg, exc=Exception, *args, **kwargs):
        try:
            raise exc(msg)
        except:
            super(Logger, self).exception(msg, stack_info=True)
            raise


logging.setLoggerClass(Logger)

import aeolis.inout
import aeolis.model
import aeolis.wind
import aeolis.shear
import aeolis.bed
import aeolis.hydro
import aeolis.threshold
import aeolis.transport
import aeolis.netcdf
import aeolis.erosion
# import aeolis.vegetation

# numba may have been imported before aeolis, in which case the environment
# default above came too late; the runtime setter still applies.
try:
    import numba
    numba.set_num_threads(min(_threads, numba.config.NUMBA_NUM_THREADS))
except Exception:   # pragma: no cover
    pass

