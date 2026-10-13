"""Process kernels of AeoLiS that can be used stand-alone.

Each sub-package contains pure functions operating on plain numpy arrays
(no model state dictionary, no parameter dictionary). The classic
``module.function(s, p)`` entry points (e.g. :func:`aeolis.bed.update`)
remain available and call these kernels internally.
"""
