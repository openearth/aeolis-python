"""
GUI Tabs package for AeoLiS GUI.

This package contains specialized tab modules for different types of data:
- domain: Domain setup visualization (bed, vegetation, etc.)
- wind: Wind input visualization (time series, wind roses)
- output_2d: 2D output visualization
- output_1d: 1D transect visualization
"""

from aeolis.gui.gui_tabs.domain import DomainVisualizer
from aeolis.gui.gui_tabs.wind import WindVisualizer
from aeolis.gui.gui_tabs.output_2d import Output2DVisualizer
from aeolis.gui.gui_tabs.output_1d import Output1DVisualizer

__all__ = ['DomainVisualizer', 'WindVisualizer', 'Output2DVisualizer', 'Output1DVisualizer']
