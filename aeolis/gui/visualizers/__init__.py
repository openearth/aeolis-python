"""
Visualizers package for AeoLiS GUI.

This package contains specialized visualizer modules for different types of data:
- domain: Domain setup visualization (bed, vegetation, etc.)
- wind: Wind input visualization (time series, wind roses)
- output_2d: 2D output visualization
- output_1d: 1D transect visualization
"""

from aeolis.gui.visualizers.domain import DomainVisualizer
from aeolis.gui.visualizers.wind import WindVisualizer
from aeolis.gui.visualizers.output_2d import Output2DVisualizer
from aeolis.gui.visualizers.output_1d import Output1DVisualizer

__all__ = ['DomainVisualizer', 'WindVisualizer', 'Output2DVisualizer', 'Output1DVisualizer']
