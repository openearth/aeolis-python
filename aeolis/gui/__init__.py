"""
AeoLiS GUI Package - Modular GUI for AeoLiS Model

This package provides a modular graphical user interface for configuring
and visualizing AeoLiS aeolian sediment transport model results.

The main entry point is launch_gui() which creates and runs the GUI application.
"""

# Import from the application module within the gui package
from aeolis.gui.application import AeolisGUI, configfile, dic
from aeolis.gui.main import launch_gui

__all__ = ['launch_gui', 'AeolisGUI', 'configfile', 'dic']
