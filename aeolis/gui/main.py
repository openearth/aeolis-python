"""
Main entry point for AeoLiS GUI.

This module provides a simple launcher for the GUI that imports
from the legacy monolithic gui.py module. In the future, this will
be refactored to use the modular package structure.
"""

from tkinter import Tk

def launch_gui():
    """
    Launch the AeoLiS GUI application.
    
    Returns
    -------
    None
    """
    # Import here to avoid circular imports
    from aeolis.gui import AeolisGUI, dic
    
    # Create the main application window
    root = Tk()
    
    # Create an instance of the AeolisGUI class
    app = AeolisGUI(root, dic)
    
    # Bring window to front and give it focus
    root.lift()
    root.attributes('-topmost', True)
    root.after_idle(root.attributes, '-topmost', False)
    root.focus_force()
    
    # Start the Tkinter event loop
    root.mainloop()


if __name__ == "__main__":
    launch_gui()
