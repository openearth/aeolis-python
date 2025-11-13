"""
AeoLiS GUI - Graphical User Interface for AeoLiS Model Configuration and Visualization

This module provides a comprehensive GUI for:
- Reading and writing configuration files
- Visualizing domain setup (topography, vegetation, etc.)
- Plotting wind input data and wind roses
- Visualizing model output (2D and 1D transects)

This is the main application module that coordinates the GUI and visualizers.
"""

import aeolis
from tkinter import *
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
import traceback
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from aeolis.constants import DEFAULT_CONFIG

# Import utilities from gui package
from aeolis.gui.utils import (
    # Constants
    HILLSHADE_AZIMUTH, HILLSHADE_ALTITUDE, HILLSHADE_AMBIENT,
    TIME_UNIT_THRESHOLDS, TIME_UNIT_DIVISORS,
    OCEAN_DEPTH_THRESHOLD, OCEAN_DISTANCE_THRESHOLD, SUBSAMPLE_RATE_DIVISOR,
    NC_COORD_VARS, VARIABLE_LABELS, VARIABLE_TITLES,
    # Utility functions
    resolve_file_path, make_relative_path, determine_time_unit,
    extract_time_slice, apply_hillshade
)

# Import visualizers
from aeolis.gui.visualizers.domain import DomainVisualizer
from aeolis.gui.visualizers.wind import WindVisualizer
from aeolis.gui.visualizers.output_2d import Output2DVisualizer
from aeolis.gui.visualizers.output_1d import Output1DVisualizer

from windrose import WindroseAxes

# Initialize with default configuration
configfile = "No file selected"
dic = DEFAULT_CONFIG.copy()

class AeolisGUI:
    """
    Main GUI class for AeoLiS model configuration and visualization.
    
    This class provides a comprehensive graphical user interface for:
    - Reading and writing AeoLiS configuration files
    - Visualizing domain setup (topography, vegetation, grid parameters)
    - Displaying wind input data (time series and wind roses)
    - Visualizing model output in 2D and 1D (transects)
    - Interactive exploration of simulation results
    
    Parameters
    ----------
    root : Tk
        The root Tkinter window
    dic : dict
        Configuration dictionary containing model parameters
        
    Attributes
    ----------
    entries : dict
        Dictionary mapping field names to Entry widgets
    nc_data_cache : dict or None
        Cached NetCDF data for 2D visualization
    nc_data_cache_1d : dict or None
        Cached NetCDF data for 1D transect visualization
    wind_data_cache : dict or None
        Cached wind data for wind visualization
    """
    def __init__(self, root, dic):
        self.root = root
        self.dic = dic
        self.root.title('Aeolis')
        
        # Initialize attributes
        self.nc_data_cache = None
        self.overlay_veg_enabled = False
        self.entries = {}  # Initialize entries dictionary
        
        self.create_widgets()

    def get_config_dir(self):
        """Get the directory of the config file, or current directory if no file selected"""
        global configfile
        if configfile and configfile != "No file selected" and os.path.exists(configfile):
            return os.path.dirname(configfile)
        elif configfile and configfile != "No file selected" and os.path.dirname(configfile):
            # configfile might be a path even if file doesn't exist yet
            return os.path.dirname(configfile)
        else:
            return os.getcwd()

    def create_widgets(self):
        # Create a tab control widget
        tab_control = ttk.Notebook(self.root)
        # Create individual tabs
        self.create_input_file_tab(tab_control)
        self.create_domain_tab(tab_control)
        self.create_wind_input_tab(tab_control)
        self.create_plot_output_2d_tab(tab_control)
        self.create_plot_output_1d_tab(tab_control)
        # Pack the tab control to expand and fill the available space
        tab_control.pack(expand=1, fill='both')
        
        # Store reference to tab control for later use
        self.tab_control = tab_control
        
        # Bind tab change event to check if domain tab is selected
        tab_control.bind('<<NotebookTabChanged>>', self.on_tab_changed)

    def on_tab_changed(self, event):
        """Handle tab change event to auto-plot domain/wind when tab is selected"""
        # Get the currently selected tab index
        selected_tab = self.tab_control.index(self.tab_control.select())
        
        # Domain tab is at index 1 (0: Input file, 1: Domain, 2: Wind Input, 3: Timeframe, etc.)
        if selected_tab == 1:
            # Check if required files are defined
            xgrid = self.entries.get('xgrid_file', None)
            ygrid = self.entries.get('ygrid_file', None)
            bed = self.entries.get('bed_file', None)
            
            if xgrid and ygrid and bed:
                xgrid_val = xgrid.get().strip()
                ygrid_val = ygrid.get().strip()
                bed_val = bed.get().strip()
                
                # Only auto-plot if all three files are specified (not empty)
                if xgrid_val and ygrid_val and bed_val:
                    try:
                        # Check if domain_visualizer exists (tab may not be created yet)
                        if hasattr(self, 'domain_visualizer'):
                            self.domain_visualizer.plot_data('bed_file', 'Bed Elevation')
                    except Exception as e:
                        # Silently fail if plotting doesn't work (e.g., files don't exist)
                        pass
        
        # Wind Input tab is at index 2 (0: Input file, 1: Domain, 2: Wind Input, 3: Timeframe, etc.)
        elif selected_tab == 2:
            # Check if wind file is defined
            wind_file_entry = self.entries.get('wind_file', None)
            
            if wind_file_entry:
                wind_file_val = wind_file_entry.get().strip()
                
                # Only auto-plot if wind file is specified and hasn't been loaded yet
                if wind_file_val and not hasattr(self, 'wind_data_cache'):
                    try:
                        self.load_and_plot_wind()
                    except Exception as e:
                        # Silently fail if plotting doesn't work (e.g., file doesn't exist)
                        pass

    def create_label_entry(self, tab, text, value, row):
        # Create a label and entry widget for a given tab
        label = ttk.Label(tab, text=text)
        label.grid(row=row, column=0, sticky=W)
        entry = ttk.Entry(tab)
        # Convert None to empty string for cleaner display
        entry.insert(0, '' if value is None else str(value))
        entry.grid(row=row, column=1, sticky=W)
        return entry

    def create_input_file_tab(self, tab_control):
        # Create the 'Read/Write Inputfile' tab
        tab0 = ttk.Frame(tab_control)
        tab_control.add(tab0, text='Read/Write Inputfile')

        # Create frame for file operations
        file_ops_frame = ttk.LabelFrame(tab0, text="Configuration File", padding=20)
        file_ops_frame.pack(padx=20, pady=20, fill=BOTH, expand=True)

        # Current config file display
        current_file_label = ttk.Label(file_ops_frame, text="Current config file:")
        current_file_label.grid(row=0, column=0, sticky=W, pady=5)
        
        self.current_config_label = ttk.Label(file_ops_frame, text=configfile, 
                                             foreground='blue', wraplength=500)
        self.current_config_label.grid(row=0, column=1, columnspan=2, sticky=W, pady=5, padx=10)

        # Read new config file
        read_label = ttk.Label(file_ops_frame, text="Read new config file:")
        read_label.grid(row=1, column=0, sticky=W, pady=10)
        
        read_button = ttk.Button(file_ops_frame, text="Browse & Load Config", 
                                command=self.load_new_config)
        read_button.grid(row=1, column=1, sticky=W, pady=10, padx=10)

        # Separator
        separator = ttk.Separator(file_ops_frame, orient='horizontal')
        separator.grid(row=2, column=0, columnspan=3, sticky=(W, E), pady=20)

        # Save config file
        save_label = ttk.Label(file_ops_frame, text="Save config file as:")
        save_label.grid(row=3, column=0, sticky=W, pady=5)
        
        self.save_config_entry = ttk.Entry(file_ops_frame, width=40)
        self.save_config_entry.grid(row=3, column=1, sticky=W, pady=5, padx=10)
        
        save_browse_button = ttk.Button(file_ops_frame, text="Browse...", 
                                       command=self.browse_save_location)
        save_browse_button.grid(row=3, column=2, sticky=W, pady=5, padx=5)

        # Save button
        save_config_button = ttk.Button(file_ops_frame, text="Save Configuration", 
                                       command=self.save_config_file)
        save_config_button.grid(row=4, column=1, sticky=W, pady=10, padx=10)

    def create_domain_tab(self, tab_control):
        # Create the 'Domain' tab
        tab1 = ttk.Frame(tab_control)
        tab_control.add(tab1, text='Domain')

        # Create frame for Domain Parameters
        params_frame = ttk.LabelFrame(tab1, text="Domain Parameters", padding=10)
        params_frame.grid(row=0, column=0, padx=10, pady=10, sticky=(N, W, E))

        # Fields to be displayed in the 'Domain Parameters' frame
        fields = ['xgrid_file', 'ygrid_file', 'bed_file', 'ne_file', 'veg_file', 'threshold_file', 'fence_file', 'wave_mask', 'tide_mask', 'threshold_mask']
        # Create label and entry widgets for each field with browse buttons
        for i, field in enumerate(fields):
            label = ttk.Label(params_frame, text=f"{field}:")
            label.grid(row=i, column=0, sticky=W, pady=2)
            entry = ttk.Entry(params_frame, width=35)
            value = self.dic.get(field, '')
            # Convert None to empty string for cleaner display
            entry.insert(0, '' if value is None else str(value))
            entry.grid(row=i, column=1, sticky=W, pady=2, padx=(0, 5))
            self.entries[field] = entry
            
            # Add browse button for each field
            browse_btn = ttk.Button(params_frame, text="Browse...", 
                                   command=lambda e=entry: self.browse_file(e))
            browse_btn.grid(row=i, column=2, sticky=W, pady=2)

        # Create frame for Domain Visualization
        viz_frame = ttk.LabelFrame(tab1, text="Domain Visualization", padding=10)
        viz_frame.grid(row=0, column=1, padx=10, pady=10, sticky=(N, S, E, W))
        
        # Configure grid weights to allow expansion
        tab1.columnconfigure(1, weight=1)
        tab1.rowconfigure(0, weight=1)
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(7, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.colorbar = None  # Initialize colorbar attribute
        self.cbar_ax = None  # Initialize colorbar axes
        
        # Create canvas for the figure
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        
        # Initialize domain visualizer
        self.domain_visualizer = DomainVisualizer(
            self.ax, self.canvas, self.fig,
            lambda: self.entries,  # get_entries function
            self.get_config_dir    # get_config_dir function
        )

        # Create a frame for buttons
        button_frame = ttk.Frame(viz_frame)
        button_frame.pack(pady=5)

        # Create plot buttons - delegate to domain visualizer
        bed_button = ttk.Button(button_frame, text="Plot Bed", 
                               command=lambda: self.domain_visualizer.plot_data('bed_file', 'Bed Elevation'))
        bed_button.grid(row=0, column=0, padx=5)
        
        ne_button = ttk.Button(button_frame, text="Plot Ne", 
                              command=lambda: self.domain_visualizer.plot_data('ne_file', 'Ne'))
        ne_button.grid(row=0, column=1, padx=5)
        
        veg_button = ttk.Button(button_frame, text="Plot Vegetation", 
                               command=lambda: self.domain_visualizer.plot_data('veg_file', 'Vegetation'))
        veg_button.grid(row=0, column=2, padx=5)
        
        combined_button = ttk.Button(button_frame, text="Bed + Vegetation", 
                                    command=self.domain_visualizer.plot_combined)
        combined_button.grid(row=0, column=3, padx=5)
        
        # Add export button for domain visualization
        export_domain_button = ttk.Button(button_frame, text="Export PNG", 
                                         command=self.domain_visualizer.export_png)
        export_domain_button.grid(row=0, column=4, padx=5)

    def browse_file(self, entry_widget):
        """
        Open file dialog to select a file and update the entry widget.
        
        Parameters
        ----------
        entry_widget : Entry
            The Entry widget to update with the selected file path
        """
        # Get initial directory from config file location
        initial_dir = self.get_config_dir()
        
        # Get current value to determine initial directory
        current_value = entry_widget.get()
        if current_value:
            current_resolved = resolve_file_path(current_value, initial_dir)
            if current_resolved and os.path.exists(current_resolved):
                initial_dir = os.path.dirname(current_resolved)
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title="Select file",
            filetypes=(("Text files", "*.txt"), 
                      ("All files", "*.*"))
        )
        
        # Update entry if a file was selected
        if file_path:
            # Try to make path relative to config file directory for portability
            config_dir = self.get_config_dir()
            file_path = make_relative_path(file_path, config_dir)
            
            entry_widget.delete(0, END)
            entry_widget.insert(0, file_path)

    def browse_nc_file(self):
        """
        Open file dialog to select a NetCDF file.
        Automatically loads and plots the data after selection.
        """
        # Get initial directory from config file location
        initial_dir = self.get_config_dir()
        
        # Get current value to determine initial directory
        current_value = self.nc_file_entry.get()
        if current_value:
            current_resolved = resolve_file_path(current_value, initial_dir)
            if current_resolved and os.path.exists(current_resolved):
                initial_dir = os.path.dirname(current_resolved)
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title="Select NetCDF output file",
            filetypes=(("NetCDF files", "*.nc"), 
                      ("All files", "*.*"))
        )
        
        # Update entry if a file was selected
        if file_path:
            # Try to make path relative to config file directory for portability
            config_dir = self.get_config_dir()
            file_path = make_relative_path(file_path, config_dir)
            
            self.nc_file_entry.delete(0, END)
            self.nc_file_entry.insert(0, file_path)
            
            # Auto-load and plot the data using visualizer
            if hasattr(self, 'output_2d_visualizer'):
                self.output_2d_visualizer.load_and_plot()

    def load_new_config(self):
        """Load a new configuration file and update all fields"""
        global configfile
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=self.get_config_dir(),
            title="Select config file",
            filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
        )
        
        if file_path:
            try:
                # Read the new configuration file (parse_files=False to get file paths, not loaded arrays)
                self.dic = aeolis.inout.read_configfile(file_path, parse_files=False)
                configfile = file_path
                
                # Update the current file label
                self.current_config_label.config(text=configfile)
                
                # Update all entry fields with new values
                for field, entry in self.entries.items():
                    value = self.dic.get(field, '')
                    # Convert None to empty string, otherwise convert to string
                    value_str = '' if value is None else str(value)
                    entry.delete(0, END)
                    entry.insert(0, value_str)
                
                # Update NC file entry if it exists
                if hasattr(self, 'nc_file_entry'):
                    self.nc_file_entry.delete(0, END)
                
                # Clear wind data cache to force reload with new config
                if hasattr(self, 'wind_data_cache'):
                    delattr(self, 'wind_data_cache')
                
                # If on Wind Input tab and wind file is defined, reload and plot
                try:
                    selected_tab = self.tab_control.index(self.tab_control.select())
                    if selected_tab == 2:  # Wind Input tab
                        wind_file = self.wind_file_entry.get()
                        if wind_file and wind_file.strip():
                            self.load_and_plot_wind()
                except:
                    pass  # Silently fail if tabs not yet initialized
                
                messagebox.showinfo("Success", f"Configuration loaded from:\n{file_path}")
                
            except Exception as e:
                import traceback
                error_msg = f"Failed to load config file: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)

    def browse_save_location(self):
        """Browse for save location for config file"""
        # Open file dialog for saving
        file_path = filedialog.asksaveasfilename(
            initialdir=self.get_config_dir(),
            title="Save config file as",
            defaultextension=".txt",
            filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
        )
        
        if file_path:
            self.save_config_entry.delete(0, END)
            self.save_config_entry.insert(0, file_path)

    def save_config_file(self):
        """Save the current configuration to a file"""
        save_path = self.save_config_entry.get()
        
        if not save_path:
            messagebox.showwarning("Warning", "Please specify a file path to save the configuration.")
            return
        
        try:
            # Update dictionary with current entry values
            for field, entry in self.entries.items():
                self.dic[field] = entry.get()
            
            # Write the configuration file
            aeolis.inout.write_configfile(save_path, self.dic)
            
            messagebox.showinfo("Success", f"Configuration saved to:\n{save_path}")
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to save config file: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)

    def toggle_color_limits(self):
        """Enable or disable colorbar limit entries based on auto limits checkbox"""
        if self.auto_limits_var.get():
            self.vmin_entry.config(state='disabled')
            self.vmax_entry.config(state='disabled')
        else:
            self.vmin_entry.config(state='normal')
            self.vmax_entry.config(state='normal')

    def toggle_y_limits(self):
        """Enable or disable Y-axis limit entries based on auto limits checkbox"""
        if self.auto_ylimits_var.get():
            self.ymin_entry_1d.config(state='disabled')
            self.ymax_entry_1d.config(state='disabled')
        else:
            self.ymin_entry_1d.config(state='normal')
            self.ymax_entry_1d.config(state='normal')
        
        # Update plot if data is loaded
        if hasattr(self, 'nc_data_cache_1d') and self.nc_data_cache_1d is not None:
            self.update_1d_plot()

    def load_and_plot_wind(self):
        """
        Load and plot wind data using the wind visualizer.
        This is a wrapper method that delegates to the wind visualizer.
        """
        if hasattr(self, 'wind_visualizer'):
            self.wind_visualizer.load_and_plot()
    
    def browse_wind_file(self):
        """
        Open file dialog to select a wind file.
        Automatically loads and plots the wind data after selection.
        """
        # Get initial directory from config file location
        initial_dir = self.get_config_dir()
        
        # Get current value to determine initial directory
        current_value = self.wind_file_entry.get()
        if current_value:
            current_resolved = resolve_file_path(current_value, initial_dir)
            if current_resolved and os.path.exists(current_resolved):
                initial_dir = os.path.dirname(current_resolved)
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title="Select wind file",
            filetypes=(("Text files", "*.txt"), 
                      ("All files", "*.*"))
        )
        
        # Update entry if a file was selected
        if file_path:
            # Try to make path relative to config file directory for portability
            config_dir = self.get_config_dir()
            file_path = make_relative_path(file_path, config_dir)
            
            self.wind_file_entry.delete(0, END)
            self.wind_file_entry.insert(0, file_path)
            
            # Clear the cache to force reload of new file
            if hasattr(self, 'wind_data_cache'):
                delattr(self, 'wind_data_cache')
            
            # Auto-load and plot the data
            self.load_and_plot_wind()

    def create_wind_input_tab(self, tab_control):
        """Create the 'Wind Input' tab with wind data visualization"""
        tab_wind = ttk.Frame(tab_control)
        tab_control.add(tab_wind, text='Wind Input')

        # Create frame for wind file selection
        file_frame = ttk.LabelFrame(tab_wind, text="Wind File Selection", padding=10)
        file_frame.grid(row=0, column=0, padx=10, pady=10, sticky=(N, W, E))

        # Wind file selection
        wind_label = ttk.Label(file_frame, text="Wind file:")
        wind_label.grid(row=0, column=0, sticky=W, pady=2)
        
        # Create entry for wind file and store it in self.entries
        self.wind_file_entry = ttk.Entry(file_frame, width=35)
        wind_file_value = self.dic.get('wind_file', '')
        self.wind_file_entry.insert(0, '' if wind_file_value is None else str(wind_file_value))
        self.wind_file_entry.grid(row=0, column=1, sticky=W, pady=2, padx=(0, 5))
        self.entries['wind_file'] = self.wind_file_entry
        
        # Browse button for wind file
        wind_browse_btn = ttk.Button(file_frame, text="Browse...", 
                                     command=self.browse_wind_file)
        wind_browse_btn.grid(row=0, column=2, sticky=W, pady=2)

        # Create frame for time series plots
        timeseries_frame = ttk.LabelFrame(tab_wind, text="Wind Time Series", padding=10)
        timeseries_frame.grid(row=0, column=1, rowspan=2, padx=10, pady=10, sticky=(N, S, E, W))
        
        # Configure grid weights for expansion
        tab_wind.columnconfigure(1, weight=2)
        tab_wind.rowconfigure(0, weight=1)
        tab_wind.rowconfigure(1, weight=1)
        
        # Create matplotlib figure for time series (2 subplots stacked)
        self.wind_ts_fig = Figure(figsize=(7, 6), dpi=100)
        self.wind_ts_fig.subplots_adjust(hspace=0.35)
        self.wind_speed_ax = self.wind_ts_fig.add_subplot(211)
        self.wind_dir_ax = self.wind_ts_fig.add_subplot(212)
        
        # Create canvas for time series
        self.wind_ts_canvas = FigureCanvasTkAgg(self.wind_ts_fig, master=timeseries_frame)
        self.wind_ts_canvas.draw()
        self.wind_ts_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Create frame for windrose
        windrose_frame = ttk.LabelFrame(tab_wind, text="Wind Rose", padding=10)
        windrose_frame.grid(row=1, column=0, padx=10, pady=(0, 10), sticky=(N, S, E, W))
        
        # Create matplotlib figure for windrose
        self.windrose_fig = Figure(figsize=(5, 5), dpi=100)
        
        # Create canvas for windrose
        self.windrose_canvas = FigureCanvasTkAgg(self.windrose_fig, master=windrose_frame)
        self.windrose_canvas.draw()
        self.windrose_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        
        # Initialize wind visualizer
        self.wind_visualizer = WindVisualizer(
            self.wind_speed_ax, self.wind_dir_ax, self.wind_ts_canvas, self.wind_ts_fig,
            self.windrose_fig, self.windrose_canvas,
            lambda: self.wind_file_entry,  # get_wind_file function
            lambda: self.entries,           # get_entries function
            self.get_config_dir,            # get_config_dir function
            lambda: self.dic                # get_dic function
        )
        
        # Now add buttons that use the visualizer
        # Load button (forces reload by clearing cache)
        wind_load_btn = ttk.Button(file_frame, text="Load & Plot", 
                                   command=self.wind_visualizer.force_reload)
        wind_load_btn.grid(row=0, column=3, sticky=W, pady=2, padx=5)

        # Export buttons for wind plots
        export_label_wind = ttk.Label(file_frame, text="Export:")
        export_label_wind.grid(row=1, column=0, sticky=W, pady=5)
        
        export_button_frame_wind = ttk.Frame(file_frame)
        export_button_frame_wind.grid(row=1, column=1, columnspan=3, sticky=W, pady=5)
        
        export_wind_ts_btn = ttk.Button(export_button_frame_wind, text="Export Time Series PNG", 
                                        command=self.wind_visualizer.export_timeseries_png)
        export_wind_ts_btn.pack(side=LEFT, padx=5)
        
        export_windrose_btn = ttk.Button(export_button_frame_wind, text="Export Wind Rose PNG", 
                                         command=self.wind_visualizer.export_windrose_png)
        export_windrose_btn.pack(side=LEFT, padx=5)

    def create_plot_output_2d_tab(self, tab_control):
        # Create the 'Plot Output 2D' tab
        tab5 = ttk.Frame(tab_control)
        tab_control.add(tab5, text='Plot Output 2D')

        # Create frame for file selection
        file_frame = ttk.LabelFrame(tab5, text="Output File & Settings", padding=10)
        file_frame.grid(row=0, column=0, padx=10, pady=10, sticky=(N, W, E))

        # NC file selection
        nc_label = ttk.Label(file_frame, text="NetCDF file:")
        nc_label.grid(row=0, column=0, sticky=W, pady=2)
        self.nc_file_entry = ttk.Entry(file_frame, width=35)
        self.nc_file_entry.grid(row=0, column=1, sticky=W, pady=2, padx=(0, 5))
        
        # Browse button for NC file
        nc_browse_btn = ttk.Button(file_frame, text="Browse...", 
                                   command=lambda: self.browse_nc_file())
        nc_browse_btn.grid(row=0, column=2, sticky=W, pady=2)

        # Variable selection dropdown
        var_label_2d = ttk.Label(file_frame, text="Variable:")
        var_label_2d.grid(row=1, column=0, sticky=W, pady=2)
        
        # Initialize with empty list - will be populated when file is loaded
        self.variable_var_2d = StringVar(value='')
        self.variable_dropdown_2d = ttk.Combobox(file_frame, textvariable=self.variable_var_2d, 
                                        values=[], state='readonly', width=13)
        self.variable_dropdown_2d.grid(row=1, column=1, sticky=W, pady=2, padx=(0, 5))
        # Binding will be set after visualizer initialization
        self.variable_dropdown_2d_needs_binding = True

        # Colorbar limits
        vmin_label = ttk.Label(file_frame, text="Color min:")
        vmin_label.grid(row=2, column=0, sticky=W, pady=2)
        self.vmin_entry = ttk.Entry(file_frame, width=15, state='disabled')
        self.vmin_entry.grid(row=2, column=1, sticky=W, pady=2, padx=(0, 5))
        
        vmax_label = ttk.Label(file_frame, text="Color max:")
        vmax_label.grid(row=3, column=0, sticky=W, pady=2)
        self.vmax_entry = ttk.Entry(file_frame, width=15, state='disabled')
        self.vmax_entry.grid(row=3, column=1, sticky=W, pady=2, padx=(0, 5))
        
        # Auto limits checkbox
        self.auto_limits_var = BooleanVar(value=True)
        auto_limits_check = ttk.Checkbutton(file_frame, text="Auto limits", 
                                           variable=self.auto_limits_var,
                                           command=self.toggle_color_limits)
        auto_limits_check.grid(row=2, column=2, rowspan=2, sticky=W, pady=2)

        # Colormap selection
        cmap_label = ttk.Label(file_frame, text="Colormap:")
        cmap_label.grid(row=4, column=0, sticky=W, pady=2)
        
        # Available colormaps
        self.colormap_options = [
            'terrain',
            'viridis',
            'plasma',
            'inferno',
            'magma',
            'cividis',
            'jet',
            'rainbow',
            'turbo',
            'coolwarm',
            'seismic',
            'RdYlBu',
            'RdYlGn',
            'Spectral',
            'Greens',
            'Blues',
            'Reds',
            'gray',
            'hot',
            'cool'
        ]
        
        self.colormap_var = StringVar(value='terrain')
        colormap_dropdown = ttk.Combobox(file_frame, textvariable=self.colormap_var, 
                                        values=self.colormap_options, state='readonly', width=13)
        colormap_dropdown.grid(row=4, column=1, sticky=W, pady=2, padx=(0, 5))

        # Overlay vegetation checkbox
        self.overlay_veg_var = BooleanVar(value=False)
        overlay_veg_check = ttk.Checkbutton(file_frame, text="Overlay vegetation", 
                                           variable=self.overlay_veg_var)
        overlay_veg_check.grid(row=5, column=1, sticky=W, pady=2)

        # Export buttons
        export_label = ttk.Label(file_frame, text="Export:")
        export_label.grid(row=6, column=0, sticky=W, pady=5)
        
        export_button_frame = ttk.Frame(file_frame)
        export_button_frame.grid(row=6, column=1, columnspan=2, sticky=W, pady=5)
        
        export_png_btn = ttk.Button(export_button_frame, text="Export PNG", 
                                    command=lambda: self.output_2d_visualizer.export_png() if hasattr(self, 'output_2d_visualizer') else None)
        export_png_btn.pack(side=LEFT, padx=5)
        
        export_mp4_btn = ttk.Button(export_button_frame, text="Export Animation (MP4)", 
                                    command=lambda: self.output_2d_visualizer.export_animation_mp4() if hasattr(self, 'output_2d_visualizer') else None)
        export_mp4_btn.pack(side=LEFT, padx=5)

        # Create frame for visualization
        plot_frame = ttk.LabelFrame(tab5, text="Output Visualization", padding=10)
        plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky=(N, S, E, W))
        
        # Configure grid weights to allow expansion
        tab5.columnconfigure(1, weight=1)
        tab5.rowconfigure(0, weight=1)
        
        # Create matplotlib figure for output
        self.output_fig = Figure(figsize=(7, 6), dpi=100)
        self.output_ax = self.output_fig.add_subplot(111)
        self.output_colorbar = None
        self.output_cbar_ax = None
        
        # Create canvas for the output figure
        self.output_canvas = FigureCanvasTkAgg(self.output_fig, master=plot_frame)
        self.output_canvas.draw()
        self.output_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Create a frame for time slider
        slider_frame = ttk.Frame(plot_frame)
        slider_frame.pack(pady=5, fill=X, padx=10)
        
        # Time slider label
        self.time_label = ttk.Label(slider_frame, text="Time step: 0")
        self.time_label.pack(side=LEFT, padx=5)
        
        # Time slider
        self.time_slider = ttk.Scale(slider_frame, from_=0, to=0, orient=HORIZONTAL,
                                     command=self.update_time_step)
        self.time_slider.pack(side=LEFT, fill=X, expand=1, padx=5)
        self.time_slider.set(0)
        
        # Initialize 2D output visualizer (after all UI components are created)
        # Use a list to allow the visualizer to update the colorbar reference
        self.output_colorbar_ref = [self.output_colorbar]
        self.output_2d_visualizer = Output2DVisualizer(
            self.output_ax, self.output_canvas, self.output_fig,
            self.output_colorbar_ref, self.time_slider, self.time_label,
            self.variable_var_2d, self.colormap_var, self.auto_limits_var,
            self.vmin_entry, self.vmax_entry, self.overlay_veg_var,
            self.nc_file_entry, self.variable_dropdown_2d,
            self.get_config_dir, self.get_variable_label, self.get_variable_title
        )
        
        # Now bind the dropdown to use the visualizer
        self.variable_dropdown_2d.bind('<<ComboboxSelected>>', 
                                      lambda e: self.output_2d_visualizer.on_variable_changed(e))
        
        # Update time slider command to use visualizer
        self.time_slider.config(command=lambda v: self.output_2d_visualizer.update_plot())

    def create_plot_output_1d_tab(self, tab_control):
        # Create the 'Plot Output 1D' tab
        tab6 = ttk.Frame(tab_control)
        tab_control.add(tab6, text='Plot Output 1D')

        # Create frame for file selection
        file_frame_1d = ttk.LabelFrame(tab6, text="Output File & Transect Selection", padding=10)
        file_frame_1d.grid(row=0, column=0, padx=10, pady=10, sticky=(N, W, E))

        # NC file selection (shared with 2D plot)
        nc_label_1d = ttk.Label(file_frame_1d, text="NetCDF file:")
        nc_label_1d.grid(row=0, column=0, sticky=W, pady=2)
        self.nc_file_entry_1d = ttk.Entry(file_frame_1d, width=35)
        self.nc_file_entry_1d.grid(row=0, column=1, sticky=W, pady=2, padx=(0, 5))
        
        # Browse button for NC file
        nc_browse_btn_1d = ttk.Button(file_frame_1d, text="Browse...", 
                                       command=lambda: self.browse_nc_file_1d())
        nc_browse_btn_1d.grid(row=0, column=2, sticky=W, pady=2)

        # Variable selection dropdown
        var_label = ttk.Label(file_frame_1d, text="Variable:")
        var_label.grid(row=1, column=0, sticky=W, pady=2)
        
        # Initialize with empty list - will be populated when file is loaded
        self.variable_var_1d = StringVar(value='')
        self.variable_dropdown_1d = ttk.Combobox(file_frame_1d, textvariable=self.variable_var_1d, 
                                        values=[], state='readonly', width=13)
        self.variable_dropdown_1d.grid(row=1, column=1, sticky=W, pady=2, padx=(0, 5))
        self.variable_dropdown_1d.bind('<<ComboboxSelected>>', self.on_variable_changed)

        # Transect direction selection
        direction_label = ttk.Label(file_frame_1d, text="Transect direction:")
        direction_label.grid(row=2, column=0, sticky=W, pady=2)
        
        self.transect_direction_var = StringVar(value='cross-shore')
        direction_frame = ttk.Frame(file_frame_1d)
        direction_frame.grid(row=2, column=1, sticky=W, pady=2)
        
        cross_shore_radio = ttk.Radiobutton(direction_frame, text="Cross-shore (fix y-index)", 
                                            variable=self.transect_direction_var, value='cross-shore',
                                            command=self.update_transect_direction)
        cross_shore_radio.pack(side=LEFT, padx=5)
        
        along_shore_radio = ttk.Radiobutton(direction_frame, text="Along-shore (fix x-index)", 
                                            variable=self.transect_direction_var, value='along-shore',
                                            command=self.update_transect_direction)
        along_shore_radio.pack(side=LEFT, padx=5)

        # Transect position slider
        self.transect_label = ttk.Label(file_frame_1d, text="Y-index: 0")
        self.transect_label.grid(row=3, column=0, sticky=W, pady=2)
        
        self.transect_slider = ttk.Scale(file_frame_1d, from_=0, to=0, orient=HORIZONTAL,
                                         command=self.update_1d_transect_position)
        self.transect_slider.grid(row=3, column=1, sticky=(W, E), pady=2, padx=(0, 5))
        self.transect_slider.set(0)

        # Y-axis limits
        ymin_label = ttk.Label(file_frame_1d, text="Y-axis min:")
        ymin_label.grid(row=4, column=0, sticky=W, pady=2)
        self.ymin_entry_1d = ttk.Entry(file_frame_1d, width=15, state='disabled')
        self.ymin_entry_1d.grid(row=4, column=1, sticky=W, pady=2, padx=(0, 5))
        
        ymax_label = ttk.Label(file_frame_1d, text="Y-axis max:")
        ymax_label.grid(row=5, column=0, sticky=W, pady=2)
        self.ymax_entry_1d = ttk.Entry(file_frame_1d, width=15, state='disabled')
        self.ymax_entry_1d.grid(row=5, column=1, sticky=W, pady=2, padx=(0, 5))
        
        # Auto Y-axis limits checkbox
        self.auto_ylimits_var = BooleanVar(value=True)
        auto_ylimits_check = ttk.Checkbutton(file_frame_1d, text="Auto Y-axis limits", 
                                            variable=self.auto_ylimits_var,
                                            command=self.toggle_y_limits)
        auto_ylimits_check.grid(row=4, column=2, rowspan=2, sticky=W, pady=2)

        # Export buttons for 1D plots
        export_label_1d = ttk.Label(file_frame_1d, text="Export:")
        export_label_1d.grid(row=6, column=0, sticky=W, pady=5)
        
        export_button_frame_1d = ttk.Frame(file_frame_1d)
        export_button_frame_1d.grid(row=6, column=1, columnspan=2, sticky=W, pady=5)
        
        export_png_btn_1d = ttk.Button(export_button_frame_1d, text="Export PNG", 
                                       command=lambda: self.output_1d_visualizer.export_png() if hasattr(self, 'output_1d_visualizer') else None)
        export_png_btn_1d.pack(side=LEFT, padx=5)
        
        export_mp4_btn_1d = ttk.Button(export_button_frame_1d, text="Export Animation (MP4)", 
                                       command=lambda: self.output_1d_visualizer.export_animation_mp4() if hasattr(self, 'output_1d_visualizer') else None)
        export_mp4_btn_1d.pack(side=LEFT, padx=5)

        # Create frame for domain overview
        overview_frame = ttk.LabelFrame(tab6, text="Domain Overview", padding=10)
        overview_frame.grid(row=1, column=0, padx=10, pady=(0, 10), sticky=(N, S, E, W))
        
        # Create matplotlib figure for domain overview (smaller size)
        self.output_1d_overview_fig = Figure(figsize=(3.5, 3.5), dpi=80)
        self.output_1d_overview_fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.12)
        self.output_1d_overview_ax = self.output_1d_overview_fig.add_subplot(111)
        
        # Create canvas for the overview figure (centered, not expanded)
        self.output_1d_overview_canvas = FigureCanvasTkAgg(self.output_1d_overview_fig, master=overview_frame)
        self.output_1d_overview_canvas.draw()
        # Center the canvas both horizontally and vertically without expanding to fill
        canvas_widget = self.output_1d_overview_canvas.get_tk_widget()
        canvas_widget.pack(expand=True)

        # Create frame for transect visualization
        plot_frame_1d = ttk.LabelFrame(tab6, text="1D Transect Visualization", padding=10)
        plot_frame_1d.grid(row=0, column=1, rowspan=2, padx=10, pady=10, sticky=(N, S, E, W))
        
        # Configure grid weights to allow expansion
        tab6.columnconfigure(1, weight=1)
        tab6.rowconfigure(0, weight=1)
        tab6.rowconfigure(1, weight=1)
        
        # Create matplotlib figure for 1D transect output
        self.output_1d_fig = Figure(figsize=(7, 6), dpi=100)
        self.output_1d_ax = self.output_1d_fig.add_subplot(111)
        
        # Create canvas for the 1D output figure
        self.output_1d_canvas = FigureCanvasTkAgg(self.output_1d_fig, master=plot_frame_1d)
        self.output_1d_canvas.draw()
        self.output_1d_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Create a frame for time slider
        slider_frame_1d = ttk.Frame(plot_frame_1d)
        slider_frame_1d.pack(pady=5, fill=X, padx=10)
        
        # Time slider label
        self.time_label_1d = ttk.Label(slider_frame_1d, text="Time step: 0")
        self.time_label_1d.pack(side=LEFT, padx=5)
        
        # Time slider
        self.time_slider_1d = ttk.Scale(slider_frame_1d, from_=0, to=0, orient=HORIZONTAL,
                                        command=self.update_1d_time_step)
        self.time_slider_1d.pack(side=LEFT, fill=X, expand=1, padx=5)
        self.time_slider_1d.set(0)
        
        # Initialize 1D output visualizer (after all UI components are created)
        self.output_1d_visualizer = Output1DVisualizer(
            self.output_1d_ax, self.output_1d_overview_ax,
            self.output_1d_canvas, self.output_1d_fig,
            self.time_slider_1d, self.time_label_1d,
            self.transect_slider, self.transect_label,
            self.variable_var_1d, self.transect_direction_var,
            self.nc_file_entry_1d, self.variable_dropdown_1d,
            self.output_1d_overview_canvas,
            self.get_config_dir, self.get_variable_label, self.get_variable_title
        )
        
        # Update slider commands to use visualizer
        self.transect_slider.config(command=self.output_1d_visualizer.update_transect_position)
        self.time_slider_1d.config(command=self.output_1d_visualizer.update_time_step)
        
        # Update dropdown binding to use visualizer
        self.variable_dropdown_1d.unbind('<<ComboboxSelected>>')
        self.variable_dropdown_1d.bind('<<ComboboxSelected>>', 
                                      lambda e: self.output_1d_visualizer.update_plot())

    def browse_nc_file_1d(self):
        """
        Open file dialog to select a NetCDF file for 1D plotting.
        Automatically loads and plots the transect data after selection.
        """
        # Get initial directory from config file location
        initial_dir = self.get_config_dir()
        
        # Get current value to determine initial directory
        current_value = self.nc_file_entry_1d.get()
        if current_value:
            current_resolved = resolve_file_path(current_value, initial_dir)
            if current_resolved and os.path.exists(current_resolved):
                initial_dir = os.path.dirname(current_resolved)
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title="Select NetCDF output file",
            filetypes=(("NetCDF files", "*.nc"), 
                      ("All files", "*.*"))
        )
        
        # Update entry if a file was selected
        if file_path:
            # Try to make path relative to config file directory for portability
            config_dir = self.get_config_dir()
            file_path = make_relative_path(file_path, config_dir)
            
            self.nc_file_entry_1d.delete(0, END)
            self.nc_file_entry_1d.insert(0, file_path)
            
            # Auto-load and plot the data using visualizer
            if hasattr(self, 'output_1d_visualizer'):
                self.output_1d_visualizer.load_and_plot()

    def on_variable_changed(self, event):
        """Update plot when variable selection changes"""
        if hasattr(self, 'output_1d_visualizer'):
            self.output_1d_visualizer.update_plot()

    def update_transect_direction(self):
        """Update transect label and slider range when direction changes"""
        # Update plot if data is loaded
        if hasattr(self, 'output_1d_visualizer') and self.output_1d_visualizer.nc_data_cache_1d is not None:
            # Reload to reconfigure slider properly
            self.output_1d_visualizer.load_and_plot()

    def update_1d_transect_position(self, value):
        """Deprecated - now handled by visualizer"""
        pass

    def update_1d_time_step(self, value):
        """Deprecated - now handled by visualizer"""
        pass
    
    def update_1d_plot(self):
        """
        Update the 1D plot.
        This is a wrapper method that delegates to the 1D output visualizer.
        """
        if hasattr(self, 'output_1d_visualizer'):
            self.output_1d_visualizer.update_plot()

    def get_variable_label(self, var_name):
        """
        Get axis label for variable.
        
        Parameters
        ----------
        var_name : str
            Variable name
            
        Returns
        -------
        str
            Formatted label with units and fraction information if applicable
        """
        base_label = VARIABLE_LABELS.get(var_name, var_name)
        
        # Special cases that don't need fraction checking
        if var_name in ['zb+rhoveg', 'ustar quiver']:
            return base_label
        
        # Check if this variable has fractions dimension
        if hasattr(self, 'nc_data_cache') and self.nc_data_cache is not None:
            if var_name in self.nc_data_cache.get('vars', {}):
                var_data = self.nc_data_cache['vars'][var_name]
                if var_data.ndim == 4:
                    n_fractions = var_data.shape[3]
                    base_label += f' (averaged over {n_fractions} fractions)'
        
        return base_label

    def get_variable_title(self, var_name):
        """
        Get title for variable.
        
        Parameters
        ----------
        var_name : str
            Variable name
            
        Returns
        -------
        str
            Formatted title with fraction information if applicable
        """
        base_title = VARIABLE_TITLES.get(var_name, var_name)
        
        # Special cases that don't need fraction checking
        if var_name in ['zb+rhoveg', 'ustar quiver']:
            return base_title
        
        # Check if this variable has fractions dimension
        if hasattr(self, 'nc_data_cache') and self.nc_data_cache is not None:
            if var_name in self.nc_data_cache.get('vars', {}):
                var_data = self.nc_data_cache['vars'][var_name]
                if var_data.ndim == 4:
                    n_fractions = var_data.shape[3]
                    base_title += f' (averaged over {n_fractions} fractions)'
        
        return base_title

    def plot_nc_bed_level(self):
        """Plot bed level from NetCDF output file"""
        try:
            # Clear the previous plot
            self.output_ax.clear()
            
            # Get the NC file path
            nc_file = self.nc_file_entry.get()
            
            if not nc_file:
                messagebox.showwarning("Warning", "No NetCDF file specified!")
                return
            
            # Get the directory of the config file to resolve relative paths
            config_dir = self.get_config_dir()
            
            # Load the NC file
            if not os.path.isabs(nc_file):
                nc_file_path = os.path.join(config_dir, nc_file)
            else:
                nc_file_path = nc_file
                
            if not os.path.exists(nc_file_path):
                messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                return
            
            # Open NetCDF file and cache data
            with netCDF4.Dataset(nc_file_path, 'r') as nc:
                # Check if zb variable exists
                if 'zb' not in nc.variables:
                    available_vars = list(nc.variables.keys())
                    messagebox.showerror("Error", 
                        f"Variable 'zb' not found in NetCDF file.\n"
                        f"Available variables: {', '.join(available_vars)}")
                    return
                
                # Read bed level data (zb)
                zb_var = nc.variables['zb']
                
                # Check if time dimension exists
                if 'time' in zb_var.dimensions:
                    # Load all time steps
                    zb_data = zb_var[:]
                    n_times = zb_data.shape[0]
                else:
                    # Single time step
                    zb_data = zb_var[:, :]
                    zb_data = np.expand_dims(zb_data, axis=0)  # Add time dimension
                    n_times = 1
                
                # Try to get x and y coordinates
                x_data = None
                y_data = None
                
                if 'x' in nc.variables:
                    x_data = nc.variables['x'][:]
                if 'y' in nc.variables:
                    y_data = nc.variables['y'][:]
                
                # Create meshgrid if we have 1D coordinates
                if x_data is not None and y_data is not None:
                    if x_data.ndim == 1 and y_data.ndim == 1:
                        x_data, y_data = np.meshgrid(x_data, y_data)
                
                # Cache data for slider updates
                self.nc_data_cache = {
                    'zb': zb_data,
                    'x': x_data,
                    'y': y_data,
                    'n_times': n_times
                }
            
            # Configure the time slider
            if n_times > 1:
                self.time_slider.configure(from_=0, to=n_times-1)
                self.time_slider.set(n_times - 1)  # Start with last time step
            else:
                self.time_slider.configure(from_=0, to=0)
                self.time_slider.set(0)
            
            # Remember current output plot state
            self.output_plot_state = {
                'key': 'zb',
                'label': 'Elevation (m)',
                'title': 'Bed Elevation'
            }

            # Plot the initial (last) time step
            self.update_time_step(n_times - 1 if n_times > 1 else 0)
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to plot NetCDF bed level: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)  # Also print to console for debugging

    def plot_nc_wind(self):
        """Plot shear velocity (ustar) from NetCDF output file (uses 'ustar' or computes from 'ustars' and 'ustarn')."""
        try:
            # Clear the previous plot
            self.output_ax.clear()

            # Resolve file path
            nc_file = self.nc_file_entry.get()
            if not nc_file:
                messagebox.showwarning("Warning", "No NetCDF file specified!")
                return
            config_dir = self.get_config_dir()
            nc_file_path = os.path.join(config_dir, nc_file) if not os.path.isabs(nc_file) else nc_file
            if not os.path.exists(nc_file_path):
                messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                return

            with netCDF4.Dataset(nc_file_path, 'r') as nc:
                vars_available = set(nc.variables.keys())

                ustar_data = None
                ustars_data = None
                ustarn_data = None
                # Prefer magnitude if available
                if 'ustar' in vars_available:
                    ustar_var = nc.variables['ustar']
                    if 'time' in ustar_var.dimensions:
                        ustar_data = ustar_var[:]
                    else:
                        ustar_data = ustar_var[:, :]
                        ustar_data = np.expand_dims(ustar_data, axis=0)
                else:
                    # Try compute magnitude from components
                    if 'ustars' in vars_available and 'ustarn' in vars_available:
                        ustars_var = nc.variables['ustars']
                        ustarn_var = nc.variables['ustarn']
                        if 'time' in ustars_var.dimensions:
                            ustars_data = ustars_var[:]
                            ustarn_data = ustarn_var[:]
                        else:
                            ustars_data = np.expand_dims(ustars_var[:, :], axis=0)
                            ustarn_data = np.expand_dims(ustarn_var[:, :], axis=0)
                        ustar_data = np.sqrt(ustars_data**2 + ustarn_data**2)
                    else:
                        messagebox.showerror(
                            "Error",
                            "No shear velocity variables found in NetCDF file.\n"
                            "Expected 'ustar' or both 'ustars' and 'ustarn'.\n"
                            f"Available: {', '.join(sorted(vars_available))}"
                        )
                        return
                
                # If we have magnitude but not components, try loading components separately for quiver
                if ustar_data is not None and ustars_data is None:
                    if 'ustars' in vars_available and 'ustarn' in vars_available:
                        ustars_var = nc.variables['ustars']
                        ustarn_var = nc.variables['ustarn']
                        if 'time' in ustars_var.dimensions:
                            ustars_data = ustars_var[:]
                            ustarn_data = ustarn_var[:]
                        else:
                            ustars_data = np.expand_dims(ustars_var[:, :], axis=0)
                            ustarn_data = np.expand_dims(ustarn_var[:, :], axis=0)

                # Get coordinates
                x_data = nc.variables['x'][:] if 'x' in vars_available else None
                y_data = nc.variables['y'][:] if 'y' in vars_available else None
                if x_data is not None and y_data is not None:
                    if x_data.ndim == 1 and y_data.ndim == 1:
                        x_data, y_data = np.meshgrid(x_data, y_data)

                n_times = ustar_data.shape[0]

                # Initialize or update cache; keep existing cached fields
                if self.nc_data_cache is None:
                    self.nc_data_cache = {}
                cache_update = {
                    'ustar': ustar_data,
                    'x': x_data,
                    'y': y_data,
                    'n_times': n_times
                }
                # Add vector components if available
                if ustars_data is not None and ustarn_data is not None:
                    cache_update['ustars'] = ustars_data
                    cache_update['ustarn'] = ustarn_data
                self.nc_data_cache.update(cache_update)

            # Configure slider range
            if n_times > 1:
                self.time_slider.configure(from_=0, to=n_times-1)
                self.time_slider.set(n_times - 1)
            else:
                self.time_slider.configure(from_=0, to=0)
                self.time_slider.set(0)

            # Set plot state for shear velocity
            self.output_plot_state = {
                'key': 'ustar',
                'label': 'Shear velocity (m/s)',
                'title': 'Shear Velocity (ustar)'
            }

            # Render
            self.update_time_step(n_times - 1 if n_times > 1 else 0)

        except Exception as e:
            import traceback
            error_msg = f"Failed to plot NetCDF shear velocity: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)

    def apply_color_limits(self):
        """Re-plot with updated colorbar limits"""
        if self.nc_data_cache is not None:
            # Get current slider value and update the plot
            current_time = int(self.time_slider.get())
            self.update_time_step(current_time)
    
    def update_time_step(self, value):
        """
        Update the 2D plot based on the time slider value.
        This is a wrapper method that delegates to the 2D output visualizer.
        """
        if hasattr(self, 'output_2d_visualizer'):
            # Set the slider to the specified value
            self.time_slider.set(value)
            # Update the plot via the visualizer
            self.output_2d_visualizer.update_plot()

    def enable_overlay_vegetation(self):
        """Enable vegetation overlay in the output plot and load vegetation data if needed"""
        # Ensure bed data is loaded and slider configured
        if self.nc_data_cache is None:
            self.plot_nc_bed_level()
            if self.nc_data_cache is None:
                return

        # Load vegetation data into cache if not present
        if 'veg' not in self.nc_data_cache:
            try:
                # Resolve file path
                nc_file = self.nc_file_entry.get()
                if not nc_file:
                    messagebox.showwarning("Warning", "No NetCDF file specified!")
                    return
                config_dir = self.get_config_dir()
                nc_file_path = os.path.join(config_dir, nc_file) if not os.path.isabs(nc_file) else nc_file
                if not os.path.exists(nc_file_path):
                    messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                    return

                # Try common vegetation variable names
                veg_candidates = ['rhoveg', 'vegetated', 'hveg', 'vegfac']
                with netCDF4.Dataset(nc_file_path, 'r') as nc:
                    available = set(nc.variables.keys())
                    veg_name = next((v for v in veg_candidates if v in available), None)
                    if veg_name is None:
                        messagebox.showerror(
                            "Error",
                            "No vegetation variable found in NetCDF file.\n"
                            f"Tried: {', '.join(veg_candidates)}\n"
                            f"Available: {', '.join(sorted(available))}"
                        )
                        return
                    veg_var = nc.variables[veg_name]
                    # Read entire time series if time dimension exists
                    if 'time' in veg_var.dimensions:
                        veg_data = veg_var[:]
                    else:
                        veg_data = veg_var[:, :]

                # Cache vegetation data and name
                self.nc_data_cache['veg'] = veg_data
                self.nc_data_cache['veg_name'] = veg_name

            except Exception as e:
                import traceback
                error_msg = f"Failed to load vegetation from NetCDF: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
                return

        # Enable overlay and refresh current time step
        self.overlay_veg_enabled = True
        current_time = int(self.time_slider.get())
        self.update_time_step(current_time)

    def save(self):
        # Save the current entries to the configuration dictionary
        for field, entry in self.entries.items():
            self.dic[field] = entry.get()
        # Write the updated configuration to a new file
        aeolis.inout.write_configfile(configfile + '2', self.dic)
        print('Saved!')

if __name__ == "__main__":
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
