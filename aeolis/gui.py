import aeolis
from tkinter import *
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from aeolis.constants import DEFAULT_CONFIG

try:
    import netCDF4
    HAVE_NETCDF = True
except ImportError:
    HAVE_NETCDF = False

# Constants
COORD_VARS = {'x', 'y', 's', 'n', 'lat', 'lon', 'time', 'layers', 'fractions',
              'x_bounds', 'y_bounds', 'lat_bounds', 'lon_bounds', 'time_bounds', 'crs', 'nv', 'nv2'}

VEG_CANDIDATES = ['rhoveg', 'vegetated', 'hveg', 'vegfac']

COLORMAP_OPTIONS = [
    'terrain', 'viridis', 'plasma', 'inferno', 'magma', 'cividis',
    'jet', 'rainbow', 'turbo', 'coolwarm', 'seismic', 'RdYlBu',
    'RdYlGn', 'Spectral', 'Greens', 'Blues', 'Reds', 'gray', 'hot', 'cool'
]

VARIABLE_LABELS = {
    'zb': 'Elevation (m)',
    'zb+rhoveg': 'Vegetation-shaded Topography',
    'ustar': 'Shear Velocity (m/s)',
    'ustar quiver': 'Shear Velocity Vectors',
    'ustars': 'Shear Velocity S-component (m/s)',
    'ustarn': 'Shear Velocity N-component (m/s)',
    'zs': 'Surface Elevation (m)',
    'zsep': 'Separation Elevation (m)',
    'Ct': 'Sediment Concentration (kg/m²)',
    'Cu': 'Equilibrium Concentration (kg/m²)',
    'q': 'Sediment Flux (kg/m/s)',
    'qs': 'Sediment Flux S-component (kg/m/s)',
    'qn': 'Sediment Flux N-component (kg/m/s)',
    'pickup': 'Sediment Entrainment (kg/m²)',
    'uth': 'Threshold Shear Velocity (m/s)',
    'w': 'Fraction Weight (-)',
}

VARIABLE_TITLES = {
    'zb': 'Bed Elevation',
    'zb+rhoveg': 'Bed Elevation with Vegetation (Shaded)',
    'ustar': 'Shear Velocity',
    'ustar quiver': 'Shear Velocity Vector Field',
    'ustars': 'Shear Velocity (S-component)',
    'ustarn': 'Shear Velocity (N-component)',
    'zs': 'Surface Elevation',
    'zsep': 'Separation Elevation',
    'Ct': 'Sediment Concentration',
    'Cu': 'Equilibrium Concentration',
    'q': 'Sediment Flux',
    'qs': 'Sediment Flux (S-component)',
    'qn': 'Sediment Flux (N-component)',
    'pickup': 'Sediment Entrainment',
    'uth': 'Threshold Shear Velocity',
    'w': 'Fraction Weight',
}

# Hillshade parameters
HILLSHADE_AZIMUTH = 155.0
HILLSHADE_ALTITUDE = 5.0
HILLSHADE_AMBIENT = 0.35

# Color definitions for combined visualization
COLOR_SAND = np.array([1.0, 239.0/255.0, 213.0/255.0])
COLOR_DARKGREEN = np.array([34/255, 139/255, 34/255])
COLOR_OCEAN = np.array([70/255, 130/255, 180/255])

def apply_hillshade(z2d, x1d, y1d, az_deg=HILLSHADE_AZIMUTH, alt_deg=HILLSHADE_ALTITUDE):
    """
    Compute a simple hillshade (0–1) for 2D elevation array.
    Uses safe gradient computation and normalization.
    Adapted from Anim2D_ShadeVeg.py
    
    Args:
        z2d: 2D elevation array
        x1d: 1D x-coordinates
        y1d: 1D y-coordinates
        az_deg: Azimuth angle in degrees (default: 155.0)
        alt_deg: Altitude angle in degrees (default: 5.0)
    
    Returns:
        2D array with hillshade values in range [0, 1]
    """
    z = np.asarray(z2d, dtype=float)
    if z.ndim != 2:
        raise ValueError("apply_hillshade expects a 2D array")

    x1 = np.asarray(x1d).ravel()
    y1 = np.asarray(y1d).ravel()

    eps = 1e-8
    dx = np.mean(np.diff(x1)) if x1.size > 1 else 1.0
    dy = np.mean(np.diff(y1)) if y1.size > 1 else 1.0
    dx = 1.0 if abs(dx) < eps else dx
    dy = 1.0 if abs(dy) < eps else dy

    dz_dy, dz_dx = np.gradient(z, dy, dx)

    nx, ny, nz = -dz_dx, -dz_dy, np.ones_like(z)
    norm = np.sqrt(nx * nx + ny * ny + nz * nz)
    norm = np.where(norm < eps, eps, norm)
    nx, ny, nz = nx / norm, ny / norm, nz / norm

    az = math.radians(az_deg)
    alt = math.radians(alt_deg)
    lx = math.cos(alt) * math.cos(az)
    ly = math.cos(alt) * math.sin(az)
    lz = math.sin(alt)

    illum = np.clip(nx * lx + ny * ly + nz * lz, 0.0, 1.0)
    shaded = HILLSHADE_AMBIENT + (1.0 - HILLSHADE_AMBIENT) * illum
    return np.clip(shaded, 0.0, 1.0)

# Function to prompt the user to select a configuration file
def prompt_file():
    """
    Prompt user to select a configuration file.
    
    Returns:
        str: Path to selected file, or empty string if canceled
    """
    file_path = filedialog.askopenfilename(
        initialdir=os.getcwd(),
        title="Select config file",
        filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
    )
    return file_path

# Prompt the user to select a configuration file or use the default
selected_file = prompt_file()

# Read the configuration file into a dictionary
if selected_file:
    # User selected a file
    configfile = selected_file
    dic = aeolis.inout.read_configfile(configfile)
else:
    # User canceled - load empty fields with defaults
    configfile = "No file selected"
    # Use the default configuration from constants
    dic = DEFAULT_CONFIG.copy()


def resolve_relative_path(file_path, base_dir):
    """
    Resolve a file path relative to a base directory.
    
    Args:
        file_path: File path to resolve (can be absolute or relative)
        base_dir: Base directory for relative paths
    
    Returns:
        Absolute path to the file
    """
    if not file_path or not base_dir:
        return file_path
    
    if os.path.isabs(file_path):
        return file_path
    else:
        return os.path.join(base_dir, file_path)


def make_relative_path(file_path, base_dir):
    """
    Convert an absolute path to a relative path if beneficial.
    
    Args:
        file_path: Absolute file path to convert
        base_dir: Base directory to make path relative to
    
    Returns:
        Relative path if beneficial, otherwise absolute path
    """
    if not file_path or not base_dir:
        return file_path
    
    try:
        rel_path = os.path.relpath(file_path, base_dir)
        # Use relative path if it doesn't go up too many levels
        parent_dir = os.pardir + os.sep + os.pardir + os.sep
        if not rel_path.startswith(parent_dir):
            return rel_path
    except (ValueError, TypeError):
        # Different drives on Windows or invalid path, keep absolute path
        pass
    
    return file_path

class AeolisGUI:
    def __init__(self, root, dic):
        self.root = root
        self.dic = dic
        self.root.title('Aeolis')
        
        # Initialize attributes
        self.nc_data_cache = None
        self.overlay_veg_enabled = False
        
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

    def _browse_file_dialog(self, entry_widget, title, filetypes, initial_value=None):
        """
        Generic file browsing method to reduce code duplication.
        
        Args:
            entry_widget: The entry widget to update with selected file
            title: Dialog title
            filetypes: File types for the dialog
            initial_value: Current value in entry widget (if None, gets from widget)
        
        Returns:
            Selected file path or None if canceled
        """
        # Get initial directory from config file location
        initial_dir = self.get_config_dir()
        
        # Get current value to determine initial directory
        current_value = initial_value if initial_value is not None else entry_widget.get()
        if current_value:
            current_abs = resolve_relative_path(current_value, initial_dir)
            if os.path.exists(current_abs):
                initial_dir = os.path.dirname(current_abs)
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=initial_dir,
            title=title,
            filetypes=filetypes
        )
        
        # Update entry if a file was selected
        if file_path:
            # Try to make path relative to config file directory for portability
            rel_path = make_relative_path(file_path, self.get_config_dir())
            entry_widget.delete(0, END)
            entry_widget.insert(0, rel_path)
            return rel_path
        
        return None

    def _handle_colorbar(self, ax, im, label, colorbar_attr='colorbar'):
        """
        Handle colorbar creation and updates to avoid plot shrinking.
        
        Args:
            ax: Matplotlib axes
            im: Image object
            label: Colorbar label
            colorbar_attr: Attribute name for storing colorbar (default: 'colorbar')
        
        Returns:
            Colorbar object
        """
        colorbar = getattr(self, colorbar_attr, None)
        
        if colorbar is not None:
            try:
                # Update existing colorbar
                colorbar.update_normal(im)
                colorbar.set_label(label)
                return colorbar
            except (AttributeError, ValueError, RuntimeError):
                # If update fails, create new one
                pass
        
        # Create new colorbar
        colorbar = ax.figure.colorbar(im, ax=ax, label=label)
        setattr(self, colorbar_attr, colorbar)
        return colorbar

    def _remove_colorbar(self, colorbar_attr='colorbar'):
        """
        Safely remove a colorbar.
        
        Args:
            colorbar_attr: Attribute name for colorbar to remove
        """
        colorbar = getattr(self, colorbar_attr, None)
        if colorbar is not None:
            try:
                colorbar.remove()
            except (AttributeError, ValueError, RuntimeError):
                # If remove() fails, try removing from figure
                try:
                    if hasattr(colorbar, 'ax'):
                        colorbar.ax.figure.delaxes(colorbar.ax)
                except (AttributeError, ValueError):
                    pass
            setattr(self, colorbar_attr, None)

    def _resolve_file_path(self, file_path):
        """
        Resolve a file path relative to the config directory.
        
        Args:
            file_path: File path to resolve
        
        Returns:
            Absolute path to the file
        """
        if not file_path:
            return None
        
        config_dir = self.get_config_dir()
        return resolve_relative_path(file_path, config_dir)

    def _load_netcdf_variables(self, nc_file_path):
        """
        Load variables from NetCDF file and return organized data structure.
        
        Args:
            nc_file_path: Path to NetCDF file
        
        Returns:
            Dictionary containing:
                - 'vars': Dictionary of variable data arrays
                - 'x', 'y': Coordinate arrays (or None)
                - 's', 'n': Grid indices (or None)
                - 'n_times': Number of time steps
                - 'available_vars': List of available variable names
        """
        with netCDF4.Dataset(nc_file_path, 'r') as nc:
            # Get available variables
            available_vars = list(nc.variables.keys())
            
            # Try to get x and y coordinates
            x_data = nc.variables['x'][:] if 'x' in nc.variables else None
            y_data = nc.variables['y'][:] if 'y' in nc.variables else None
            
            # Get s and n coordinates (grid indices)
            s_data = nc.variables['s'][:] if 's' in nc.variables else None
            n_data = nc.variables['n'][:] if 'n' in nc.variables else None
            
            # Find all available 2D/3D variables (potential plot candidates)
            candidate_vars = []
            var_data_dict = {}
            n_times = 1
            
            for var_name in available_vars:
                if var_name in COORD_VARS:
                    continue
                
                var = nc.variables[var_name]
                
                # Check if time dimension exists
                if 'time' in var.dimensions:
                    # Load all time steps
                    var_data = var[:]
                    # Need at least 3 dimensions: (time, n, s)
                    if var_data.ndim < 3:
                        continue  # Skip variables without spatial dimensions
                    n_times = max(n_times, var_data.shape[0])
                else:
                    # Single time step - validate shape
                    if var.ndim < 2:
                        continue  # Skip variables without spatial dimensions
                    
                    if var.ndim == 2:
                        var_data = var[:, :]
                        var_data = np.expand_dims(var_data, axis=0)  # Add time dimension
                    elif var.ndim == 3:  # (n, s, fractions)
                        var_data = var[:, :, :]
                        var_data = np.expand_dims(var_data, axis=0)  # Add time dimension
                    else:
                        # Skip variables with more than 3 spatial dimensions
                        continue
                
                var_data_dict[var_name] = var_data
                candidate_vars.append(var_name)
            
            return {
                'vars': var_data_dict,
                'x': x_data,
                'y': y_data,
                's': s_data,
                'n': n_data,
                'n_times': n_times,
                'available_vars': candidate_vars
            }

    def create_widgets(self):
        # Create a tab control widget
        tab_control = ttk.Notebook(self.root)
        # Create individual tabs
        self.create_input_file_tab(tab_control)
        self.create_domain_tab(tab_control)
        self.create_timeframe_tab(tab_control)
        self.create_boundary_conditions_tab(tab_control)
        self.create_sediment_transport_tab(tab_control)
        self.create_plot_output_2d_tab(tab_control)
        self.create_plot_output_1d_tab(tab_control)
        # Pack the tab control to expand and fill the available space
        tab_control.pack(expand=1, fill='both')
        
        # Store reference to tab control for later use
        self.tab_control = tab_control
        
        # Bind tab change event to check if domain tab is selected
        tab_control.bind('<<NotebookTabChanged>>', self.on_tab_changed)

    def on_tab_changed(self, event):
        """Handle tab change event to auto-plot domain when tab is selected"""
        # Get the currently selected tab index
        selected_tab = self.tab_control.index(self.tab_control.select())
        
        # Domain tab is at index 1 (0: Input file, 1: Domain, 2: Timeframe, etc.)
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
                        self.plot_data('bed_file', 'Bed Elevation')
                    except Exception as e:
                        # Silently fail if plotting doesn't work (e.g., files don't exist)
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
        self.entries = {}
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

        # Create a frame for buttons
        button_frame = ttk.Frame(viz_frame)
        button_frame.pack(pady=5)

        # Create plot buttons
        bed_button = ttk.Button(button_frame, text="Plot Bed", command=lambda: self.plot_data('bed_file', 'Bed Elevation'))
        bed_button.grid(row=0, column=0, padx=5)
        
        ne_button = ttk.Button(button_frame, text="Plot Ne", command=lambda: self.plot_data('ne_file', 'Ne'))
        ne_button.grid(row=0, column=1, padx=5)
        
        veg_button = ttk.Button(button_frame, text="Plot Vegetation", command=lambda: self.plot_data('veg_file', 'Vegetation'))
        veg_button.grid(row=0, column=2, padx=5)
        
        combined_button = ttk.Button(button_frame, text="Bed + Vegetation", command=self.plot_combined)
        combined_button.grid(row=0, column=3, padx=5)

    def browse_file(self, entry_widget):
        """Open file dialog to select a file and update the entry widget"""
        self._browse_file_dialog(
            entry_widget,
            "Select file",
            (("Text files", "*.txt"), ("All files", "*.*"))
        )

    def browse_nc_file(self):
        """Open file dialog to select a NetCDF file"""
        result = self._browse_file_dialog(
            self.nc_file_entry,
            "Select NetCDF output file",
            (("NetCDF files", "*.nc"), ("All files", "*.*"))
        )
        
        # Auto-load and plot the data if a file was selected
        if result:
            self.plot_nc_2d()

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
                # Read the new configuration file
                self.dic = aeolis.inout.read_configfile(file_path)
                configfile = file_path
                
                # Update the current file label
                self.current_config_label.config(text=configfile)
                
                # Update all entry fields with new values
                for field, entry in self.entries.items():
                    entry.delete(0, END)
                    entry.insert(0, str(self.dic.get(field, '')))
                
                # Update NC file entry if it exists
                if hasattr(self, 'nc_file_entry'):
                    self.nc_file_entry.delete(0, END)
                
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

    def create_timeframe_tab(self, tab_control):
        # Create the 'Timeframe' tab
        tab2 = ttk.Frame(tab_control)
        tab_control.add(tab2, text='Timeframe')

        # Fields to be displayed in the 'Timeframe' tab
        fields = ['tstart', 'tstop', 'dt', 'restart', 'refdate']
        # Create label and entry widgets for each field
        self.entries.update({field: self.create_label_entry(tab2, f"{field}:", self.dic.get(field, ''), i) for i, field in enumerate(fields)})

    def create_boundary_conditions_tab(self, tab_control):
        # Create the 'Boundary Conditions' tab
        tab3 = ttk.Frame(tab_control)
        tab_control.add(tab3, text='Boundary Conditions')

        # Fields to be displayed in the 'Boundary Conditions' tab
        fields = ['boundary1', 'boundary2', 'boundary3']
        # Create label and entry widgets for each field
        self.entries.update({field: self.create_label_entry(tab3, f"{field}:", self.dic.get(field, ''), i) for i, field in enumerate(fields)})

    def create_sediment_transport_tab(self, tab_control):
        # Create the 'Sediment Transport' tab
        tab4 = ttk.Frame(tab_control)
        tab_control.add(tab4, text='Sediment Transport')

        # Create a 'Save' button
        save_button = ttk.Button(tab4, text='Save', command=self.save)
        save_button.pack()

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
        self.variable_dropdown_2d.bind('<<ComboboxSelected>>', self.on_variable_changed_2d)

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
        
        self.colormap_var = StringVar(value='terrain')
        colormap_dropdown = ttk.Combobox(file_frame, textvariable=self.colormap_var, 
                                        values=COLORMAP_OPTIONS, state='readonly', width=13)
        colormap_dropdown.grid(row=4, column=1, sticky=W, pady=2, padx=(0, 5))

        # Overlay vegetation checkbox
        self.overlay_veg_var = BooleanVar(value=False)
        overlay_veg_check = ttk.Checkbutton(file_frame, text="Overlay vegetation", 
                                           variable=self.overlay_veg_var)
        overlay_veg_check.grid(row=5, column=1, sticky=W, pady=2)

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

        # Create frame for visualization
        plot_frame_1d = ttk.LabelFrame(tab6, text="1D Transect Visualization", padding=10)
        plot_frame_1d.grid(row=0, column=1, padx=10, pady=10, sticky=(N, S, E, W))
        
        # Configure grid weights to allow expansion
        tab6.columnconfigure(1, weight=1)
        tab6.rowconfigure(0, weight=1)
        
        # Create matplotlib figure for 1D output
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

    def browse_nc_file_1d(self):
        """Open file dialog to select a NetCDF file for 1D plotting"""
        result = self._browse_file_dialog(
            self.nc_file_entry_1d,
            "Select NetCDF output file",
            (("NetCDF files", "*.nc"), ("All files", "*.*"))
        )
        
        # Auto-load and plot the data if a file was selected
        if result:
            self.plot_1d_transect()

    def on_variable_changed(self, event):
        """Update plot when variable selection changes"""
        if hasattr(self, 'nc_data_cache_1d') and self.nc_data_cache_1d is not None:
            self.update_1d_plot()

    def update_transect_direction(self):
        """Update transect label and slider range when direction changes"""
        # Update plot if data is loaded
        if hasattr(self, 'nc_data_cache_1d') and self.nc_data_cache_1d is not None:
            # Reconfigure slider range based on new direction
            first_var = list(self.nc_data_cache_1d['vars'].values())[0]
            
            if self.transect_direction_var.get() == 'cross-shore':
                # Fix y-index, vary along x (s dimension)
                max_idx = first_var.shape[1] - 1  # n dimension
                self.transect_slider.configure(from_=0, to=max_idx)
                # Set to middle or constrain current value
                current_val = int(self.transect_slider.get())
                if current_val > max_idx:
                    self.transect_slider.set(max_idx // 2)
                self.transect_label.config(text=f"Y-index: {int(self.transect_slider.get())}")
            else:
                # Fix x-index, vary along y (n dimension)
                max_idx = first_var.shape[2] - 1  # s dimension
                self.transect_slider.configure(from_=0, to=max_idx)
                # Set to middle or constrain current value
                current_val = int(self.transect_slider.get())
                if current_val > max_idx:
                    self.transect_slider.set(max_idx // 2)
                self.transect_label.config(text=f"X-index: {int(self.transect_slider.get())}")
            
            self.update_1d_plot()
        else:
            # Just update the label if no data loaded yet
            idx = int(self.transect_slider.get())
            if self.transect_direction_var.get() == 'cross-shore':
                self.transect_label.config(text=f"Y-index: {idx}")
            else:
                self.transect_label.config(text=f"X-index: {idx}")

    def update_1d_transect_position(self, value):
        """Update the transect position label"""
        idx = int(float(value))
        if self.transect_direction_var.get() == 'cross-shore':
            self.transect_label.config(text=f"Y-index: {idx}")
        else:
            self.transect_label.config(text=f"X-index: {idx}")
        
        # Update plot if data is loaded
        if hasattr(self, 'nc_data_cache_1d') and self.nc_data_cache_1d is not None:
            self.update_1d_plot()

    def update_1d_time_step(self, value):
        """Update the 1D plot based on the time slider value"""
        if not hasattr(self, 'nc_data_cache_1d') or self.nc_data_cache_1d is None:
            return
        
        # Get time index from slider
        time_idx = int(float(value))
        
        # Update label
        self.time_label_1d.config(text=f"Time step: {time_idx}")
        
        # Update plot
        self.update_1d_plot()

    def plot_1d_transect(self):
        """Load NetCDF file and plot 1D transect"""
        if not HAVE_NETCDF:
            messagebox.showerror("Error", "netCDF4 library is not available!")
            return
            
        try:
            # Get and resolve the NC file path
            nc_file = self.nc_file_entry_1d.get()
            if not nc_file:
                messagebox.showwarning("Warning", "No NetCDF file specified!")
                return
            
            nc_file_path = self._resolve_file_path(nc_file)
            if not os.path.exists(nc_file_path):
                messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                return
            
            # Load NetCDF data using helper method
            nc_data = self._load_netcdf_variables(nc_file_path)
            
            # Check if any variables were loaded
            if not nc_data['vars']:
                messagebox.showerror("Error", "No valid variables found in NetCDF file!")
                return
            
            # Update variable dropdown with available variables
            candidate_vars = sorted(nc_data['available_vars'])
            self.variable_dropdown_1d['values'] = candidate_vars
            
            # Set default to first variable (prefer 'zb' if available)
            if 'zb' in candidate_vars:
                self.variable_var_1d.set('zb')
            else:
                self.variable_var_1d.set(candidate_vars[0])
            
            # Cache data for slider updates
            self.nc_data_cache_1d = nc_data
            
            # Configure the time slider
            n_times = nc_data['n_times']
            if n_times > 1:
                self.time_slider_1d.configure(from_=0, to=n_times-1)
                self.time_slider_1d.set(n_times - 1)  # Start with last time step
            else:
                self.time_slider_1d.configure(from_=0, to=0)
                self.time_slider_1d.set(0)
            
            # Configure transect slider based on data shape
            first_var = list(nc_data['vars'].values())[0]
            if self.transect_direction_var.get() == 'cross-shore':
                # Fix y-index, vary along x (s dimension)
                max_idx = first_var.shape[1] - 1  # n dimension
                self.transect_slider.configure(from_=0, to=max_idx)
                self.transect_slider.set(max_idx // 2)  # Middle
            else:
                # Fix x-index, vary along y (n dimension)
                max_idx = first_var.shape[2] - 1  # s dimension
                self.transect_slider.configure(from_=0, to=max_idx)
                self.transect_slider.set(max_idx // 2)  # Middle
            
            # Plot the initial (last) time step
            self.update_1d_plot()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to plot 1D transect: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)  # Also print to console for debugging

    def update_1d_plot(self):
        """Update the 1D plot with current settings"""
        if not hasattr(self, 'nc_data_cache_1d') or self.nc_data_cache_1d is None:
            return
        
        try:
            # Clear the previous plot
            self.output_1d_ax.clear()
            
            # Get time index from slider
            time_idx = int(self.time_slider_1d.get())
            
            # Get transect index from slider
            transect_idx = int(self.transect_slider.get())
            
            # Get selected variable
            var_name = self.variable_var_1d.get()
            
            # Check if variable exists in cache
            if var_name not in self.nc_data_cache_1d['vars']:
                messagebox.showwarning("Warning", f"Variable '{var_name}' not found in NetCDF file!")
                return
            
            # Get the data
            var_data = self.nc_data_cache_1d['vars'][var_name]
            
            # Check if variable has fractions dimension (4D: time, n, s, fractions)
            has_fractions = var_data.ndim == 4
            
            # Extract transect based on direction
            if self.transect_direction_var.get() == 'cross-shore':
                # Fix y-index (n), vary along x (s)
                if has_fractions:
                    # Extract all fractions for this transect: (fractions,)
                    transect_data = var_data[time_idx, transect_idx, :, :]  # (s, fractions)
                    # Average or select first fraction
                    transect_data = transect_data.mean(axis=1)  # Average across fractions
                else:
                    transect_data = var_data[time_idx, transect_idx, :]
                
                # Get x-coordinates
                if self.nc_data_cache_1d['x'] is not None:
                    x_data = self.nc_data_cache_1d['x']
                    if x_data.ndim == 2:
                        x_coords = x_data[transect_idx, :]
                    else:
                        x_coords = x_data
                    xlabel = 'X (m)'
                elif self.nc_data_cache_1d['s'] is not None:
                    x_coords = self.nc_data_cache_1d['s']
                    xlabel = 'S-index'
                else:
                    x_coords = np.arange(len(transect_data))
                    xlabel = 'Grid Index'
            else:
                # Fix x-index (s), vary along y (n)
                if has_fractions:
                    # Extract all fractions for this transect: (fractions,)
                    transect_data = var_data[time_idx, :, transect_idx, :]  # (n, fractions)
                    # Average or select first fraction
                    transect_data = transect_data.mean(axis=1)  # Average across fractions
                else:
                    transect_data = var_data[time_idx, :, transect_idx]
                
                # Get y-coordinates
                if self.nc_data_cache_1d['y'] is not None:
                    y_data = self.nc_data_cache_1d['y']
                    if y_data.ndim == 2:
                        x_coords = y_data[:, transect_idx]
                    else:
                        x_coords = y_data
                    xlabel = 'Y (m)'
                elif self.nc_data_cache_1d['n'] is not None:
                    x_coords = self.nc_data_cache_1d['n']
                    xlabel = 'N-index'
                else:
                    x_coords = np.arange(len(transect_data))
                    xlabel = 'Grid Index'
            
            # Plot the transect
            self.output_1d_ax.plot(x_coords, transect_data, 'b-', linewidth=2)
            self.output_1d_ax.set_xlabel(xlabel)
            
            # Set ylabel based on variable (use same labels as 2D plots)
            ylabel = VARIABLE_LABELS.get(var_name, var_name)
            
            # Add indication if variable has fractions dimension
            if has_fractions:
                n_fractions = var_data.shape[3]
                ylabel += f' (averaged over {n_fractions} fractions)'
            
            self.output_1d_ax.set_ylabel(ylabel)
            
            # Set title
            direction = 'Cross-shore' if self.transect_direction_var.get() == 'cross-shore' else 'Along-shore'
            idx_label = 'Y' if self.transect_direction_var.get() == 'cross-shore' else 'X'
            
            # Get variable title (use same titles as 2D plots)
            var_title = VARIABLE_TITLES.get(var_name, var_name)
            if has_fractions:
                n_fractions = var_data.shape[3]
                var_title += f' (averaged over {n_fractions} fractions)'
            
            self.output_1d_ax.set_title(f'{direction} Transect: {var_title} ({idx_label}-index={transect_idx}, Time={time_idx})')
            
            # Apply Y-axis limits if specified
            if not self.auto_ylimits_var.get():
                try:
                    ymin_str = self.ymin_entry_1d.get().strip()
                    ymax_str = self.ymax_entry_1d.get().strip()
                    if ymin_str and ymax_str:
                        ymin = float(ymin_str)
                        ymax = float(ymax_str)
                        self.output_1d_ax.set_ylim(ymin, ymax)
                    elif ymin_str:
                        ymin = float(ymin_str)
                        self.output_1d_ax.set_ylim(bottom=ymin)
                    elif ymax_str:
                        ymax = float(ymax_str)
                        self.output_1d_ax.set_ylim(top=ymax)
                except ValueError:
                    pass  # Use auto limits if conversion fails
            
            # Add grid
            self.output_1d_ax.grid(True, alpha=0.3)
            
            # Redraw the canvas
            self.output_1d_canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to update 1D plot: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)  # Print to console for debugging

    def on_variable_changed_2d(self, event):
        """Update plot when variable selection changes in 2D tab"""
        if hasattr(self, 'nc_data_cache') and self.nc_data_cache is not None:
            self.update_2d_plot()

    def plot_nc_2d(self):
        """Load NetCDF file and plot 2D data"""
        if not HAVE_NETCDF:
            messagebox.showerror("Error", "netCDF4 library is not available!")
            return
            
        try:
            # Get and resolve the NC file path
            nc_file = self.nc_file_entry.get()
            if not nc_file:
                messagebox.showwarning("Warning", "No NetCDF file specified!")
                return
            
            nc_file_path = self._resolve_file_path(nc_file)
            if not os.path.exists(nc_file_path):
                messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                return
            
            # Load NetCDF data using helper method (loads only non-2D vars for 2D plotting)
            with netCDF4.Dataset(nc_file_path, 'r') as nc:
                # Get available variables
                available_vars = list(nc.variables.keys())
                
                # Get coordinates
                x_data = nc.variables['x'][:] if 'x' in nc.variables else None
                y_data = nc.variables['y'][:] if 'y' in nc.variables else None
                
                # Find all available 2D/3D variables
                candidate_vars = []
                var_data_dict = {}
                n_times = 1
                veg_data = None
                
                for var_name in available_vars:
                    if var_name in COORD_VARS:
                        continue
                    
                    var = nc.variables[var_name]
                    
                    # Check if time dimension exists
                    if 'time' in var.dimensions:
                        var_data = var[:]
                        # Need at least 3 dimensions: (time, n, s)
                        if var_data.ndim < 3:
                            continue
                        n_times = max(n_times, var_data.shape[0])
                    else:
                        # Single time step - need exactly 2 spatial dimensions
                        if var.ndim != 2:
                            continue
                        var_data = var[:, :]
                        var_data = np.expand_dims(var_data, axis=0)
                    
                    var_data_dict[var_name] = var_data
                    candidate_vars.append(var_name)
                
                # Load vegetation data if requested
                if self.overlay_veg_var.get():
                    for veg_name in VEG_CANDIDATES:
                        if veg_name in available_vars:
                            veg_var = nc.variables[veg_name]
                            if 'time' in veg_var.dimensions:
                                veg_data = veg_var[:]
                            else:
                                veg_data = veg_var[:, :]
                                veg_data = np.expand_dims(veg_data, axis=0)
                            break
                
                # Check if any variables were loaded
                if not var_data_dict:
                    messagebox.showerror("Error", "No valid variables found in NetCDF file!")
                    return
                
                # Add special combined options if components are available
                if 'zb' in var_data_dict and 'rhoveg' in var_data_dict:
                    candidate_vars.append('zb+rhoveg')
                if 'ustarn' in var_data_dict and 'ustars' in var_data_dict:
                    candidate_vars.append('ustar quiver')
                
                # Cache data for slider updates
                self.nc_data_cache = {
                    'vars': var_data_dict,
                    'x': x_data,
                    'y': y_data,
                    'n_times': n_times,
                    'available_vars': candidate_vars,
                    'veg': veg_data
                }
            
            # Update variable dropdown with available variables
            self.variable_dropdown_2d['values'] = sorted(candidate_vars)
            # Set default to 'zb' if available, otherwise first variable
            if 'zb' in candidate_vars:
                self.variable_var_2d.set('zb')
            else:
                self.variable_var_2d.set(sorted(candidate_vars)[0])
            
            # Configure the time slider
            if n_times > 1:
                self.time_slider.configure(from_=0, to=n_times-1)
                self.time_slider.set(n_times - 1)  # Start with last time step
            else:
                self.time_slider.configure(from_=0, to=0)
                self.time_slider.set(0)
            
            # Remember current output plot state
            self.output_plot_state = {
                'key': self.variable_var_2d.get(),
                'label': self.get_variable_label(self.variable_var_2d.get()),
                'title': self.get_variable_title(self.variable_var_2d.get())
            }

            # Plot the initial (last) time step
            self.update_2d_plot()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to plot 2D data: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)  # Also print to console for debugging

    def get_variable_label(self, var_name):
        """Get axis label for variable"""
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
        """Get title for variable"""
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

    def update_2d_plot(self):
        """Update the 2D plot with current settings"""
        if not hasattr(self, 'nc_data_cache') or self.nc_data_cache is None:
            return
        
        try:
            # Clear the previous plot
            self.output_ax.clear()
            
            # Get time index from slider
            time_idx = int(self.time_slider.get())
            
            # Get selected variable
            var_name = self.variable_var_2d.get()
            
            # Special handling for zb+rhoveg combined visualization
            if var_name == 'zb+rhoveg':
                self.render_zb_rhoveg_shaded(time_idx)
                return
            
            # Special handling for ustar quiver plot
            if var_name == 'ustar quiver':
                self.render_ustar_quiver(time_idx)
                return
            
            # Check if variable exists in cache
            if var_name not in self.nc_data_cache['vars']:
                messagebox.showwarning("Warning", f"Variable '{var_name}' not found in NetCDF file!")
                return
            
            # Get the data
            var_data = self.nc_data_cache['vars'][var_name]
            
            # Check if variable has fractions dimension (4D: time, n, s, fractions)
            if var_data.ndim == 4:
                # Average across fractions or select first fraction
                z_data = var_data[time_idx, :, :, :].mean(axis=2)  # Average across fractions
            else:
                z_data = var_data[time_idx, :, :]
            
            x_data = self.nc_data_cache['x']
            y_data = self.nc_data_cache['y']
            
            # Get colorbar limits
            vmin = None
            vmax = None
            if not self.auto_limits_var.get():
                try:
                    vmin_str = self.vmin_entry.get().strip()
                    vmax_str = self.vmax_entry.get().strip()
                    if vmin_str:
                        vmin = float(vmin_str)
                    if vmax_str:
                        vmax = float(vmax_str)
                except ValueError:
                    pass  # Use auto limits if conversion fails
            
            # Get selected colormap
            cmap = self.colormap_var.get()
            
            # Create the plot
            if x_data is not None and y_data is not None:
                # Use pcolormesh for 2D grid data with coordinates
                im = self.output_ax.pcolormesh(x_data, y_data, z_data, shading='auto', 
                                              cmap=cmap, vmin=vmin, vmax=vmax)
                self.output_ax.set_xlabel('X (m)')
                self.output_ax.set_ylabel('Y (m)')
            else:
                # Use imshow if no coordinate data available
                im = self.output_ax.imshow(z_data, cmap=cmap, origin='lower', 
                                          aspect='auto', vmin=vmin, vmax=vmax)
                self.output_ax.set_xlabel('Grid X Index')
                self.output_ax.set_ylabel('Grid Y Index')
            
            # Set title with time step
            title = self.get_variable_title(var_name)
            self.output_ax.set_title(f'{title} (Time step: {time_idx})')
            
            # Handle colorbar properly to avoid shrinking
            cbar_label = self.get_variable_label(var_name)
            self.output_colorbar = self._handle_colorbar(self.output_ax, im, cbar_label, 'output_colorbar')

            # Overlay vegetation if enabled and available
            if self.overlay_veg_var.get() and self.nc_data_cache['veg'] is not None:
                veg_slice = self.nc_data_cache['veg']
                if veg_slice.ndim == 3:
                    veg_data = veg_slice[time_idx, :, :]
                else:
                    veg_data = veg_slice[:, :]

                # Choose plotting method consistent with base plot
                if x_data is not None and y_data is not None:
                    self.output_ax.pcolormesh(x_data, y_data, veg_data, shading='auto', 
                                              cmap='Greens', vmin=0, vmax=1, alpha=0.4)
                else:
                    self.output_ax.imshow(veg_data, cmap='Greens', origin='lower', 
                                          aspect='auto', vmin=0, vmax=1, alpha=0.4)
            
            # Redraw the canvas
            self.output_canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to update 2D plot: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)  # Print to console for debugging

    def render_zb_rhoveg_shaded(self, time_idx):
        """
        Render zb+rhoveg combined visualization with hillshading and vegetation blending.
        Inspired by Anim2D_ShadeVeg.py
        """
        try:
            # Get zb and rhoveg data - check if they exist
            if 'zb' not in self.nc_data_cache['vars']:
                raise ValueError("Variable 'zb' not found in NetCDF cache")
            if 'rhoveg' not in self.nc_data_cache['vars']:
                raise ValueError("Variable 'rhoveg' not found in NetCDF cache")
            
            zb_data = self.nc_data_cache['vars']['zb']
            veg_data = self.nc_data_cache['vars']['rhoveg']
            
            # Extract time slice
            if zb_data.ndim == 4:
                zb = zb_data[time_idx, :, :, :].mean(axis=2)
            else:
                zb = zb_data[time_idx, :, :]
            
            if veg_data.ndim == 4:
                veg = veg_data[time_idx, :, :, :].mean(axis=2)
            else:
                veg = veg_data[time_idx, :, :]
            
            # Ensure zb and veg have the same shape
            if zb.shape != veg.shape:
                raise ValueError(f"Shape mismatch: zb={zb.shape}, veg={veg.shape}")
            
            # Get coordinates
            x_data = self.nc_data_cache['x']
            y_data = self.nc_data_cache['y']
            
            # Convert x, y to 1D arrays if needed
            if x_data is not None and y_data is not None:
                if x_data.ndim == 2:
                    x1d = x_data[0, :].astype(float)
                    y1d = y_data[:, 0].astype(float)
                else:
                    x1d = np.asarray(x_data, dtype=float).ravel()
                    y1d = np.asarray(y_data, dtype=float).ravel()
            else:
                # Use indices if no coordinate data
                x1d = np.arange(zb.shape[1], dtype=float)
                y1d = np.arange(zb.shape[0], dtype=float)
            
            # Normalize vegetation to [0,1]
            veg_max = np.nanmax(veg)
            if veg_max is not None and veg_max > 0:
                veg_norm = np.clip(veg / veg_max, 0.0, 1.0)
            else:
                veg_norm = np.clip(veg, 0.0, 1.0)
            
            # Replace any NaNs with 0
            veg_norm = np.nan_to_num(veg_norm, nan=0.0)
            
            # Apply hillshade to topography
            shaded = apply_hillshade(zb, x1d, y1d)
            
            # Create base color by blending sand and vegetation
            # rgb shape: (ny, nx, 3)
            rgb = COLOR_SAND[None, None, :] * (1.0 - veg_norm[..., None]) + COLOR_DARKGREEN[None, None, :] * veg_norm[..., None]
            
            # Apply ocean mask: zb < -0.5 and x < 200
            if x_data is not None:
                X2d, _ = np.meshgrid(x1d, y1d)
                ocean_mask = (zb < -0.5) & (X2d < 200)
                rgb[ocean_mask] = COLOR_OCEAN
            
            # Apply hillshade to modulate colors
            rgb *= shaded[..., None]
            
            # Clip to valid range
            rgb = np.clip(rgb, 0.0, 1.0)
            
            # Plot the RGB image
            if x_data is not None and y_data is not None:
                extent = [x1d.min(), x1d.max(), y1d.min(), y1d.max()]
                self.output_ax.imshow(rgb, origin='lower', extent=extent, interpolation='nearest', aspect='auto')
                self.output_ax.set_xlabel('X (m)')
                self.output_ax.set_ylabel('Y (m)')
            else:
                self.output_ax.imshow(rgb, origin='lower', interpolation='nearest', aspect='auto')
                self.output_ax.set_xlabel('Grid X Index')
                self.output_ax.set_ylabel('Grid Y Index')
            
            # Set title
            title = self.get_variable_title('zb+rhoveg')
            self.output_ax.set_title(f'{title} (Time step: {time_idx})')
            
            # Remove colorbar for RGB visualization
            self._remove_colorbar('output_colorbar')
            
            # Redraw the canvas
            self.output_canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to render zb+rhoveg: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)
            messagebox.showerror("Error", f"Failed to render zb+rhoveg visualization:\n{str(e)}")

    def render_ustar_quiver(self, time_idx):
        """
        Render quiver plot of shear velocity vectors (ustars, ustarn) overlaid on ustar magnitude.
        Background: color plot of ustar magnitude
        Arrows: black vectors showing direction and magnitude
        """
        try:
            # Get ustar component data - check if they exist
            if 'ustars' not in self.nc_data_cache['vars']:
                raise ValueError("Variable 'ustars' not found in NetCDF cache")
            if 'ustarn' not in self.nc_data_cache['vars']:
                raise ValueError("Variable 'ustarn' not found in NetCDF cache")
            
            ustars_data = self.nc_data_cache['vars']['ustars']
            ustarn_data = self.nc_data_cache['vars']['ustarn']
            
            # Extract time slice
            if ustars_data.ndim == 4:
                ustars = ustars_data[time_idx, :, :, :].mean(axis=2)
            else:
                ustars = ustars_data[time_idx, :, :]
            
            if ustarn_data.ndim == 4:
                ustarn = ustarn_data[time_idx, :, :, :].mean(axis=2)
            else:
                ustarn = ustarn_data[time_idx, :, :]
            
            # Calculate ustar magnitude from components
            ustar = np.sqrt(ustars**2 + ustarn**2)
            
            # Get coordinates
            x_data = self.nc_data_cache['x']
            y_data = self.nc_data_cache['y']
            
            # Get colorbar limits
            vmin = None
            vmax = None
            if not self.auto_limits_var.get():
                try:
                    vmin_str = self.vmin_entry.get().strip()
                    vmax_str = self.vmax_entry.get().strip()
                    if vmin_str:
                        vmin = float(vmin_str)
                    if vmax_str:
                        vmax = float(vmax_str)
                except ValueError:
                    pass  # Use auto limits if conversion fails
            
            # Get selected colormap
            cmap = self.colormap_var.get()
            
            # Plot the background ustar magnitude
            if x_data is not None and y_data is not None:
                # Use pcolormesh for 2D grid data with coordinates
                im = self.output_ax.pcolormesh(x_data, y_data, ustar, shading='auto', 
                                              cmap=cmap, vmin=vmin, vmax=vmax)
                self.output_ax.set_xlabel('X (m)')
                self.output_ax.set_ylabel('Y (m)')
            else:
                # Use imshow if no coordinate data available
                im = self.output_ax.imshow(ustar, cmap=cmap, origin='lower', 
                                          aspect='auto', vmin=vmin, vmax=vmax)
                self.output_ax.set_xlabel('Grid X Index')
                self.output_ax.set_ylabel('Grid Y Index')
            
            # Handle colorbar
            self.output_colorbar = self._handle_colorbar(
                self.output_ax, im, 'Shear Velocity (m/s)', 'output_colorbar'
            )
            
            # Create coordinate arrays for quiver
            if x_data is not None and y_data is not None:
                if x_data.ndim == 2:
                    X = x_data
                    Y = y_data
                else:
                    X, Y = np.meshgrid(x_data, y_data)
            else:
                # Use indices if no coordinate data
                X, Y = np.meshgrid(np.arange(ustars.shape[1]), np.arange(ustars.shape[0]))
            
            # Filter out invalid vectors (NaN, zero magnitude)
            valid = np.isfinite(ustars) & np.isfinite(ustarn)
            magnitude = np.sqrt(ustars**2 + ustarn**2)
            valid = valid & (magnitude > 1e-10)
            
            # Subsample for better visibility (every nth point)
            subsample = max(1, min(ustars.shape[0], ustars.shape[1]) // 25)
            
            X_sub = X[::subsample, ::subsample]
            Y_sub = Y[::subsample, ::subsample]
            ustars_sub = ustars[::subsample, ::subsample]
            ustarn_sub = ustarn[::subsample, ::subsample]
            valid_sub = valid[::subsample, ::subsample]
            
            # Apply mask
            X_plot = X_sub[valid_sub]
            Y_plot = Y_sub[valid_sub]
            U_plot = ustars_sub[valid_sub]
            V_plot = ustarn_sub[valid_sub]
            
            # Overlay quiver plot with black arrows
            if len(X_plot) > 0:
                q = self.output_ax.quiver(X_plot, Y_plot, U_plot, V_plot,
                                          color='black', scale=None, scale_units='xy',
                                          angles='xy', pivot='mid', width=0.003)
                
                # Calculate reference vector magnitude for quiver key
                magnitude_all = np.sqrt(U_plot**2 + V_plot**2)
                if magnitude_all.max() > 0:
                    ref_magnitude = magnitude_all.max() * 0.5
                    qk = self.output_ax.quiverkey(q, 0.9, 0.95, ref_magnitude,
                                                 f'{ref_magnitude:.3f} m/s',
                                                 labelpos='E', coordinates='figure',
                                                 color='black')
            
            # Set title
            title = self.get_variable_title('ustar quiver')
            self.output_ax.set_title(f'{title} (Time step: {time_idx})')
            
            # Redraw the canvas
            self.output_canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to render ustar quiver: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)
            messagebox.showerror("Error", f"Failed to render ustar quiver visualization:\n{str(e)}")

    def plot_data(self, file_key, title):
        """Plot data from specified file (bed_file, ne_file, or veg_file)"""
        try:
            # Clear the previous plot
            self.ax.clear()
            
            # Get the file paths from the entries
            xgrid_file = self.entries['xgrid_file'].get()
            ygrid_file = self.entries['ygrid_file'].get()
            data_file = self.entries[file_key].get()
            
            # Check if files are specified
            if not data_file:
                messagebox.showwarning("Warning", f"No {file_key} specified!")
                return
            
            # Resolve file paths using helper
            data_file_path = self._resolve_file_path(data_file)
            if not os.path.exists(data_file_path):
                messagebox.showerror("Error", f"File not found: {data_file_path}")
                return
            
            # Load data
            z_data = np.loadtxt(data_file_path)
            
            # Try to load x and y grid data if available
            x_data = None
            y_data = None
            
            if xgrid_file:
                xgrid_file_path = self._resolve_file_path(xgrid_file)
                if os.path.exists(xgrid_file_path):
                    x_data = np.loadtxt(xgrid_file_path)
            
            if ygrid_file:
                ygrid_file_path = self._resolve_file_path(ygrid_file)
                if os.path.exists(ygrid_file_path):
                    y_data = np.loadtxt(ygrid_file_path)
            
            # Choose colormap based on data type
            if file_key == 'bed_file':
                cmap = 'terrain'
                label = 'Elevation (m)'
            elif file_key == 'ne_file':
                cmap = 'viridis'
                label = 'Ne'
            elif file_key == 'veg_file':
                cmap = 'Greens'
                label = 'Vegetation'
            else:
                cmap = 'viridis'
                label = 'Value'
            
            # Create the plot
            if x_data is not None and y_data is not None:
                # Use pcolormesh for 2D grid data with coordinates
                im = self.ax.pcolormesh(x_data, y_data, z_data, shading='auto', cmap=cmap)
                self.ax.set_xlabel('X (m)')
                self.ax.set_ylabel('Y (m)')
            else:
                # Use imshow if no coordinate data available
                im = self.ax.imshow(z_data, cmap=cmap, origin='lower', aspect='auto')
                self.ax.set_xlabel('Grid X Index')
                self.ax.set_ylabel('Grid Y Index')
            
            self.ax.set_title(title)
            
            # Handle colorbar properly to avoid shrinking
            self.colorbar = self._handle_colorbar(self.ax, im, label, 'colorbar')

            # Enforce equal aspect ratio in domain visualization
            self.ax.set_aspect('equal', adjustable='box')
            
            # Redraw the canvas
            self.canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to plot {file_key}: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)  # Also print to console for debugging

    def plot_combined(self):
        """Plot bed elevation with vegetation overlay"""
        try:
            # Clear the previous plot
            self.ax.clear()
            
            # Get the file paths from the entries
            xgrid_file = self.entries['xgrid_file'].get()
            ygrid_file = self.entries['ygrid_file'].get()
            bed_file = self.entries['bed_file'].get()
            veg_file = self.entries['veg_file'].get()
            
            # Check if files are specified
            if not bed_file:
                messagebox.showwarning("Warning", "No bed_file specified!")
                return
            if not veg_file:
                messagebox.showwarning("Warning", "No veg_file specified!")
                return
            
            # Resolve file paths using helper
            bed_file_path = self._resolve_file_path(bed_file)
            if not os.path.exists(bed_file_path):
                messagebox.showerror("Error", f"Bed file not found: {bed_file_path}")
                return
            
            veg_file_path = self._resolve_file_path(veg_file)
            if not os.path.exists(veg_file_path):
                messagebox.showerror("Error", f"Vegetation file not found: {veg_file_path}")
                return
            
            # Load data
            bed_data = np.loadtxt(bed_file_path)
            veg_data = np.loadtxt(veg_file_path)
            
            # Try to load x and y grid data if available
            x_data = None
            y_data = None
            
            if xgrid_file:
                xgrid_file_path = self._resolve_file_path(xgrid_file)
                if os.path.exists(xgrid_file_path):
                    x_data = np.loadtxt(xgrid_file_path)
            
            if ygrid_file:
                ygrid_file_path = self._resolve_file_path(ygrid_file)
                if os.path.exists(ygrid_file_path):
                    y_data = np.loadtxt(ygrid_file_path)
            
            # Create the bed elevation plot
            if x_data is not None and y_data is not None:
                # Use pcolormesh for 2D grid data with coordinates
                im = self.ax.pcolormesh(x_data, y_data, bed_data, shading='auto', cmap='terrain')
                self.ax.set_xlabel('X (m)')
                self.ax.set_ylabel('Y (m)')
                
                # Overlay vegetation as contours where vegetation exists
                veg_mask = veg_data > 0
                if np.any(veg_mask):
                    # Create contour lines for vegetation
                    contour = self.ax.contour(x_data, y_data, veg_data, levels=[0.5], 
                                             colors='darkgreen', linewidths=2)
                    # Fill vegetation areas with semi-transparent green
                    contourf = self.ax.contourf(x_data, y_data, veg_data, levels=[0.5, veg_data.max()], 
                                               colors=['green'], alpha=0.3)
            else:
                # Use imshow if no coordinate data available
                im = self.ax.imshow(bed_data, cmap='terrain', origin='lower', aspect='auto')
                self.ax.set_xlabel('Grid X Index')
                self.ax.set_ylabel('Grid Y Index')
                
                # Overlay vegetation
                veg_mask = veg_data > 0
                if np.any(veg_mask):
                    # Create a masked array for vegetation overlay
                    veg_overlay = np.ma.masked_where(~veg_mask, veg_data)
                    self.ax.imshow(veg_overlay, cmap='Greens', origin='lower', aspect='auto', alpha=0.5)
            
            self.ax.set_title('Bed Elevation with Vegetation')
            
            # Handle colorbar properly to avoid shrinking
            self.colorbar = self._handle_colorbar(self.ax, im, 'Elevation (m)', 'colorbar')

            # Enforce equal aspect ratio in domain visualization
            self.ax.set_aspect('equal', adjustable='box')
            
            # Redraw the canvas
            self.canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to plot combined view: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)  # Also print to console for debugging

    def plot_nc_bed_level(self):
        """Plot bed level from NetCDF output file"""
        if not HAVE_NETCDF:
            messagebox.showerror("Error", "netCDF4 library is not available!")
            return
            
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

    def update_time_step(self, value):
        """Update the plot based on the time slider value"""
        if self.nc_data_cache is None:
            return
        
        # Get time index from slider
        time_idx = int(float(value))
        
        # Update label
        self.time_label.config(text=f"Time step: {time_idx}")
        
        # Update the 2D plot
        self.update_2d_plot()
    def plot_nc_wind(self):
        """Plot shear velocity (ustar) from NetCDF output file (uses 'ustar' or computes from 'ustars' and 'ustarn')."""
        if not HAVE_NETCDF:
            messagebox.showerror("Error", "netCDF4 library is not available!")
            return
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

    def enable_overlay_vegetation(self):
        """Enable vegetation overlay in the output plot and load vegetation data if needed"""
        if not HAVE_NETCDF:
            messagebox.showerror("Error", "netCDF4 library is not available!")
            return

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
                
                nc_file_path = self._resolve_file_path(nc_file)
                if not os.path.exists(nc_file_path):
                    messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                    return

                # Try common vegetation variable names (use constant)
                with netCDF4.Dataset(nc_file_path, 'r') as nc:
                    available = set(nc.variables.keys())
                    veg_name = next((v for v in VEG_CANDIDATES if v in available), None)
                    if veg_name is None:
                        messagebox.showerror(
                            "Error",
                            "No vegetation variable found in NetCDF file.\n"
                            f"Tried: {', '.join(VEG_CANDIDATES)}\n"
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
        """Save the current entries to the configuration dictionary"""
        global configfile
        # Save the current entries to the configuration dictionary
        for field, entry in self.entries.items():
            self.dic[field] = entry.get()
        # Write the updated configuration to a new file
        # Use original configfile if available, otherwise use current directory
        if configfile and configfile != "No file selected":
            output_path = configfile + '2'
        else:
            output_path = os.path.join(os.getcwd(), 'aeolis_config2.txt')
        aeolis.inout.write_configfile(output_path, self.dic)
        print(f'Saved to: {output_path}')

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
