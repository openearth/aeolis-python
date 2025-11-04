import aeolis
from tkinter import *
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

try:
    import netCDF4
    HAVE_NETCDF = True
except ImportError:
    HAVE_NETCDF = False

# Default configuration file path
configfile = r'C:\Users\svries\Documents\GitHub\OE_aeolis-python\aeolis\examples\2D\Barchan_dune\aeolis.txt'

# Function to prompt the user to select a configuration file
def prompt_file():
    file_path = filedialog.askopenfilename(
        initialdir=os.path.dirname(configfile),
        title="Select config file",
        filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
    )
    return file_path if file_path else configfile

# Prompt the user to select a configuration file or use the default
configfile = prompt_file()
# Read the configuration file into a dictionary
dic = aeolis.inout.read_configfile(configfile)

class AeolisGUI:
    def __init__(self, root, dic):
        self.root = root
        self.dic = dic
        self.root.title('Aeolis')
        
        # Initialize attributes
        self.nc_data_cache = None
        self.overlay_veg_enabled = False
        
        self.create_widgets()

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

    def create_label_entry(self, tab, text, value, row):
        # Create a label and entry widget for a given tab
        label = ttk.Label(tab, text=text)
        label.grid(row=row, column=0, sticky=W)
        entry = ttk.Entry(tab)
        entry.insert(0, str(value))
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
            entry.insert(0, str(self.dic.get(field, '')))
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
        # Get initial directory from config file location
        initial_dir = os.path.dirname(configfile)
        
        # Get current value to determine initial directory
        current_value = entry_widget.get()
        if current_value:
            if os.path.isabs(current_value):
                initial_dir = os.path.dirname(current_value)
            else:
                full_path = os.path.join(initial_dir, current_value)
                if os.path.exists(full_path):
                    initial_dir = os.path.dirname(full_path)
        
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
            config_dir = os.path.dirname(configfile)
            try:
                rel_path = os.path.relpath(file_path, config_dir)
                # Use relative path if it doesn't go up too many levels
                if not rel_path.startswith('..\\..\\'):
                    file_path = rel_path
            except ValueError:
                # Different drives on Windows, keep absolute path
                pass
            
            entry_widget.delete(0, END)
            entry_widget.insert(0, file_path)

    def browse_nc_file(self):
        """Open file dialog to select a NetCDF file"""
        # Get initial directory from config file location
        initial_dir = os.path.dirname(configfile)
        
        # Get current value to determine initial directory
        current_value = self.nc_file_entry.get()
        if current_value:
            if os.path.isabs(current_value):
                initial_dir = os.path.dirname(current_value)
            else:
                full_path = os.path.join(initial_dir, current_value)
                if os.path.exists(full_path):
                    initial_dir = os.path.dirname(full_path)
        
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
            config_dir = os.path.dirname(configfile)
            try:
                rel_path = os.path.relpath(file_path, config_dir)
                # Use relative path if it doesn't go up too many levels
                if not rel_path.startswith('..\\..\\'):
                    file_path = rel_path
            except ValueError:
                # Different drives on Windows, keep absolute path
                pass
            
            self.nc_file_entry.delete(0, END)
            self.nc_file_entry.insert(0, file_path)

    def load_new_config(self):
        """Load a new configuration file and update all fields"""
        global configfile
        
        # Open file dialog
        file_path = filedialog.askopenfilename(
            initialdir=os.path.dirname(configfile),
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
            initialdir=os.path.dirname(configfile),
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
        file_frame = ttk.LabelFrame(tab5, text="Output File", padding=10)
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

        # Colorbar limits
        vmin_label = ttk.Label(file_frame, text="Color min:")
        vmin_label.grid(row=1, column=0, sticky=W, pady=2)
        self.vmin_entry = ttk.Entry(file_frame, width=15, state='disabled')
        self.vmin_entry.grid(row=1, column=1, sticky=W, pady=2, padx=(0, 5))
        
        vmax_label = ttk.Label(file_frame, text="Color max:")
        vmax_label.grid(row=2, column=0, sticky=W, pady=2)
        self.vmax_entry = ttk.Entry(file_frame, width=15, state='disabled')
        self.vmax_entry.grid(row=2, column=1, sticky=W, pady=2, padx=(0, 5))
        
        # Auto limits checkbox
        self.auto_limits_var = BooleanVar(value=True)
        auto_limits_check = ttk.Checkbutton(file_frame, text="Auto limits", 
                                           variable=self.auto_limits_var,
                                           command=self.toggle_color_limits)
        auto_limits_check.grid(row=1, column=2, rowspan=2, sticky=W, pady=2)

        # Colormap selection
        cmap_label = ttk.Label(file_frame, text="Colormap:")
        cmap_label.grid(row=3, column=0, sticky=W, pady=2)
        
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
        colormap_dropdown.grid(row=3, column=1, sticky=W, pady=2, padx=(0, 5))

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

        # Create a frame for buttons
        output_button_frame = ttk.Frame(plot_frame)
        output_button_frame.pack(pady=5)

        # Create plot button
        plot_bed_button = ttk.Button(output_button_frame, text="Plot Bed Level", 
                                     command=self.plot_nc_bed_level)
        plot_bed_button.grid(row=0, column=0, padx=5)
        
        # Create plot shear velocity button
        plot_wind_button = ttk.Button(output_button_frame, text="Plot Shear Velocity", 
                                      command=self.plot_nc_wind)
        plot_wind_button.grid(row=0, column=1, padx=5)

        # Create apply limits button
        apply_button = ttk.Button(output_button_frame, text="Apply Limits", 
                                  command=self.apply_color_limits)
        apply_button.grid(row=0, column=2, padx=5)

        # Overlay vegetation button
        overlay_button = ttk.Button(output_button_frame, text="Overlay Vegetation", 
                                    command=self.enable_overlay_vegetation)
        overlay_button.grid(row=0, column=3, padx=5)

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
        
        self.variable_options_1d = ['zb', 'ustar', 'ustars', 'ustarn', 'zs', 'zsep']
        self.variable_var_1d = StringVar(value='zb')
        variable_dropdown = ttk.Combobox(file_frame_1d, textvariable=self.variable_var_1d, 
                                        values=self.variable_options_1d, state='readonly', width=13)
        variable_dropdown.grid(row=1, column=1, sticky=W, pady=2, padx=(0, 5))
        variable_dropdown.bind('<<ComboboxSelected>>', self.on_variable_changed)

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

        # Create a frame for buttons
        output_button_frame_1d = ttk.Frame(plot_frame_1d)
        output_button_frame_1d.pack(pady=5)

        # Create plot button
        plot_button_1d = ttk.Button(output_button_frame_1d, text="Load & Plot", 
                                     command=self.plot_1d_transect)
        plot_button_1d.grid(row=0, column=0, padx=5)

    def browse_nc_file_1d(self):
        """Open file dialog to select a NetCDF file for 1D plotting"""
        # Get initial directory from config file location
        initial_dir = os.path.dirname(configfile)
        
        # Get current value to determine initial directory
        current_value = self.nc_file_entry_1d.get()
        if current_value:
            if os.path.isabs(current_value):
                initial_dir = os.path.dirname(current_value)
            else:
                full_path = os.path.join(initial_dir, current_value)
                if os.path.exists(full_path):
                    initial_dir = os.path.dirname(full_path)
        
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
            config_dir = os.path.dirname(configfile)
            try:
                rel_path = os.path.relpath(file_path, config_dir)
                # Use relative path if it doesn't go up too many levels
                if not rel_path.startswith('..\\..\\'):
                    file_path = rel_path
            except ValueError:
                # Different drives on Windows, keep absolute path
                pass
            
            self.nc_file_entry_1d.delete(0, END)
            self.nc_file_entry_1d.insert(0, file_path)

    def on_variable_changed(self, event):
        """Update plot when variable selection changes"""
        if hasattr(self, 'nc_data_cache_1d') and self.nc_data_cache_1d is not None:
            self.update_1d_plot()

    def update_transect_direction(self):
        """Update transect label when direction changes"""
        if self.transect_direction_var.get() == 'cross-shore':
            idx = int(self.transect_slider.get())
            self.transect_label.config(text=f"Y-index: {idx}")
        else:
            idx = int(self.transect_slider.get())
            self.transect_label.config(text=f"X-index: {idx}")
        
        # Update plot if data is loaded
        if hasattr(self, 'nc_data_cache_1d') and self.nc_data_cache_1d is not None:
            self.update_1d_plot()

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
            # Get the NC file path
            nc_file = self.nc_file_entry_1d.get()
            
            if not nc_file:
                messagebox.showwarning("Warning", "No NetCDF file specified!")
                return
            
            # Get the directory of the config file to resolve relative paths
            config_dir = os.path.dirname(configfile)
            
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
                # Get available variables
                available_vars = list(nc.variables.keys())
                
                # Try to get x and y coordinates
                x_data = None
                y_data = None
                
                if 'x' in nc.variables:
                    x_data = nc.variables['x'][:]
                if 'y' in nc.variables:
                    y_data = nc.variables['y'][:]
                
                # Get s and n coordinates (grid indices)
                s_data = None
                n_data = None
                if 's' in nc.variables:
                    s_data = nc.variables['s'][:]
                if 'n' in nc.variables:
                    n_data = nc.variables['n'][:]
                
                # Load all available variables from the dropdown
                var_data_dict = {}
                n_times = 1
                
                for var_name in self.variable_options_1d:
                    if var_name in available_vars:
                        var = nc.variables[var_name]
                        
                        # Check if time dimension exists
                        if 'time' in var.dimensions:
                            # Load all time steps
                            var_data = var[:]
                            n_times = max(n_times, var_data.shape[0])
                        else:
                            # Single time step
                            var_data = var[:, :]
                            var_data = np.expand_dims(var_data, axis=0)  # Add time dimension
                        
                        var_data_dict[var_name] = var_data
                
                # Cache data for slider updates
                self.nc_data_cache_1d = {
                    'vars': var_data_dict,
                    'x': x_data,
                    'y': y_data,
                    's': s_data,
                    'n': n_data,
                    'n_times': n_times,
                    'available_vars': available_vars
                }
            
            # Configure the time slider
            if n_times > 1:
                self.time_slider_1d.configure(from_=0, to=n_times-1)
                self.time_slider_1d.set(n_times - 1)  # Start with last time step
            else:
                self.time_slider_1d.configure(from_=0, to=0)
                self.time_slider_1d.set(0)
            
            # Configure transect slider based on data shape
            # Get shape from first available variable
            first_var = next(iter(var_data_dict.values()))
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
            
            # Extract transect based on direction
            if self.transect_direction_var.get() == 'cross-shore':
                # Fix y-index (n), vary along x (s)
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
            
            # Set ylabel based on variable
            ylabel_dict = {
                'zb': 'Bed Elevation (m)',
                'ustar': 'Shear Velocity (m/s)',
                'ustars': 'Shear Velocity S-component (m/s)',
                'ustarn': 'Shear Velocity N-component (m/s)',
                'zs': 'Surface Elevation (m)',
                'zsep': 'Separation Elevation (m)'
            }
            ylabel = ylabel_dict.get(var_name, var_name)
            self.output_1d_ax.set_ylabel(ylabel)
            
            # Set title
            direction = 'Cross-shore' if self.transect_direction_var.get() == 'cross-shore' else 'Along-shore'
            idx_label = 'Y' if self.transect_direction_var.get() == 'cross-shore' else 'X'
            self.output_1d_ax.set_title(f'{direction} Transect: {var_name} ({idx_label}-index={transect_idx}, Time={time_idx})')
            
            # Add grid
            self.output_1d_ax.grid(True, alpha=0.3)
            
            # Redraw the canvas
            self.output_1d_canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to update 1D plot: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)  # Print to console for debugging

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
            
            # Get the directory of the config file to resolve relative paths
            config_dir = os.path.dirname(configfile)
            
            # Load the data file
            if not os.path.isabs(data_file):
                data_file_path = os.path.join(config_dir, data_file)
            else:
                data_file_path = data_file
                
            if not os.path.exists(data_file_path):
                messagebox.showerror("Error", f"File not found: {data_file_path}")
                return
            
            # Load data
            z_data = np.loadtxt(data_file_path)
            
            # Try to load x and y grid data if available
            x_data = None
            y_data = None
            
            if xgrid_file:
                xgrid_file_path = os.path.join(config_dir, xgrid_file) if not os.path.isabs(xgrid_file) else xgrid_file
                if os.path.exists(xgrid_file_path):
                    x_data = np.loadtxt(xgrid_file_path)
            
            if ygrid_file:
                ygrid_file_path = os.path.join(config_dir, ygrid_file) if not os.path.isabs(ygrid_file) else ygrid_file
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
            if self.colorbar is not None:
                # Update existing colorbar
                self.colorbar.update_normal(im)
                self.colorbar.set_label(label)
            else:
                # Create new colorbar only on first run
                self.colorbar = self.fig.colorbar(im, ax=self.ax, label=label)

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
            
            # Get the directory of the config file to resolve relative paths
            config_dir = os.path.dirname(configfile)
            
            # Load the bed file
            if not os.path.isabs(bed_file):
                bed_file_path = os.path.join(config_dir, bed_file)
            else:
                bed_file_path = bed_file
                
            if not os.path.exists(bed_file_path):
                messagebox.showerror("Error", f"Bed file not found: {bed_file_path}")
                return
            
            # Load the vegetation file
            if not os.path.isabs(veg_file):
                veg_file_path = os.path.join(config_dir, veg_file)
            else:
                veg_file_path = veg_file
                
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
                xgrid_file_path = os.path.join(config_dir, xgrid_file) if not os.path.isabs(xgrid_file) else xgrid_file
                if os.path.exists(xgrid_file_path):
                    x_data = np.loadtxt(xgrid_file_path)
            
            if ygrid_file:
                ygrid_file_path = os.path.join(config_dir, ygrid_file) if not os.path.isabs(ygrid_file) else ygrid_file
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
            if self.colorbar is not None:
                # Update existing colorbar
                self.colorbar.update_normal(im)
                self.colorbar.set_label('Elevation (m)')
            else:
                # Create new colorbar only on first run
                self.colorbar = self.fig.colorbar(im, ax=self.ax, label='Elevation (m)')

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
            config_dir = os.path.dirname(configfile)
            
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
        
        try:
            # Get time index from slider
            time_idx = int(float(value))
            
            # Update label
            self.time_label.config(text=f"Time step: {time_idx}")
            
            # Clear the previous plot
            self.output_ax.clear()
            
            # Get data from cache
            # Determine which variable to plot (default to 'zb')
            plot_key = getattr(self, 'output_plot_state', {}).get('key', 'zb')
            z_data = self.nc_data_cache.get(plot_key)
            if z_data is None:
                # Fallback to bed if desired key missing
                plot_key = 'zb'
                z_data = self.nc_data_cache['zb']
            # Select time slice
            z_data = z_data[time_idx, :, :]
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
            title_base = getattr(self, 'output_plot_state', {}).get('title', 'Bed Elevation')
            self.output_ax.set_title(f'{title_base} (Time step: {time_idx})')
            
            # Handle colorbar properly to avoid shrinking
            if self.output_colorbar is not None:
                # Update existing colorbar
                self.output_colorbar.update_normal(im)
                cbar_label = getattr(self, 'output_plot_state', {}).get('label', 'Elevation (m)')
                self.output_colorbar.set_label(cbar_label)
            else:
                # Create new colorbar only on first run
                cbar_label = getattr(self, 'output_plot_state', {}).get('label', 'Elevation (m)')
                self.output_colorbar = self.output_fig.colorbar(im, ax=self.output_ax, label=cbar_label)

            # Optionally overlay vegetation if enabled and available in cache
            if getattr(self, 'overlay_veg_enabled', False) and 'veg' in self.nc_data_cache:
                veg_slice = self.nc_data_cache['veg']
                # veg_slice may be 3D (time,y,x) or 2D (y,x)
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
            
            # Add quiver overlay for shear velocity vectors if plotting ustar and components available
            if plot_key == 'ustar' and 'ustars' in self.nc_data_cache and 'ustarn' in self.nc_data_cache:
                ustars_slice = self.nc_data_cache['ustars'][time_idx, :, :]
                ustarn_slice = self.nc_data_cache['ustarn'][time_idx, :, :]
                
                # Subsample for cleaner quiver plot
                skip = max(1, min(ustars_slice.shape) // 20)  # ~20 arrows per dimension
                
                # Filter out invalid values (zeros, NaNs, infs) to avoid quiver warnings
                ustars_sub = ustars_slice[::skip, ::skip]
                ustarn_sub = ustarn_slice[::skip, ::skip]
                
                # Create mask for valid (non-zero, finite) vectors
                valid_mask = (
                    np.isfinite(ustars_sub) & 
                    np.isfinite(ustarn_sub) & 
                    ((np.abs(ustars_sub) > 1e-10) | (np.abs(ustarn_sub) > 1e-10))
                )
                
                if np.any(valid_mask):
                    if x_data is not None and y_data is not None:
                        # Use actual coordinates
                        x_sub = x_data[::skip, ::skip][valid_mask]
                        y_sub = y_data[::skip, ::skip][valid_mask]
                        us_sub = ustars_sub[valid_mask]
                        un_sub = ustarn_sub[valid_mask]
                        self.output_ax.quiver(
                            x_sub, y_sub, us_sub, un_sub,
                            color='black', alpha=0.6, scale_units='xy', width=0.003
                        )
                    else:
                        # Use indices
                        ny, nx = ustars_slice.shape
                        Y, X = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
                        x_sub = X[::skip, ::skip][valid_mask]
                        y_sub = Y[::skip, ::skip][valid_mask]
                        us_sub = ustars_sub[valid_mask]
                        un_sub = ustarn_sub[valid_mask]
                        self.output_ax.quiver(
                            x_sub, y_sub, us_sub, un_sub,
                            color='black', alpha=0.6, scale_units='xy', width=0.003
                        )
            
            # Redraw the canvas
            self.output_canvas.draw()
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to update time step: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)  # Print to console for debugging

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
            config_dir = os.path.dirname(configfile)
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
                config_dir = os.path.dirname(configfile)
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
    # Start the Tkinter event loop
    root.mainloop()
