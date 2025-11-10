"""
Domain Visualizer Module

Handles visualization of domain setup including:
- Bed elevation
- Vegetation distribution
- Ne (erodibility) parameter
- Combined bed + vegetation views
"""

import os
import numpy as np
import traceback
from tkinter import messagebox
from aeolis.gui.utils import resolve_file_path


class DomainVisualizer:
    """
    Visualizer for domain setup data (bed elevation, vegetation, etc.).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axes to plot on
    canvas : FigureCanvasTkAgg
        The canvas to draw on
    fig : matplotlib.figure.Figure
        The figure containing the axes
    get_entries_func : callable
        Function to get entry widgets dictionary
    get_config_dir_func : callable
        Function to get configuration directory
    """
    
    def __init__(self, ax, canvas, fig, get_entries_func, get_config_dir_func):
        self.ax = ax
        self.canvas = canvas
        self.fig = fig
        self.get_entries = get_entries_func
        self.get_config_dir = get_config_dir_func
        self.colorbar = None
    
    def _load_grid_data(self, xgrid_file, ygrid_file, config_dir):
        """
        Load x and y grid data if available.
        
        Parameters
        ----------
        xgrid_file : str
            Path to x-grid file (may be relative or absolute)
        ygrid_file : str
            Path to y-grid file (may be relative or absolute)
        config_dir : str
            Base directory for resolving relative paths
            
        Returns
        -------
        tuple
            (x_data, y_data) numpy arrays or (None, None) if not available
        """
        x_data = None
        y_data = None
        
        if xgrid_file:
            xgrid_file_path = resolve_file_path(xgrid_file, config_dir)
            if xgrid_file_path and os.path.exists(xgrid_file_path):
                x_data = np.loadtxt(xgrid_file_path)
        
        if ygrid_file:
            ygrid_file_path = resolve_file_path(ygrid_file, config_dir)
            if ygrid_file_path and os.path.exists(ygrid_file_path):
                y_data = np.loadtxt(ygrid_file_path)
        
        return x_data, y_data
    
    def _get_colormap_and_label(self, file_key):
        """
        Get appropriate colormap and label for a given file type.
        
        Parameters
        ----------
        file_key : str
            File type key ('bed_file', 'ne_file', 'veg_file', etc.)
            
        Returns
        -------
        tuple
            (colormap_name, label_text)
        """
        colormap_config = {
            'bed_file': ('terrain', 'Elevation (m)'),
            'ne_file': ('viridis', 'Ne'),
            'veg_file': ('Greens', 'Vegetation'),
        }
        return colormap_config.get(file_key, ('viridis', 'Value'))
    
    def _update_or_create_colorbar(self, im, label):
        """
        Update existing colorbar or create a new one.
        
        Parameters
        ----------
        im : mappable
            The image/mesh object returned by pcolormesh or imshow
        label : str
            Colorbar label
            
        Returns
        -------
        Colorbar
            The updated or newly created colorbar
        """
        if self.colorbar is not None:
            try:
                # Update existing colorbar
                self.colorbar.update_normal(im)
                self.colorbar.set_label(label)
                return self.colorbar
            except Exception:
                # If update fails, create new one
                pass
        
        # Create new colorbar
        self.colorbar = self.fig.colorbar(im, ax=self.ax, label=label)
        return self.colorbar
    
    def plot_data(self, file_key, title):
        """
        Plot data from specified file (bed_file, ne_file, or veg_file).
        
        Parameters
        ----------
        file_key : str
            Key for the file entry (e.g., 'bed_file', 'ne_file', 'veg_file')
        title : str
            Plot title
        """
        try:
            # Clear the previous plot
            self.ax.clear()
            
            # Get the file paths from the entries
            entries = self.get_entries()
            xgrid_file = entries['xgrid_file'].get()
            ygrid_file = entries['ygrid_file'].get()
            data_file = entries[file_key].get()
            
            # Check if files are specified
            if not data_file:
                messagebox.showwarning("Warning", f"No {file_key} specified!")
                return
            
            # Get the directory of the config file to resolve relative paths
            config_dir = self.get_config_dir()
            
            # Load the data file
            data_file_path = resolve_file_path(data_file, config_dir)
            if not data_file_path or not os.path.exists(data_file_path):
                messagebox.showerror("Error", f"File not found: {data_file_path}")
                return
            
            # Load data
            z_data = np.loadtxt(data_file_path)
            
            # Try to load x and y grid data if available
            x_data, y_data = self._load_grid_data(xgrid_file, ygrid_file, config_dir)
            
            # Choose colormap based on data type
            cmap, label = self._get_colormap_and_label(file_key)
            
            # Use pcolormesh for 2D grid data with coordinates
            im = self.ax.pcolormesh(x_data, y_data, z_data, shading='auto', cmap=cmap)
            self.ax.set_xlabel('X (m)')
            self.ax.set_ylabel('Y (m)')
            
            self.ax.set_title(title)
            
            # Handle colorbar properly to avoid shrinking
            self.colorbar = self._update_or_create_colorbar(im, label)

            # Enforce equal aspect ratio in domain visualization
            self.ax.set_aspect('equal', adjustable='box')
            
            # Redraw the canvas
            self.canvas.draw()
            
        except Exception as e:
            error_msg = f"Failed to plot {file_key}: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)
    
    def plot_combined(self):
        """Plot bed elevation with vegetation overlay."""
        try:
            # Clear the previous plot
            self.ax.clear()
            
            # Get the file paths from the entries
            entries = self.get_entries()
            xgrid_file = entries['xgrid_file'].get()
            ygrid_file = entries['ygrid_file'].get()
            bed_file = entries['bed_file'].get()
            veg_file = entries['veg_file'].get()
            
            # Check if files are specified
            if not bed_file:
                messagebox.showwarning("Warning", "No bed_file specified!")
                return
            if not veg_file:
                messagebox.showwarning("Warning", "No veg_file specified!")
                return
            
            # Get the directory of the config file to resolve relative paths
            config_dir = self.get_config_dir()
            
            # Load the bed file
            bed_file_path = resolve_file_path(bed_file, config_dir)
            if not bed_file_path or not os.path.exists(bed_file_path):
                messagebox.showerror("Error", f"Bed file not found: {bed_file_path}")
                return
            
            # Load the vegetation file
            veg_file_path = resolve_file_path(veg_file, config_dir)
            if not veg_file_path or not os.path.exists(veg_file_path):
                messagebox.showerror("Error", f"Vegetation file not found: {veg_file_path}")
                return
            
            # Load data
            bed_data = np.loadtxt(bed_file_path)
            veg_data = np.loadtxt(veg_file_path)
            
            # Try to load x and y grid data if available
            x_data, y_data = self._load_grid_data(xgrid_file, ygrid_file, config_dir)
            
            # Use pcolormesh for 2D grid data with coordinates
            im = self.ax.pcolormesh(x_data, y_data, bed_data, shading='auto', cmap='terrain')
            self.ax.set_xlabel('X (m)')
            self.ax.set_ylabel('Y (m)')
            
            # Overlay vegetation as contours where vegetation exists
            veg_mask = veg_data > 0
            if np.any(veg_mask):
                # Create contour lines for vegetation
                self.ax.contour(x_data, y_data, veg_data, levels=[0.5], 
                                colors='darkgreen', linewidths=2)
                # Fill vegetation areas with semi-transparent green
                self.ax.contourf(x_data, y_data, veg_data, levels=[0.5, veg_data.max()], 
                                 colors=['green'], alpha=0.3)
            
            self.ax.set_title('Bed Elevation with Vegetation')
            
            # Handle colorbar properly to avoid shrinking
            self.colorbar = self._update_or_create_colorbar(im, 'Elevation (m)')

            # Enforce equal aspect ratio in domain visualization
            self.ax.set_aspect('equal', adjustable='box')
            
            # Redraw the canvas
            self.canvas.draw()
            
        except Exception as e:
            error_msg = f"Failed to plot combined view: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)
    
    def export_png(self, default_filename="domain_plot.png"):
        """
        Export the current domain plot as PNG.
        
        Parameters
        ----------
        default_filename : str
            Default filename for the export dialog
            
        Returns
        -------
        str or None
            Path to saved file, or None if cancelled/failed
        """
        from tkinter import filedialog
        
        # Open file dialog for saving
        file_path = filedialog.asksaveasfilename(
            initialdir=self.get_config_dir(),
            title="Save plot as PNG",
            defaultextension=".png",
            initialfile=default_filename,
            filetypes=(("PNG files", "*.png"), ("All files", "*.*"))
        )
        
        if file_path:
            try:
                # Ensure canvas is drawn before saving
                self.canvas.draw()
                # Use tight layout to ensure everything fits
                self.fig.tight_layout()
                # Save the figure
                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Plot exported to:\n{file_path}")
                return file_path
            except Exception as e:
                error_msg = f"Failed to export plot: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
        
        return None
