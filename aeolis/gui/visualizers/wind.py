"""
Wind Visualizer Module

Handles visualization of wind input data including:
- Wind speed time series
- Wind direction time series
- Wind rose diagrams
- PNG export for wind plots
"""

import os
import numpy as np
import traceback
from tkinter import messagebox, filedialog
import matplotlib.patches as mpatches
from windrose import WindroseAxes
from aeolis.gui.utils import resolve_file_path, determine_time_unit


class WindVisualizer:
    """
    Visualizer for wind input data (time series and wind rose).
    
    Parameters
    ----------
    wind_speed_ax : matplotlib.axes.Axes
        Axes for wind speed time series
    wind_dir_ax : matplotlib.axes.Axes
        Axes for wind direction time series
    wind_ts_canvas : FigureCanvasTkAgg
        Canvas for time series plots
    wind_ts_fig : matplotlib.figure.Figure
        Figure containing time series
    windrose_fig : matplotlib.figure.Figure
        Figure for wind rose
    windrose_canvas : FigureCanvasTkAgg
        Canvas for wind rose
    get_wind_file_func : callable
        Function to get wind file entry widget
    get_entries_func : callable
        Function to get all entry widgets
    get_config_dir_func : callable
        Function to get configuration directory
    get_dic_func : callable
        Function to get configuration dictionary
    """
    
    def __init__(self, wind_speed_ax, wind_dir_ax, wind_ts_canvas, wind_ts_fig,
                 windrose_fig, windrose_canvas, get_wind_file_func, get_entries_func,
                 get_config_dir_func, get_dic_func):
        self.wind_speed_ax = wind_speed_ax
        self.wind_dir_ax = wind_dir_ax
        self.wind_ts_canvas = wind_ts_canvas
        self.wind_ts_fig = wind_ts_fig
        self.windrose_fig = windrose_fig
        self.windrose_canvas = windrose_canvas
        self.get_wind_file = get_wind_file_func
        self.get_entries = get_entries_func
        self.get_config_dir = get_config_dir_func
        self.get_dic = get_dic_func
        self.wind_data_cache = None
    
    def load_and_plot(self):
        """Load wind file and plot time series and wind rose."""
        try:
            # Get the wind file path
            wind_file = self.get_wind_file().get()
            
            if not wind_file:
                messagebox.showwarning("Warning", "No wind file specified!")
                return
            
            # Get the directory of the config file to resolve relative paths
            config_dir = self.get_config_dir()
            
            # Resolve wind file path
            wind_file_path = resolve_file_path(wind_file, config_dir)
            if not wind_file_path or not os.path.exists(wind_file_path):
                messagebox.showerror("Error", f"Wind file not found: {wind_file_path}")
                return
            
            # Check if we already loaded this file (avoid reloading)
            if self.wind_data_cache and self.wind_data_cache.get('file_path') == wind_file_path:
                # Data already loaded, just return (don't reload)
                return
            
            # Load wind data (time, speed, direction)
            wind_data = np.loadtxt(wind_file_path)
            
            # Check data format
            if wind_data.ndim != 2 or wind_data.shape[1] < 3:
                messagebox.showerror("Error", "Wind file must have at least 3 columns: time, speed, direction")
                return
            
            time = wind_data[:, 0]
            speed = wind_data[:, 1]
            direction = wind_data[:, 2]
            
            # Get wind convention from config
            dic = self.get_dic()
            wind_convention = dic.get('wind_convention', 'nautical')
            
            # Cache the wind data along with file path and convention
            self.wind_data_cache = {
                'file_path': wind_file_path,
                'time': time,
                'speed': speed,
                'direction': direction,
                'convention': wind_convention
            }
            
            # Determine appropriate time unit based on simulation time (tstart and tstop)
            tstart = 0
            tstop = 0
            use_sim_limits = False
            
            try:
                entries = self.get_entries()
                tstart_entry = entries.get('tstart')
                tstop_entry = entries.get('tstop')
                
                if tstart_entry and tstop_entry:
                    tstart = float(tstart_entry.get() or 0)
                    tstop = float(tstop_entry.get() or 0)
                    if tstop > tstart:
                        sim_duration = tstop - tstart  # in seconds
                        use_sim_limits = True
                    else:
                        sim_duration = time[-1] - time[0] if len(time) > 0 else 0
                else:
                    sim_duration = time[-1] - time[0] if len(time) > 0 else 0
            except (ValueError, AttributeError, TypeError):
                sim_duration = time[-1] - time[0] if len(time) > 0 else 0
            
            # Choose appropriate time unit and convert using utility function
            time_unit, time_divisor = determine_time_unit(sim_duration)
            time_converted = time / time_divisor
            
            # Plot wind speed time series
            self.wind_speed_ax.clear()
            self.wind_speed_ax.plot(time_converted, speed, 'b-', linewidth=1.5, zorder=2, label='Wind Speed')
            self.wind_speed_ax.set_xlabel(f'Time ({time_unit})')
            self.wind_speed_ax.set_ylabel('Wind Speed (m/s)')
            self.wind_speed_ax.set_title('Wind Speed Time Series')
            self.wind_speed_ax.grid(True, alpha=0.3, zorder=1)
            
            # Calculate axis limits with 10% padding and add shading
            if use_sim_limits:
                tstart_converted = tstart / time_divisor
                tstop_converted = tstop / time_divisor
                axis_range = tstop_converted - tstart_converted
                padding = 0.1 * axis_range
                xlim_min = tstart_converted - padding
                xlim_max = tstop_converted + padding
                
                self.wind_speed_ax.set_xlim([xlim_min, xlim_max])
                self.wind_speed_ax.axvspan(xlim_min, tstart_converted, alpha=0.15, color='gray', zorder=3)
                self.wind_speed_ax.axvspan(tstop_converted, xlim_max, alpha=0.15, color='gray', zorder=3)
                
                shaded_patch = mpatches.Patch(color='gray', alpha=0.15, label='Outside simulation time')
                self.wind_speed_ax.legend(handles=[shaded_patch], loc='upper right', fontsize=8)
            
            # Plot wind direction time series
            self.wind_dir_ax.clear()
            self.wind_dir_ax.plot(time_converted, direction, 'r-', linewidth=1.5, zorder=2, label='Wind Direction')
            self.wind_dir_ax.set_xlabel(f'Time ({time_unit})')
            self.wind_dir_ax.set_ylabel('Wind Direction (degrees)')
            self.wind_dir_ax.set_title(f'Wind Direction Time Series ({wind_convention} convention)')
            self.wind_dir_ax.set_ylim([0, 360])
            self.wind_dir_ax.grid(True, alpha=0.3, zorder=1)
            
            if use_sim_limits:
                self.wind_dir_ax.set_xlim([xlim_min, xlim_max])
                self.wind_dir_ax.axvspan(xlim_min, tstart_converted, alpha=0.15, color='gray', zorder=3)
                self.wind_dir_ax.axvspan(tstop_converted, xlim_max, alpha=0.15, color='gray', zorder=3)
                
                shaded_patch = mpatches.Patch(color='gray', alpha=0.15, label='Outside simulation time')
                self.wind_dir_ax.legend(handles=[shaded_patch], loc='upper right', fontsize=8)
            
            # Redraw time series canvas
            self.wind_ts_canvas.draw()
            
            # Plot wind rose
            self.plot_windrose(speed, direction, wind_convention)
            
        except Exception as e:
            error_msg = f"Failed to load and plot wind data: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)
    
    def force_reload(self):
        """Force reload of wind data by clearing cache."""
        self.wind_data_cache = None
        self.load_and_plot()
    
    def plot_windrose(self, speed, direction, convention='nautical'):
        """
        Plot wind rose diagram.
        
        Parameters
        ----------
        speed : array
            Wind speed values
        direction : array
            Wind direction values in degrees
        convention : str
            'nautical' or 'cartesian'
        """
        try:
            # Clear the windrose figure
            self.windrose_fig.clear()
            
            # Convert direction based on convention to meteorological standard
            if convention == 'cartesian':
                direction_met = (270 - direction) % 360
            else:
                direction_met = direction
            
            # Create windrose axes
            ax = WindroseAxes.from_ax(fig=self.windrose_fig)
            ax.bar(direction_met, speed, normed=True, opening=0.8, edgecolor='white')
            ax.set_legend(title='Wind Speed (m/s)')
            ax.set_title(f'Wind Rose ({convention} convention)', fontsize=14, fontweight='bold')
            
            # Redraw windrose canvas
            self.windrose_canvas.draw()
            
        except Exception as e:
            error_msg = f"Failed to plot wind rose: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)
            # Create a simple text message instead
            self.windrose_fig.clear()
            ax = self.windrose_fig.add_subplot(111)
            ax.text(0.5, 0.5, 'Wind rose plot failed.\nSee console for details.', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            self.windrose_canvas.draw()
    
    def export_timeseries_png(self, default_filename="wind_timeseries.png"):
        """
        Export the wind time series plot as PNG.
        
        Parameters
        ----------
        default_filename : str
            Default filename for the export dialog
            
        Returns
        -------
        str or None
            Path to saved file, or None if cancelled/failed
        """
        if self.wind_ts_fig is None:
            messagebox.showwarning("Warning", "No wind plot to export. Please load wind data first.")
            return None
        
        file_path = filedialog.asksaveasfilename(
            initialdir=self.get_config_dir(),
            title="Save wind time series as PNG",
            defaultextension=".png",
            initialfile=default_filename,
            filetypes=(("PNG files", "*.png"), ("All files", "*.*"))
        )
        
        if file_path:
            try:
                self.wind_ts_fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Wind time series exported to:\n{file_path}")
                return file_path
            except Exception as e:
                error_msg = f"Failed to export plot: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
        
        return None
    
    def export_windrose_png(self, default_filename="wind_rose.png"):
        """
        Export the wind rose plot as PNG.
        
        Parameters
        ----------
        default_filename : str
            Default filename for the export dialog
            
        Returns
        -------
        str or None
            Path to saved file, or None if cancelled/failed
        """
        if self.windrose_fig is None:
            messagebox.showwarning("Warning", "No wind rose plot to export. Please load wind data first.")
            return None
        
        file_path = filedialog.asksaveasfilename(
            initialdir=self.get_config_dir(),
            title="Save wind rose as PNG",
            defaultextension=".png",
            initialfile=default_filename,
            filetypes=(("PNG files", "*.png"), ("All files", "*.*"))
        )
        
        if file_path:
            try:
                self.windrose_fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Wind rose exported to:\n{file_path}")
                return file_path
            except Exception as e:
                error_msg = f"Failed to export plot: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
        
        return None
