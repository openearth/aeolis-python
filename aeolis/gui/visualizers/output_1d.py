"""
1D Output Visualizer Module

Handles visualization of 1D transect data from NetCDF output including:
- Cross-shore and along-shore transects
- Time evolution with slider control
- Domain overview with transect indicator
- PNG and MP4 animation export
"""

import os
import numpy as np
import traceback
import netCDF4
from tkinter import messagebox, filedialog, Toplevel
from tkinter import ttk


from aeolis.gui.utils import (
    NC_COORD_VARS, VARIABLE_LABELS, VARIABLE_TITLES,
    resolve_file_path, extract_time_slice
)


class Output1DVisualizer:
    """
    Visualizer for 1D transect data from NetCDF output.
    
    Handles loading, plotting, and exporting 1D transect visualizations
    with support for time evolution and domain overview.
    """
    
    def __init__(self, transect_ax, overview_ax, transect_canvas, transect_fig,
                 time_slider_1d, time_label_1d, transect_slider, transect_label,
                 variable_var_1d, direction_var, nc_file_entry_1d,
                 variable_dropdown_1d, overview_canvas, get_config_dir_func, 
                 get_variable_label_func, get_variable_title_func):
        """Initialize the 1D output visualizer."""
        self.transect_ax = transect_ax
        self.overview_ax = overview_ax
        self.transect_canvas = transect_canvas
        self.transect_fig = transect_fig
        self.overview_canvas = overview_canvas
        self.time_slider_1d = time_slider_1d
        self.time_label_1d = time_label_1d
        self.transect_slider = transect_slider
        self.transect_label = transect_label
        self.variable_var_1d = variable_var_1d
        self.direction_var = direction_var
        self.nc_file_entry_1d = nc_file_entry_1d
        self.variable_dropdown_1d = variable_dropdown_1d
        self.get_config_dir = get_config_dir_func
        self.get_variable_label = get_variable_label_func
        self.get_variable_title = get_variable_title_func
        
        self.nc_data_cache_1d = None
    
    def load_and_plot(self):
        """Load NetCDF file and plot 1D transect data."""
        try:
            nc_file = self.nc_file_entry_1d.get()
            if not nc_file:
                messagebox.showwarning("Warning", "No NetCDF file specified!")
                return
            
            config_dir = self.get_config_dir()
            nc_file_path = resolve_file_path(nc_file, config_dir)
            if not nc_file_path or not os.path.exists(nc_file_path):
                messagebox.showerror("Error", f"NetCDF file not found: {nc_file_path}")
                return
            
            # Open NetCDF file and cache data
            with netCDF4.Dataset(nc_file_path, 'r') as nc:
                available_vars = list(nc.variables.keys())
                
                # Get coordinates
                x_data = nc.variables['x'][:] if 'x' in nc.variables else None
                y_data = nc.variables['y'][:] if 'y' in nc.variables else None
                
                # Load variables
                var_data_dict = {}
                n_times = 1
                
                for var_name in available_vars:
                    if var_name in NC_COORD_VARS:
                        continue
                    
                    var = nc.variables[var_name]
                    if 'time' in var.dimensions:
                        var_data = var[:]
                        if var_data.ndim < 3:
                            continue
                        n_times = max(n_times, var_data.shape[0])
                    else:
                        if var.ndim != 2:
                            continue
                        var_data = np.expand_dims(var[:, :], axis=0)
                    
                    var_data_dict[var_name] = var_data
                
                if not var_data_dict:
                    messagebox.showerror("Error", "No valid variables found in NetCDF file!")
                    return
                
                # Update UI
                candidate_vars = list(var_data_dict.keys())
                self.variable_dropdown_1d['values'] = sorted(candidate_vars)
                if candidate_vars:
                    self.variable_var_1d.set(candidate_vars[0])
                
                # Cache data
                self.nc_data_cache_1d = {
                    'file_path': nc_file_path,
                    'vars': var_data_dict,
                    'x': x_data,
                    'y': y_data,
                    'n_times': n_times
                }
                
                # Get grid dimensions
                first_var = list(var_data_dict.values())[0]
                n_transects = first_var.shape[1] if self.direction_var.get() == 'cross-shore' else first_var.shape[2]
                
                # Setup sliders
                self.time_slider_1d.config(to=n_times - 1)
                self.time_slider_1d.set(0)
                self.time_label_1d.config(text=f"Time step: 0 / {n_times-1}")
                
                self.transect_slider.config(to=n_transects - 1)
                self.transect_slider.set(n_transects // 2)
                self.transect_label.config(text=f"Transect: {n_transects // 2} / {n_transects-1}")
                
                # Plot initial data
                self.update_plot()
                
        except Exception as e:
            error_msg = f"Failed to load NetCDF: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)
    
    def update_transect_position(self, value):
        """Update transect position from slider."""
        if not self.nc_data_cache_1d:
            return
        
        transect_idx = int(float(value))
        first_var = list(self.nc_data_cache_1d['vars'].values())[0]
        n_transects = first_var.shape[1] if self.direction_var.get() == 'cross-shore' else first_var.shape[2]
        self.transect_label.config(text=f"Transect: {transect_idx} / {n_transects-1}")
        self.update_plot()
    
    def update_time_step(self, value):
        """Update time step from slider."""
        if not self.nc_data_cache_1d:
            return
        
        time_idx = int(float(value))
        n_times = self.nc_data_cache_1d['n_times']
        self.time_label_1d.config(text=f"Time step: {time_idx} / {n_times-1}")
        self.update_plot()
    
    def update_plot(self):
        """Update the 1D transect plot with current settings."""
        if not self.nc_data_cache_1d:
            return
        
        try:
            self.transect_ax.clear()
            
            time_idx = int(self.time_slider_1d.get())
            transect_idx = int(self.transect_slider.get())
            var_name = self.variable_var_1d.get()
            direction = self.direction_var.get()
            
            if var_name not in self.nc_data_cache_1d['vars']:
                messagebox.showwarning("Warning", f"Variable '{var_name}' not found!")
                return
            
            # Get data
            var_data = self.nc_data_cache_1d['vars'][var_name]
            z_data = extract_time_slice(var_data, time_idx)
            
            # Extract transect
            if direction == 'cross-shore':
                transect_data = z_data[transect_idx, :]
                x_data = self.nc_data_cache_1d['x'][transect_idx, :] if self.nc_data_cache_1d['x'].ndim == 2 else self.nc_data_cache_1d['x']
                xlabel = 'Cross-shore distance (m)'
            else:  # along-shore
                transect_data = z_data[:, transect_idx]
                x_data = self.nc_data_cache_1d['y'][:, transect_idx] if self.nc_data_cache_1d['y'].ndim == 2 else self.nc_data_cache_1d['y']
                xlabel = 'Along-shore distance (m)'
            
            # Plot transect
            if x_data is not None:
                self.transect_ax.plot(x_data, transect_data, 'b-', linewidth=2)
                self.transect_ax.set_xlabel(xlabel)
            else:
                self.transect_ax.plot(transect_data, 'b-', linewidth=2)
                self.transect_ax.set_xlabel('Grid Index')
            
            ylabel = self.get_variable_label(var_name)
            self.transect_ax.set_ylabel(ylabel)
            
            title = self.get_variable_title(var_name)
            self.transect_ax.set_title(f'{title} - {direction.capitalize()} (Time: {time_idx}, Transect: {transect_idx})')
            self.transect_ax.grid(True, alpha=0.3)
            
            # Update overview
            self.update_overview(transect_idx)
            
            self.transect_canvas.draw()
            
        except Exception as e:
            error_msg = f"Failed to update 1D plot: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)
    
    def update_overview(self, transect_idx):
        """Update the domain overview showing transect position."""
        if not self.nc_data_cache_1d:
            return
        
        try:
            self.overview_ax.clear()
            
            time_idx = int(self.time_slider_1d.get())
            var_name = self.variable_var_1d.get()
            direction = self.direction_var.get()
            
            if var_name not in self.nc_data_cache_1d['vars']:
                return
            
            # Get data for overview
            var_data = self.nc_data_cache_1d['vars'][var_name]
            z_data = extract_time_slice(var_data, time_idx)
            
            x_data = self.nc_data_cache_1d['x']
            y_data = self.nc_data_cache_1d['y']
            
            # Plot domain overview
            if x_data is not None and y_data is not None:
                im = self.overview_ax.pcolormesh(x_data, y_data, z_data, shading='auto', cmap='terrain')
                
                # Draw transect line
                if direction == 'cross-shore':
                    if x_data.ndim == 2:
                        x_line = x_data[transect_idx, :]
                        y_line = y_data[transect_idx, :]
                    else:
                        x_line = x_data
                        y_line = np.full_like(x_data, y_data[transect_idx] if y_data.ndim == 1 else y_data[transect_idx, 0])
                else:  # along-shore
                    if y_data.ndim == 2:
                        x_line = x_data[:, transect_idx]
                        y_line = y_data[:, transect_idx]
                    else:
                        y_line = y_data
                        x_line = np.full_like(y_data, x_data[transect_idx] if x_data.ndim == 1 else x_data[0, transect_idx])
                
                self.overview_ax.plot(x_line, y_line, 'r-', linewidth=2, label='Transect')
                self.overview_ax.set_xlabel('X (m)')
                self.overview_ax.set_ylabel('Y (m)')
            else:
                im = self.overview_ax.imshow(z_data, cmap='terrain', origin='lower', aspect='auto')
                
                # Draw transect line
                if direction == 'cross-shore':
                    self.overview_ax.axhline(y=transect_idx, color='r', linewidth=2, label='Transect')
                else:
                    self.overview_ax.axvline(x=transect_idx, color='r', linewidth=2, label='Transect')
                
                self.overview_ax.set_xlabel('Grid X')
                self.overview_ax.set_ylabel('Grid Y')
            
            self.overview_ax.set_title('Domain Overview')
            self.overview_ax.legend()
            
            # Redraw the overview canvas
            self.overview_canvas.draw()
            
        except Exception as e:
            error_msg = f"Failed to update overview: {str(e)}"
            print(error_msg)
    
    def export_png(self, default_filename="output_1d.png"):
        """Export current 1D plot as PNG."""
        if not self.transect_fig:
            messagebox.showwarning("Warning", "No plot to export.")
            return None
        
        file_path = filedialog.asksaveasfilename(
            initialdir=self.get_config_dir(),
            title="Save plot as PNG",
            defaultextension=".png",
            initialfile=default_filename,
            filetypes=(("PNG files", "*.png"), ("All files", "*.*"))
        )
        
        if file_path:
            try:
                self.transect_fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Plot exported to:\n{file_path}")
                return file_path
            except Exception as e:
                error_msg = f"Failed to export: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
        return None
    
    def export_animation_mp4(self, default_filename="output_1d_animation.mp4"):
        """Export 1D transect animation as MP4."""
        if not self.nc_data_cache_1d or self.nc_data_cache_1d['n_times'] <= 1:
            messagebox.showwarning("Warning", "Need multiple time steps for animation.")
            return None
        
        file_path = filedialog.asksaveasfilename(
            initialdir=self.get_config_dir(),
            title="Save animation as MP4",
            defaultextension=".mp4",
            initialfile=default_filename,
            filetypes=(("MP4 files", "*.mp4"), ("All files", "*.*"))
        )
        
        if file_path:
            try:
                from matplotlib.animation import FuncAnimation, FFMpegWriter
                
                n_times = self.nc_data_cache_1d['n_times']
                progress_window = Toplevel()
                progress_window.title("Exporting Animation")
                progress_window.geometry("300x100")
                progress_label = ttk.Label(progress_window, text="Creating animation...\nThis may take a few minutes.")
                progress_label.pack(pady=20)
                progress_bar = ttk.Progressbar(progress_window, mode='determinate', maximum=n_times)
                progress_bar.pack(pady=10, padx=20, fill='x')
                progress_window.update()
                
                original_time = int(self.time_slider_1d.get())
                
                def update_frame(frame_num):
                    self.time_slider_1d.set(frame_num)
                    self.update_plot()
                    try:
                        if progress_window.winfo_exists():
                            progress_bar['value'] = frame_num + 1
                            progress_window.update()
                    except:
                        pass  # Window may have been closed
                    return []
                
                ani = FuncAnimation(self.transect_fig, update_frame, frames=n_times,
                                   interval=200, blit=False, repeat=False)
                writer = FFMpegWriter(fps=5, bitrate=1800)
                ani.save(file_path, writer=writer)
                
                # Stop the animation by deleting the animation object
                del ani
                
                self.time_slider_1d.set(original_time)
                self.update_plot()
                
                try:
                    if progress_window.winfo_exists():
                        progress_window.destroy()
                except:
                    pass  # Window already destroyed
                
                messagebox.showinfo("Success", f"Animation exported to:\n{file_path}")
                return file_path
                
            except ImportError:
                messagebox.showerror("Error", "Animation export requires ffmpeg.")
            except Exception as e:
                error_msg = f"Failed to export animation: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
            finally:
                try:
                    if 'progress_window' in locals() and progress_window.winfo_exists():
                        progress_window.destroy()
                except:
                    pass  # Window already destroyed
        return None
