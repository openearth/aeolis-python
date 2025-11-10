"""
2D Output Visualizer Module

Handles visualization of 2D NetCDF output data including:
- Variable selection and plotting
- Time slider control
- Colorbar customization
- Special renderings (hillshade, quiver plots)
- PNG and MP4 export
"""

import os
import numpy as np
import traceback
import netCDF4
from tkinter import messagebox, filedialog, Toplevel
from tkinter import ttk
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from aeolis.gui.utils import (
    NC_COORD_VARS,
    resolve_file_path, extract_time_slice, apply_hillshade
)


class Output2DVisualizer:
    """
    Visualizer for 2D NetCDF output data.
    
    Handles loading, plotting, and exporting 2D output visualizations with
    support for multiple variables, time evolution, and special renderings.
    """
    
    def __init__(self, output_ax, output_canvas, output_fig,
                 output_colorbar_ref, time_slider, time_label,
                 variable_var_2d, colormap_var, auto_limits_var,
                 vmin_entry, vmax_entry, overlay_veg_var,
                 nc_file_entry, variable_dropdown_2d,
                 get_config_dir_func, get_variable_label_func, get_variable_title_func):
        """Initialize the 2D output visualizer."""
        self.output_ax = output_ax
        self.output_canvas = output_canvas
        self.output_fig = output_fig
        self.output_colorbar_ref = output_colorbar_ref
        self.time_slider = time_slider
        self.time_label = time_label
        self.variable_var_2d = variable_var_2d
        self.colormap_var = colormap_var
        self.auto_limits_var = auto_limits_var
        self.vmin_entry = vmin_entry
        self.vmax_entry = vmax_entry
        self.overlay_veg_var = overlay_veg_var
        self.nc_file_entry = nc_file_entry
        self.variable_dropdown_2d = variable_dropdown_2d
        self.get_config_dir = get_config_dir_func
        self.get_variable_label = get_variable_label_func
        self.get_variable_title = get_variable_title_func
        
        self.nc_data_cache = None

    def on_variable_changed(self, event=None):
        """Handle variable selection change."""
        self.update_plot()

    def load_and_plot(self):
        """Load NetCDF file and plot 2D data."""
        try:
            nc_file = self.nc_file_entry.get()
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
                veg_data = None
                
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
                
                # Load vegetation if requested
                if self.overlay_veg_var.get():
                    for veg_name in ['rhoveg', 'vegetated', 'hveg', 'vegfac']:
                        if veg_name in available_vars:
                            veg_var = nc.variables[veg_name]
                            veg_data = veg_var[:] if 'time' in veg_var.dimensions else np.expand_dims(veg_var[:, :], axis=0)
                            break
                
                if not var_data_dict:
                    messagebox.showerror("Error", "No valid variables found in NetCDF file!")
                    return
                
                # Add special options
                candidate_vars = list(var_data_dict.keys())
                if 'zb' in var_data_dict and 'rhoveg' in var_data_dict:
                    candidate_vars.append('zb+rhoveg')
                if 'ustarn' in var_data_dict and 'ustars' in var_data_dict:
                    candidate_vars.append('ustar quiver')
                
                # Update UI
                self.variable_dropdown_2d['values'] = sorted(candidate_vars)
                if candidate_vars:
                    self.variable_var_2d.set(candidate_vars[0])
                
                # Cache data
                self.nc_data_cache = {
                    'file_path': nc_file_path,
                    'vars': var_data_dict,
                    'x': x_data,
                    'y': y_data,
                    'n_times': n_times,
                    'veg': veg_data
                }
                
                # Setup time slider
                self.time_slider.config(to=n_times - 1)
                self.time_slider.set(0)
                self.time_label.config(text=f"Time step: 0 / {n_times-1}")
                
                # Plot initial data
                self.update_plot()
                
        except Exception as e:
            error_msg = f"Failed to load NetCDF: {str(e)}\n\n{traceback.format_exc()}"
            messagebox.showerror("Error", error_msg)
            print(error_msg)

    def update_plot(self):
        """Update the 2D plot with current settings."""
        if not self.nc_data_cache:
            return
        
        try:
            self.output_ax.clear()
            time_idx = int(self.time_slider.get())
            var_name = self.variable_var_2d.get()
            
            # Update time label
            n_times = self.nc_data_cache.get('n_times', 1)
            self.time_label.config(text=f"Time step: {time_idx} / {n_times-1}")
            
            # Special renderings
            if var_name == 'zb+rhoveg':
                self._render_zb_rhoveg_shaded(time_idx)
                return
            if var_name == 'ustar quiver':
                self._render_ustar_quiver(time_idx)
                return
            
            if var_name not in self.nc_data_cache['vars']:
                messagebox.showwarning("Warning", f"Variable '{var_name}' not found!")
                return
            
            # Get data
            var_data = self.nc_data_cache['vars'][var_name]
            z_data = extract_time_slice(var_data, time_idx)
            x_data = self.nc_data_cache['x']
            y_data = self.nc_data_cache['y']
            
            # Get colorbar limits
            vmin, vmax = None, None
            if not self.auto_limits_var.get():
                try:
                    vmin_str = self.vmin_entry.get().strip()
                    vmax_str = self.vmax_entry.get().strip()
                    vmin = float(vmin_str) if vmin_str else None
                    vmax = float(vmax_str) if vmax_str else None
                except ValueError:
                    messagebox.showwarning(
                        "Invalid Input",
                        "Colorbar limits must be valid numbers. Using automatic limits instead."
                    )
            
            cmap = self.colormap_var.get()
            
            # Plot with pcolormesh (x and y always exist in AeoLiS NetCDF files)
            im = self.output_ax.pcolormesh(x_data, y_data, z_data, shading='auto',
                                          cmap=cmap, vmin=vmin, vmax=vmax)
            self.output_ax.set_xlabel('X (m)')
            self.output_ax.set_ylabel('Y (m)')
            
            title = self.get_variable_title(var_name)
            self.output_ax.set_title(f'{title} (Time step: {time_idx})')
            
            # Update colorbar
            self._update_colorbar(im, var_name)
            
            # Overlay vegetation
            if self.overlay_veg_var.get() and self.nc_data_cache['veg'] is not None:
                veg_slice = self.nc_data_cache['veg']
                veg_data = veg_slice[time_idx, :, :] if veg_slice.ndim == 3 else veg_slice[:, :]
                self.output_ax.pcolormesh(x_data, y_data, veg_data, shading='auto',
                                        cmap='Greens', vmin=0, vmax=1, alpha=0.4)
            
            self.output_canvas.draw_idle()
            
        except Exception as e:
            error_msg = f"Failed to update 2D plot: {str(e)}\n\n{traceback.format_exc()}"
            print(error_msg)
    
    def _update_colorbar(self, im, var_name):
        """Update or create colorbar."""
        cbar_label = self.get_variable_label(var_name)
        if self.output_colorbar_ref[0] is not None:
            try:
                self.output_colorbar_ref[0].update_normal(im)
                self.output_colorbar_ref[0].set_label(cbar_label)
            except Exception:
                self.output_colorbar_ref[0] = self.output_fig.colorbar(im, ax=self.output_ax, label=cbar_label)
        else:
            self.output_colorbar_ref[0] = self.output_fig.colorbar(im, ax=self.output_ax, label=cbar_label)

    def export_png(self, default_filename="output_2d.png"):
        """Export current 2D plot as PNG."""
        if not self.output_fig:
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
                self.output_fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Plot exported to:\n{file_path}")
                return file_path
            except Exception as e:
                error_msg = f"Failed to export: {str(e)}\n\n{traceback.format_exc()}"
                messagebox.showerror("Error", error_msg)
                print(error_msg)
        return None
    
    def export_animation_mp4(self, default_filename="output_2d_animation.mp4"):
        """Export 2D plot animation as MP4."""
        if not self.nc_data_cache or self.nc_data_cache['n_times'] <= 1:
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
                
                n_times = self.nc_data_cache['n_times']
                progress_window = Toplevel()
                progress_window.title("Exporting Animation")
                progress_window.geometry("300x100")
                progress_label = ttk.Label(progress_window, text="Creating animation...\nThis may take a few minutes.")
                progress_label.pack(pady=20)
                progress_bar = ttk.Progressbar(progress_window, mode='determinate', maximum=n_times)
                progress_bar.pack(pady=10, padx=20, fill='x')
                progress_window.update()
                
                original_time = int(self.time_slider.get())
                
                def update_frame(frame_num):
                    self.time_slider.set(frame_num)
                    self.update_plot()
                    try:
                        if progress_window.winfo_exists():
                            progress_bar['value'] = frame_num + 1
                            progress_window.update()
                    except:
                        pass  # Window may have been closed
                    return []
                
                ani = FuncAnimation(self.output_fig, update_frame, frames=n_times,
                                   interval=200, blit=False, repeat=False)
                writer = FFMpegWriter(fps=5, bitrate=1800)
                ani.save(file_path, writer=writer)
                
                # Stop the animation by deleting the animation object
                del ani
                
                self.time_slider.set(original_time)
                self.update_plot()
                
                try:
                    if progress_window.winfo_exists():
                        progress_window.destroy()
                except Exception:
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
                except Exception:
                    pass  # Window already destroyed
        return None

    def _render_zb_rhoveg_shaded(self, time_idx):
        """Render combined bed + vegetation with hillshading matching Anim2D_ShadeVeg.py."""
        try:
            zb_data = extract_time_slice(self.nc_data_cache['vars']['zb'], time_idx)
            rhoveg_data = extract_time_slice(self.nc_data_cache['vars']['rhoveg'], time_idx)
            x_data = self.nc_data_cache['x']
            y_data = self.nc_data_cache['y']
            
            # Normalize vegetation to [0,1]
            veg_max = np.nanmax(rhoveg_data)
            veg_norm = rhoveg_data / veg_max if (veg_max is not None and veg_max > 0) else np.clip(rhoveg_data, 0.0, 1.0)
            veg_norm = np.clip(veg_norm, 0.0, 1.0)
            
            # Apply hillshade
            x1d = x_data[0, :] if x_data.ndim == 2 else x_data
            y1d = y_data[:, 0] if y_data.ndim == 2 else y_data
            hillshade = apply_hillshade(zb_data, x1d, y1d, az_deg=155.0, alt_deg=5.0)
            
            # Color definitions
            sand = np.array([1.0, 239.0/255.0, 213.0/255.0])  # light sand
            darkgreen = np.array([34/255, 139/255, 34/255])
            ocean = np.array([70/255, 130/255, 180/255])  # steelblue
            
            # Create RGB array (ny, nx, 3)
            ny, nx = zb_data.shape
            rgb = np.zeros((ny, nx, 3), dtype=float)
            
            # Base color: blend sand and vegetation
            for i in range(3):  # R, G, B channels
                rgb[:, :, i] = sand[i] * (1.0 - veg_norm) + darkgreen[i] * veg_norm
            
            # Apply ocean mask: zb < -0.5 and x < 200
            if x_data is not None:
                X2d = x_data if x_data.ndim == 2 else np.meshgrid(x1d, y1d)[0]
                ocean_mask = (zb_data < -0.5) & (X2d < 200)
                rgb[ocean_mask] = ocean
            
            # Apply shading to all RGB channels
            rgb *= hillshade[:, :, np.newaxis]
            rgb = np.clip(rgb, 0.0, 1.0)
            
            # Plot RGB image
            extent = [x1d.min(), x1d.max(), y1d.min(), y1d.max()]
            self.output_ax.imshow(rgb, origin='lower', extent=extent, 
                                 interpolation='nearest', aspect='auto')
            self.output_ax.set_xlabel('X (m)')
            self.output_ax.set_ylabel('Y (m)')
            
            self.output_ax.set_title(f'Bed + Vegetation (Time step: {time_idx})')
            
            # Get colorbar limits for vegetation
            vmin, vmax = 0, veg_max
            if not self.auto_limits_var.get():
                try:
                    vmin_str = self.vmin_entry.get().strip()
                    vmax_str = self.vmax_entry.get().strip()
                    vmin = float(vmin_str) if vmin_str else 0
                    vmax = float(vmax_str) if vmax_str else veg_max
                except ValueError:
                    pass  # Use default limits if invalid input
            
            # Create a ScalarMappable for the colorbar (showing vegetation density)
            norm = Normalize(vmin=vmin, vmax=vmax)
            sm = ScalarMappable(cmap='Greens', norm=norm)
            sm.set_array(rhoveg_data)
            
            # Add colorbar for vegetation density
            self._update_colorbar(sm, 'rhoveg')
            
            self.output_canvas.draw_idle()
        except Exception as e:
            print(f"Failed to render zb+rhoveg: {e}")
            traceback.print_exc()
    
    def _render_ustar_quiver(self, time_idx):
        """Render quiver plot of shear velocity with magnitude background."""
        try:
            ustarn = extract_time_slice(self.nc_data_cache['vars']['ustarn'], time_idx)
            ustars = extract_time_slice(self.nc_data_cache['vars']['ustars'], time_idx)
            x_data = self.nc_data_cache['x']
            y_data = self.nc_data_cache['y']
            
            # Calculate magnitude for background coloring
            ustar_mag = np.sqrt(ustarn**2 + ustars**2)
            
            # Subsample for quiver
            step = max(1, min(ustarn.shape) // 25)
            
            # Get colormap and limits
            cmap = self.colormap_var.get()
            vmin, vmax = None, None
            if not self.auto_limits_var.get():
                try:
                    vmin_str = self.vmin_entry.get().strip()
                    vmax_str = self.vmax_entry.get().strip()
                    vmin = float(vmin_str) if vmin_str else None
                    vmax = float(vmax_str) if vmax_str else None
                except ValueError:
                    pass  # Use auto limits
            
            # Plot background field (magnitude)
            im = self.output_ax.pcolormesh(x_data, y_data, ustar_mag, 
                                          shading='auto', cmap=cmap, 
                                          vmin=vmin, vmax=vmax, alpha=0.7)
            
            # Calculate appropriate scaling for arrows
            x1d = x_data[0, :] if x_data.ndim == 2 else x_data
            y1d = y_data[:, 0] if y_data.ndim == 2 else y_data
            x_range = x1d.max() - x1d.min()
            y_range = y1d.max() - y1d.min()
            
            # Calculate typical velocity magnitude (handle masked arrays)
            valid_mag = np.asarray(ustar_mag[ustar_mag > 0])
            typical_vel = np.percentile(valid_mag, 75) if valid_mag.size > 0 else 1.0
            arrow_scale = typical_vel * 20  # Scale factor to make arrows visible
            
            # Add quiver plot with black arrows
            Q = self.output_ax.quiver(x_data[::step, ::step], y_data[::step, ::step],
                                     ustars[::step, ::step], ustarn[::step, ::step],
                                     scale=arrow_scale, color='black', width=0.004,
                                     headwidth=3, headlength=4, headaxislength=3.5,
                                     zorder=10)
            
            # Add quiver key (legend for arrow scale) - placed to the right, above colorbar
            self.output_ax.quiverkey(Q, 1.1, 1.05, typical_vel,
                                    f'{typical_vel:.2f} m/s',
                                    labelpos='N', coordinates='axes',
                                    color='black', labelcolor='black',
                                    fontproperties={'size': 9})
            
            self.output_ax.set_xlabel('X (m)')
            self.output_ax.set_ylabel('Y (m)')
            self.output_ax.set_title(f'Shear Velocity (Time step: {time_idx})')
            
            # Update colorbar for magnitude
            self._update_colorbar(im, 'ustar magnitude')
            
            self.output_canvas.draw_idle()
        except Exception as e:
            print(f"Failed to render ustar quiver: {e}")
            traceback.print_exc()
