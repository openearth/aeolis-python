"""
Model Runner Module

Handles running AeoLiS model simulations from the GUI including:
- Model execution in separate thread
- Real-time logging output capture
- Start/stop controls
- Progress indication
"""

import os
import threading
import logging
import traceback
from tkinter import messagebox, END, NORMAL, DISABLED


class ModelRunner:
    """
    Model runner for executing AeoLiS simulations from GUI.
    
    Handles model execution in a separate thread with real-time logging
    output and user controls for starting/stopping the model.
    """
    
    def __init__(self, start_btn, stop_btn, progress_bar, status_label, 
                 output_text, config_label, root, get_config_func):
        """Initialize the model runner."""
        self.start_btn = start_btn
        self.stop_btn = stop_btn
        self.progress_bar = progress_bar
        self.status_label = status_label
        self.output_text = output_text
        self.config_label = config_label
        self.root = root
        self.get_config = get_config_func
        
        self.model_runner = None
        self.model_thread = None
        self.model_running = False
        
    def start_model(self):
        """Start the AeoLiS model run in a separate thread"""
        configfile = self.get_config()
        
        # Check if config file is selected
        if not configfile or configfile == "No file selected":
            messagebox.showerror("Error", "Please select a configuration file first in the 'Read/Write Inputfile' tab.")
            return
        
        if not os.path.exists(configfile):
            messagebox.showerror("Error", f"Configuration file not found:\n{configfile}")
            return
        
        # Update UI
        self.config_label.config(text=os.path.basename(configfile), foreground="black")
        self.status_label.config(text="Initializing model...", foreground="orange")
        self.start_btn.config(state=DISABLED)
        self.stop_btn.config(state=NORMAL)
        self.progress_bar.start(10)
        
        # Clear output text
        self.output_text.delete(1.0, END)
        self.append_output("="*60 + "\n")
        self.append_output(f"Starting AeoLiS model\n")
        self.append_output(f"Config file: {configfile}\n")
        self.append_output("="*60 + "\n\n")
        
        # Run model in separate thread to prevent GUI freezing
        self.model_running = True
        self.model_thread = threading.Thread(target=self.run_model_thread, 
                                             args=(configfile,), daemon=True)
        self.model_thread.start()
        
    def stop_model(self):
        """Stop the running model"""
        if self.model_running:
            self.model_running = False
            self.status_label.config(text="Stopping model...", foreground="red")
            self.append_output("\n" + "="*60 + "\n")
            self.append_output("STOP requested by user\n")
            self.append_output("="*60 + "\n")
            
    def run_model_thread(self, configfile):
        """Run the model in a separate thread"""
        try:
            # Import here to avoid issues if aeolis.model is not available
            from aeolis.model import AeoLiSRunner
            
            # Create custom logging handler to capture output
            class TextHandler(logging.Handler):
                def __init__(self, text_widget, gui_callback):
                    super().__init__()
                    self.text_widget = text_widget
                    self.gui_callback = gui_callback
                    
                def emit(self, record):
                    msg = self.format(record)
                    # Schedule GUI update from main thread
                    self.gui_callback(msg + "\n")
            
            # Update status
            self.root.after(0, lambda: self.status_label.config(
                text="Running model...", foreground="green"))
            
            # Create model runner
            self.model_runner = AeoLiSRunner(configfile=configfile)
            
            # Set up logging to capture to GUI
            logger = logging.getLogger('aeolis')
            text_handler = TextHandler(self.output_text, self.append_output_threadsafe)
            text_handler.setLevel(logging.INFO)
            text_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s',
                                                       datefmt='%H:%M:%S'))
            logger.addHandler(text_handler)
            
            # Run the model with a callback to check for stop requests
            def check_stop(model):
                if not self.model_running:
                    raise KeyboardInterrupt("Model stopped by user")
            
            try:
                self.model_runner.run(callback=check_stop)
                
                # Model completed successfully
                self.root.after(0, lambda: self.status_label.config(
                    text="Model completed successfully!", foreground="green"))
                self.append_output_threadsafe("\n" + "="*60 + "\n")
                self.append_output_threadsafe("Model run completed successfully!\n")
                self.append_output_threadsafe("="*60 + "\n")
                
            except KeyboardInterrupt:
                self.root.after(0, lambda: self.status_label.config(
                    text="Model stopped by user", foreground="red"))
            except Exception as e:
                error_msg = f"Model error: {str(e)}"
                self.append_output_threadsafe(f"\nERROR: {error_msg}\n")
                self.append_output_threadsafe(traceback.format_exc())
                self.root.after(0, lambda: self.status_label.config(
                    text="Model failed - see output", foreground="red"))
            finally:
                # Clean up
                logger.removeHandler(text_handler)
                
        except Exception as e:
            error_msg = f"Failed to start model: {str(e)}\n{traceback.format_exc()}"
            self.append_output_threadsafe(error_msg)
            self.root.after(0, lambda: self.status_label.config(
                text="Failed to start model", foreground="red"))
        
        finally:
            # Reset UI
            self.model_running = False
            self.root.after(0, self.reset_ui)
            
    def append_output(self, text):
        """Append text to the output widget (must be called from main thread)"""
        self.output_text.insert(END, text)
        self.output_text.see(END)
        self.output_text.update_idletasks()
        
    def append_output_threadsafe(self, text):
        """Thread-safe version of append_output"""
        self.root.after(0, lambda: self.append_output(text))
        
    def reset_ui(self):
        """Reset the UI elements after model run"""
        self.start_btn.config(state=NORMAL)
        self.stop_btn.config(state=DISABLED)
        self.progress_bar.stop()
        
    def update_config_display(self, configfile):
        """Update the config file display label"""
        if configfile and configfile != "No file selected":
            self.config_label.config(text=os.path.basename(configfile), foreground="black")
        else:
            self.config_label.config(text="No file selected", foreground="gray")
