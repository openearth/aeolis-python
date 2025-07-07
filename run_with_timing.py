#!/usr/bin/env python3
"""
Script to run AeoLiS model with timing analysis
"""

import os
import time
import sys
import logging
from datetime import datetime, timedelta
import traceback
import cProfile
import pstats
import io

# Try to import psutil for memory monitoring (optional)
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False
    print("Note: psutil not available, memory monitoring will be skipped")

# Add the aeolis package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'aeolis'))

try:
    from aeolis.model import AeoLiSRunner
    print("✓ AeoLiS imported successfully")
except ImportError as e:
    print(f"✗ Error importing AeoLiS: {e}")
    sys.exit(1)

def format_time(seconds):
    """Format seconds into readable time string"""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    elif seconds < 3600:
        return f"{seconds/60:.2f} minutes"
    else:
        return f"{seconds/3600:.2f} hours"

def get_memory_usage():
    """Get current memory usage in MB"""
    if HAS_PSUTIL:
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    else:
        return 0.0  # Return 0 if psutil not available

def run_aeolis_with_timing(config_file):
    """Run AeoLiS model with detailed timing analysis using cProfile"""
    
    print("=" * 60)
    print("AeoLiS Model Execution with cProfile Analysis")
    print("=" * 60)
    print(f"Configuration file: {config_file}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if HAS_PSUTIL:
        print(f"Initial memory usage: {get_memory_usage():.2f} MB")
    print("=" * 60)
    
    # Change to the directory containing the configuration file
    config_dir = os.path.dirname(os.path.abspath(config_file))
    original_dir = os.getcwd()
    os.chdir(config_dir)
    
    try:
        # Record start time
        start_time = time.time()
        
        # Initialize model
        print("\n📊 Initializing AeoLiS model...")
        init_start = time.time()
        model = AeoLiSRunner(configfile=os.path.basename(config_file))
        init_end = time.time()
        print(f"   ✓ Model initialization completed in {format_time(init_end - init_start)}")
        if HAS_PSUTIL:
            print(f"   Memory usage after init: {get_memory_usage():.2f} MB")
        
        # Run model with cProfile
        print("\n🚀 Running AeoLiS simulation with cProfile...")
        run_start = time.time()
        
        # Create profiler
        profiler = cProfile.Profile()
        profiler.enable()
        
        # Run the model
        model.run()
        
        # Stop profiling
        profiler.disable()
        run_end = time.time()
        
        # Calculate total time
        total_time = run_end - start_time
        
        print("\n" + "=" * 60)
        print("TIMING SUMMARY")
        print("=" * 60)
        print(f"Model initialization: {format_time(init_end - init_start)}")
        print(f"Model execution:      {format_time(run_end - run_start)}")
        print(f"Total runtime:        {format_time(total_time)}")
        if HAS_PSUTIL:
            print(f"Peak memory usage:    {get_memory_usage():.2f} MB")
        print(f"End time:             {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 60)
        
        # Performance metrics
        print("\n📈 PERFORMANCE METRICS")
        print("=" * 60)
        if hasattr(model, 'p') and hasattr(model.p, 'tstop') and hasattr(model.p, 'tstart'):
            sim_duration = model.p['tstop'] - model.p['tstart']
            if sim_duration > 0:
                sim_speed = sim_duration / total_time
                print(f"Simulation duration:  {format_time(sim_duration)} (model time)")
                print(f"Simulation speed:     {sim_speed:.2f}x real-time")
                print(f"Time step efficiency: {sim_duration/3600:.1f} model hours per wall-clock hour")
        
        # Generate profiling report
        print("\n🔍 PROFILING REPORT")
        print("=" * 60)
        
        # Save detailed profile to file
        profile_file = "aeolis_profile.prof"
        profiler.dump_stats(profile_file)
        print(f"Detailed profile saved to: {profile_file}")
        
        # Create in-memory string buffer for profile output
        s = io.StringIO()
        ps = pstats.Stats(profiler, stream=s)
        
        # Sort by cumulative time and show top functions
        ps.sort_stats('cumulative')
        ps.print_stats(20)  # Top 20 functions by cumulative time
        
        profile_output = s.getvalue()
        print("\nTop 20 functions by cumulative time:")
        print(profile_output)
        
        # Sort by internal time and show top functions
        s = io.StringIO()
        ps = pstats.Stats(profiler, stream=s)
        ps.sort_stats('tottime')
        ps.print_stats(20)  # Top 20 functions by total time
        
        profile_output = s.getvalue()
        print("\nTop 20 functions by total time:")
        print(profile_output)
        
        # Save human-readable profile report
        report_file = "aeolis_profile_report.txt"
        with open(report_file, 'w') as f:
            f.write("AeoLiS Model Profiling Report\n")
            f.write("=" * 60 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Configuration: {config_file}\n")
            f.write(f"Total runtime: {format_time(total_time)}\n\n")
            
            # Cumulative time report
            f.write("TOP FUNCTIONS BY CUMULATIVE TIME:\n")
            f.write("-" * 40 + "\n")
            ps = pstats.Stats(profiler, stream=f)
            ps.sort_stats('cumulative')
            ps.print_stats(50)
            
            # Total time report
            f.write("\n\nTOP FUNCTIONS BY TOTAL TIME:\n")
            f.write("-" * 40 + "\n")
            ps.sort_stats('tottime')
            ps.print_stats(50)
            
        print(f"\nDetailed profile report saved to: {report_file}")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Error running model: {e}")
        print(f"Error type: {type(e).__name__}")
        print("\nFull traceback:")
        traceback.print_exc()
        return False
        
    finally:
        # Return to original directory
        os.chdir(original_dir)

def main():
    """Main function"""
    config_file = r"c:\Users\svries\Documents\GitHub\OE_aeolis-python\aeolis\examples\vanWesten2024\blowout\blowout_short.txt"
    
    if not os.path.exists(config_file):
        print(f"❌ Configuration file not found: {config_file}")
        return 1
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    success = run_aeolis_with_timing(config_file)
    
    if success:
        print("\n✅ Model execution completed successfully!")
        return 0
    else:
        print("\n❌ Model execution failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
