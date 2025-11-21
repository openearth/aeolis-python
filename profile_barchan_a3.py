"""
Profile AeoLiS run with barchan_b3.txt configuration.

This script runs cProfile on the AeoLiS model and saves the profiling results.
"""

import cProfile
import pstats
import os
from datetime import datetime
from aeolis.model import AeoLiSRunner

# Configuration file path
config_file = r"aeolis\examples\vanWesten2024\barchan_rotate\barchan_b3.txt"

# Output profile file with timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
profile_output = f"aeolis_profile_barchan_b3_{timestamp}.prof"
print(f"Starting profiling of AeoLiS with configuration: {config_file}")
print(f"Profile will be saved to: {profile_output}")
print("-" * 80)

# Create and run the profiler
profiler = cProfile.Profile()
profiler.enable()

try:
    # Run the model
    model = AeoLiSRunner(configfile=config_file)
    model.run()
finally:
    profiler.disable()

# Save profile data
profiler.dump_stats(profile_output)
print("-" * 80)
print(f"Profiling complete! Profile saved to: {profile_output}")

# Print top 30 time-consuming functions
print("\n" + "=" * 80)
print("TOP 30 TIME-CONSUMING FUNCTIONS")
print("=" * 80)
stats = pstats.Stats(profile_output)
stats.strip_dirs()
stats.sort_stats('cumulative')
stats.print_stats(30)

print("\n" + "=" * 80)
print("TOP 30 FUNCTIONS BY TOTAL TIME")
print("=" * 80)
stats.sort_stats('tottime')
stats.print_stats(30)

print(f"\n\nTo analyze the profile further, you can use:")
print(f"  python -m pstats {profile_output}")
print("\nOr load it in a Python script:")
print(f"  import pstats")
print(f"  stats = pstats.Stats('{profile_output}')")
print(f"  stats.sort_stats('cumulative').print_stats(50)")
