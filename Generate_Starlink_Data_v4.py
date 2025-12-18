# Generate_Starlink_Data_v4.py (FINAL FIX)
import numpy as np
import os
import scipy.io as sio
from skyfield.api import load, wgs84
from skyfield.framelib import itrs # This is ECEF

def generate_data():
    print("--- 1. Cleaning Cache ---")
    # Force clean start
    for f in ['starlink_fresh.txt']:
        if os.path.exists(f):
            try:
                os.remove(f)
            except: pass

    print("--- 2. Downloading Starlink TLEs ---")
    url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle'
    try:
        sats = load.tle_file(url, filename='starlink_fresh.txt')
    except Exception as e:
        print(f"CRITICAL: Download failed. {e}")
        return

    # Select NEW satellites (to avoid decayed ones)
    sats = sats[-1000:]
    num_sats = len(sats)
    print(f"   Processing {num_sats} satellites.")

    print("--- 3. Propagating Orbits (SGP4) ---")
    # Force download of Earth rotation data (needed for ECEF conversion)
    try:
        ts = load.timescale(builtin=True)
    except:
        ts = load.timescale()
        
    duration_min = 60
    step_sec = 10
    num_steps = int((duration_min * 60) / step_sec)
    t0 = ts.now()

    pos_data = np.zeros((num_steps, num_sats, 3))
    vel_data = np.zeros((num_steps, num_sats, 3))

    success_count = 0

    for i in range(num_steps):
        offset_sec = i * step_sec
        t_current = ts.utc(t0.utc.year, t0.utc.month, t0.utc.day, 
                           t0.utc.hour, t0.utc.minute, t0.utc.second + offset_sec)
        
        if i % 20 == 0:
            print(f"   > Step {i}/{num_steps}...")

        for k, sat in enumerate(sats):
            try:
                geom = sat.at(t_current)
                
                # --- THE FIX IS HERE ---
                # We use 'frame_xyz_and_velocity' to get both vectors at once
                p_obj, v_obj = geom.frame_xyz_and_velocity(itrs)
                
                pos_data[i, k, :] = p_obj.km
                vel_data[i, k, :] = v_obj.km_per_s
                
                if i == 0: success_count += 1
            except Exception as e:
                if i == 0 and k == 0:
                    print(f"   ! Error on sat {k}: {e}")
                pos_data[i, k, :] = [np.nan, np.nan, np.nan]

    print("--- 4. Saving Data ---")
    if success_count == 0:
        print("FATAL: Zero satellites propagated.")
    else:
        print(f"Success! Valid data for {success_count} satellites.")
        data = {
            'pos': pos_data, 
            'vel': vel_data, 
            'dt': step_sec, 
            'N': num_sats,
            'steps': num_steps
        }
        sio.savemat('starlink_data.mat', data)
        print("Done! 'starlink_data.mat' is ready for MATLAB.")

if __name__ == "__main__":
    generate_data()