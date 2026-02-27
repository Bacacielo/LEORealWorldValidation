# LEO Real-World Validation Engine

A high-fidelity visualization and simulation engine for validating LEO (Low Earth Orbit) satellite links under real-world conditions. This project integrates SGP4 propagation with atmospheric physics models to assess link viability within the Starlink network.

## 🚀 Key Features

* **Real-Time Simulation:** SGP4 propagation for 2000+ Starlink satellites using live TLE (Two-Line Element) data from NORAD.
* **Weather Modeling (ITU-R P.838):** Rain attenuation calculation for Ka-band (26 GHz) frequencies.
* **Link Telemetry:** Real-time computation of dynamic Path Loss, Doppler Shift, and Latency.
* **Interactive 3D Globe:** Select ground stations (London, New York, Tokyo, etc.) and visualize routing paths.
* **Timeline Control:** Playback and time-scrubbing capabilities for the simulation timeline.

## 📂 File Structure

The project consists of three main components:

1.  **`Generate_Starlink_Data_v4.py` (Python Backend)**
    * Downloads the latest Starlink TLEs from CelesTrak.
    * Uses the `skyfield` library to propagate orbits.
    * Exports satellite positions and velocities to `starlink_data.mat`.

2.  **`AppRealWorldValidation.m` (MATLAB Frontend)**
    * The main GUI application.
    * Handles visualization, user interaction, and real-time metric display.

3.  **`SimUtils2.m` (Physics Engine)**
    * Core utility class handling physical calculations.
    * Includes coordinate conversions (LLA to ECEF), Link Budget calculations, and routing algorithms (K-Shortest Paths).

## 🛠️ Prerequisites

### Python (for data generation)
* Python 3.x
* Libraries: `numpy`, `scipy`, `skyfield`
    ```bash
    pip install numpy scipy skyfield
    ```

### MATLAB (for simulation)
* MATLAB (R2020b or later recommended for `uifigure` support).
* Standard toolboxes for math and GUI components.

## ⚙️ Installation & Usage

1.  **Generate Orbit Data:**
    Run the Python script first to download satellite data and generate the required `.mat` file.
    ```bash
    python Generate_Starlink_Data_v4.py
    ```
    *Note: This step is mandatory as MATLAB will report an error if `starlink_data.mat` is missing.*

2.  **Launch Application:**
    Open MATLAB, navigate to the project folder, and run:
    ```matlab
    AppRealWorldValidation
    ```

## 📊 Physics Models & Parameters

The simulator uses specific models for realism:

* **Frequency:** 26 GHz (Ka-Band) for Inter-Satellite Links (ISL) and Ground-to-Space.
* **Rain Model:** **ITU-R P.838** with coefficients $k=0.187$ and $\alpha=1.154$.
* **Doppler:** Relative velocity and frequency shift calculation with a stability warning threshold at 100 MHz.
* **Routing:** Yen’s Algorithm for finding k-shortest paths, weighing distance and network load.

## 👥 Credits

* **Author:** chrisvasill
* **Repository:** LEO-RealWorld-Validation
