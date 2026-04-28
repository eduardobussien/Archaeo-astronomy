# Archaeo-astronomy Alignment Tool 

A Python-based digital archaeology tool designed to calculate historical stellar alignments for any terrestrial location and date. This project bridges the gap between ancient history and modern software development, allowing researchers and explorers to mathematically reconstruct the night sky as it appeared above ancient temples and monuments millennia ago.

## Features
* **Time-Travel Calculations:** Computes Julian Dates and Local Sidereal Time for deep historical dates (e.g., 10,500 BC).
* **Geospatial Precision:** Uses exact latitude and longitude coordinates to determine the altitude and azimuth of major navigational stars.
* **Targeted Stellar Database:** Built-in tracking for archaeo-astronomically significant stars like Sirius and the stars of Orion's Belt.

## Tech Stack
* **Python 3**
* **NumPy** for high-performance astronomical trigonometry (no heavy external astronomy libraries required).

## Usage
Simply define the geographical coordinates (latitude/longitude) and the historical date. The script calculates the exact positions of key stars relative to that specific location on Earth.

```python
# Example: Calculating alignments for the Great Pyramid of Giza in 10,500 BC at Midnight
r = calculate_alignments(lat=29.9792, lon=31.1342, year=-10499, month=3, day=20, hour=0)