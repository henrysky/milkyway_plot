"""
Warm-up simbad query due to astroquery#3204 if a new IP address is used
"""

from astroquery.simbad import Simbad


simbad = Simbad()
try:
    simbad.query_objects(['LMC', 'M31'])
except Exception as e:
    pass