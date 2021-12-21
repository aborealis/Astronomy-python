"""
This file explains API of Directions class.
"""
from primary_directions import Directions

# To set initial data, use
# dictionary with following fields
DATA = dict(
    year=1948,
    month=11,
    day=14,
    hour=21,
    min=14 + 39/60,
    tzone=0,
    lat=51 + 30/60,
    lon=0 - 10/60
)

# Create an instance of Direction class
dr = Directions(DATA)


# Now you can access astronomical properties
# of a particular planet/house cusp in the
# Placidus house system, for houses (houses
# indexed 0..11), or for planets (planets are
# indexed 0..9, Sun..Pluto). Full list of
# properties is given below.

print('Declination of 5th house:',
      dr.house[4].dec)
print('Ascention difference (AC) of Merciry',
      dr.planet[3].AD)

# You can create a custom object which contains
# astrolomical data of a particular ecliptical
# point with givel zodiacal longitude and
# celestial latitude

zero_degree_of_cancer = dr.get_equatorial_data(
    lon=90.0,  # 0 degree of Cancer
    lat=0.0
)
print('Semi-arc (SA) of 0° Cancer:',
      zero_degree_of_cancer.SA)

# You can also access sidereal sphere related
# properties, f.e. oblique ascention of AS.
# Full list of spherical peoperties is given
# below.

print('Oblique ascention (OA) of ASC:',
      dr.sphere.OA_asc)

# Instance of Direction class also contains
# methods to calculate mundane directions
# from the promissor to the significator

print('Moon conjuncts MC (mundi)', dr.aspect_mundi(
    promissor=dr.planet[1],  # The Moon
    significator=dr.house[9]  # 10th house cusp
))

print('Jup conjuncts DC (mundi):',  dr.aspect_mundi(
    promissor=dr.planet[5],  # Jupiter
    significator=dr.house[6]  # 7th house cusp
))

print('Sun conjuncts Merc (mundi):', dr.aspect_mundi(
    promissor=dr.planet[0],  # The Sun
    significator=dr.planet[2]  # Mercury
))

# You can specify options for mundane
# conjunctions:
# option=1: calc arc to singificator's parallel
# option=2: calc arc to significator's contrparallel
# option=3: calc arc to

print('Sat || Mer (mundi):', dr.aspect_mundi(
    promissor=dr.planet[6],
    significator=dr.planet[2],
    option=1
))

print('Pluto conjuncts VTX/AVTX (mundi):', dr.aspect_mundi(
    promissor=dr.planet[9],
    significator=None,
    option=3))

# If you specify aspect in degree (30, 60, 90...),
# you'll get according mundane aspect, f.e.

print('Mon Tri Sat (mundi):', dr.aspect_mundi(
    promissor=dr.planet[1],
    significator=dr.planet[6],
    aspect=120
))

# You can also calculate zodiac directions.
# It accepsts the same options as .aspect_mundi()
# method above

print('Sun Trin Merc (zodiaco):', dr.aspect_zodiaco(
    promissor=dr.planet[1],
    significator=dr.planet[6],
    aspect=120
))

# And you also can get field-plane directions
# It accepsts the same options as .aspect_mundi()
# method above

print('Mon Trin Sat (field plane):', dr.aspect_field_plane(
    promissor=dr.planet[1],
    significator=dr.planet[6],
    planetary_orbit=1,
    aspect=120
))


# Eventally you can print a summary table on terminal
# with major astronomical data
dr.print_coordinates()

# MAIN PROPERTIES OF SIDERIAL SPHERE
# Accessed with .sphere.PROPERTY_NAME
# RAMC: MC Right Ascension
# RAIC: IC Right Ascension
# OA_asc: ASC Oblique Ascension
# OD_dsc: DSC Oblique Descension

# MAIN PROPERTIES OF A GIVEN PLANET/HOUSE
# Accessed with .planet[N].PROPERTY_NAME
# .house[N].PROPERTY_NAME, where N = 0..9
# for planets (Sun, Moon, Mercury...Pluto),
# and N = 0..11 for 1st, 2nd..12th houses.
# lon: zodiac (ecliptical) longitude in degrees
# lat: celstial latitude (ecliptical coordinate)
# RA: Right Ascension
# dec: Declination
# AD: Ascension difference (for raising over horizon)
# OA: Oblique Ascension = RA - AD
# OD: Oblique Descension = RA + AD
# AD2: Ascension Difference type 2
#      (for crossing the VTX/AVTX)
# OA2: Oblique Ascension type 2 = RA + AD2
# OD2: Oblique Descension type 2 = RA - AD2
# DSA: Day Semi-arc 90 + AD
# NSA: Night Semi-arc 90 - AD
# UMD: Upper Meridian Distance, |RAMC - RA|
# LMD: Lower Meridian Distance, |RAIC - RA|
# MD: Meridian Distance = min(UMD,LMD)
# quadrant: true mundane position of a body.
#          1 - east & under horizon,
#          2 - west & under horizon,
#          3 - west & above horizon,
#          4 - east & above horizon.
# SA: Current Semi-arc
#          - NSA for planets below horizon
#          - DSA for planets above horizon
# O_AD:    Ascension Difference for eastern planets or
#          Ascension Difference for western planets
# O_AD2:   Ascension Difference (type 2) for eastern planets or
#          Ascension Difference (type 2) for western planets
# RA_zod:  Right Ascention of ecliptic longitude of a given body
#          (zodiac longitude), ignoring its celestial latitude
# dec_zod: Declination of zodiac longitude
# AD_zod: Ascension Difference of zodiac longitude
# OA_zod: Oblique Ascension of zodiac longitude = RA_zod - AD_zod
# OD_zod  Oblique Descension of zodiac longitude = RA_zod + AD_zod
# AD2_zod: Ascension Difference (type 2) of zodiac longitude
# OA2_zod: Oblique Ascension of zodiac longitude (type 2)
#          = RA_zod + AD2_zod
# OD2_zod: Oblique Descension of zodiac longitude (type 2)
#          = RA_zod - AD2_zod
