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
      dr.house[4].D)
print('Ascention difference (AD) of Merciry',
      dr.planet[3].AD)

# You can create a custom object which contains
# astrolomical data of a particular ecliptical
# point with givel zodiacal longitude and
# celestial latitude

zero_degree_of_cancer = dr.set_object(
    '0° Cancer',
    lon_ecl=90.0,  # 0 degree of Cancer
    lat_ecl=0.0
)
print('Night Semi-arc (NSA) of 0° Cancer:',
      zero_degree_of_cancer.NSA)

# Alternatively you may create an object from
# equatorial coordinates (real ascen. & decl)

custom_obj = dr.set_object_eqt(
    'My Favorite Star',
    r_asc=125.3,
    dec=-12.4
)

print('Ecliptical coordinate of',
      custom_obj.id, 'is',
      custom_obj.lon
      )

# You can also access RAMC and RAIC property.

print('Right ascension of Medium Coeli:',
      dr.RAMC)

# Instance of Direction class also contains
# methods to calculate mundane directions
# from the promissor to the significator

print('Moon conjuncts MC (mundi)', dr.aspect_mundi(
    promissor=dr.planet[1],  # The Moon
    significator=dr.house[9],  # 10th house cusp
    aspect=120
))

parallel_of_mercury = dr.create_parallel(dr.planet[2])
print(
    'Sat || Mer:',
    dr.conjunction_placidus(
        dr.planet[6],  # Saturn, 6th planet
        parallel_of_mercury
    )
)

vertex = dr.create_vertex(dr.planet[5])
print(
    'Jup-VTX:',
    dr.conjunction_placidus(
        dr.planet[5], vertex
    )
)

print(
    'Mo Tri Sat (Field Plane direction):',
    dr.aspect_field_plane(
        dr.planet[1],
        dr.planet[6],
        aspect=120
    )
)

print(
    'Mo Tri Sat (Zodiac direction):',
    dr.aspect_zodiaco(
        dr.planet[1],
        dr.planet[6],
        aspect=120
    )
)


# Property of the particular planet/house:
# lon: zodiac (ecliptical) longitude in degrees
# lat: celstial latitude (ecliptical coordinate)
# RA: Right Ascension
# D: Declination
# AD: Ascension difference (for raising over horizon)
# AD2: Ascension Difference type 2
#      (for crossing the VTX/AVTX)
# DSA: Day Semi-arc 90 + AD
# NSA: Night Semi-arc 90 - AD
# UMD: Upper Meridian Distance, |RAMC - RA|
# LMD: Lower Meridian Distance, |RAIC - RA|
# quadrant: true mundane position of a body.
#          1 - east & under horizon,
#          2 - west & under horizon,
#          3 - west & above horizon,
#          4 - east & above horizon.

# Public methods:
# .shortest_distance(deg1, deg2) - Returns
#      the shortest arc between two angles

# .set_object_eqt(name, right ascention: float
#      declination: float). Returns an object
#      with celestial properties (see above)
# .set_obj() - does the same based on a given
#      longitude and latitude in the ecliptic plane
# .create_vertex() create_antivertex() - creates
#      a celestial object of the cross-point of a
#      given body with prime meridian on the west
#      or east. Requires celestial object as argument
#
# EX:
#      dr = Directions(...)
#      my_custom_object = dr.set_obj(...)
#      .create_certex(dr.house[1])
#      .create_vertex(my_cistom_object)
#
# .create_parallel(), .create_cont_parallel() creates
#      mundane parallel for a give object. The object
#      is set up as in example above.
# .conjunction_placidus(obj1, obj2) - calculates
#      directions for seni-arc proportional conjunction
#      between to celestial objects. It is equivalent
#      to .aspect_mundi() with aspect=0
# .aspect_zodiaco(), .aspect_mundi(), aspect_field_plane()
#      Calculates the distance (degree) the significator
#      should pass till meets aspect point of promissor
#      Requires clestial object1 (promissor) and celestial
#      object2 (significator)
# .get_years_ptolemey(arc), get_years_simmonite(arc)
# .get_years_naibod(arc), .get_years_placidus(arc)
# .get_years_ascendant_arc(arc), get_years_vertical_arc(arc)
#      Transfer arc of direction into years accorsing
#      different methods.
