# Python packages for astronomers
## 1. Primary direction (ancient astronomy)

A valuable class to calculate primary directions in Placidus house system. Requires swisseph library to be installed (`pip[env] install pyswisseph`)

### 1.1. Basic usage
#### 1.1.1. Import and Declaration

Declare an object with the local time and position of the observer

```python
from primary_directions import Directions

dr = Directions({
    'year': 1948,
    'month': 11,
    'day': 14,
    'hour': 21,
    'min': 14.65,
    'tzone': 0,
    # Geoposition
    'lat': 51.5, # positive for northern latitudes
    'lon': -0.6  # negative for western longitudes
})
```
#### 1.1.2. Access Astronomical Data

You can now access astronomical data of any particular planet from the list:
[Sun, Mon, Mer, Ven, Mar, Jup, Sat, Urn, Nep, Plt], for instance:

```python
dr.planet[1].AD  # Ascension difference of the Moon
dr.planet[2].OA  # Oblique Ascension of Mercury
dr.planet[0].lat # Celestial latitude of the Sun
...
```

You can pick up any point on the sky with particular ecliptic coordinates and access its astronomical data, for instance:

```python

# Let's pick up the point which corresponds
# to 0° of Leo zodiac sign
leo = dr.get_equatorial_data(
    lon=120.0,  # 0 degree of Leo
    lat=0.0
)

leo.RA # Right Ascension of 0°Leo
leo.MD # Meridian Distance of 0°Leo
...
```

A complete list of properties is described in `directions_API_explained.py` file

#### 1.1.3. Primary directions

You can also easily access methods to calculate primary directions:

```python
# Mundane direction of the Moon to square Ascendant
dr.aspect_mundi(
    promissor=dr.planet[1], # The 1st planet in the list (the Moon)
    significator=dr.house[0], # the 0th house in [0..11] (Ascendent)
    aspect=90 # Square aspect
)

# Zodiac direction of the Moon to square Ascendant
dr.aspect_zodiaco(
    promissor=dr.planet[1],
    significator=dr.house[0],
    aspect=90
)

# Arc of progressed Moon trines Saturn in a field-plane direction
aspect_field_plane(
    promissor=dr.planet[1], 
    significator=dr.planet[6], 
    planetary_orbit=1,
    aspect=120
)
```

A complete list of extra options for primary direction methods is described in `directions_API_explained.py` file

#### 1.1.3. Astrological Calculations

Finally, you may translate the arc of directions to the year of life (start point is the date, used to declare the object).

```python
dr.get_years_naibod(33.54)
# gives 34.03 years
```

You may add this time delta to initial data to see the timing:

```python
dr.years_to_date(34.03)
# gives sting output '1982, November'
```
Full list of methods to express arc in yesrs is described in `directions_API_explained.py` file