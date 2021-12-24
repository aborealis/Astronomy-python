from primary_directions import Directions

DATA = dict(
    year=1984,
    month=9,
    day=15,
    hour=16,
    min=20 + 00/60,
    tzone=0,
    lat=51 + 32/60,
    lon=0 - 12/60
)

dr = Directions(DATA)

Antares = dr.set_object_eqt(
    'Antares',
    r_asc=(16 + 29/60) * 15,
    decl=-(26 + 26/60)
)

Vindemiatrix = dr.set_object_eqt(
    'Vindemiatrix',
    r_asc=(13 + 2/60) * 15,
    decl=-(10 + 57/60)
)
