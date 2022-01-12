"""
A useful class to calculate primary
directions in Placidus house system.
It requires the swisseph library to
be installed: pip install pyswisseph
"""
from collections import namedtuple
from math import tan, sin, cos, asin, atan, pi
from datetime import datetime
import swisseph as swe


class Directions:
    """
    Defines methods to express coordinates
    of planets and Placidus house cusps in
    different coordinate systems. Contains
    methods to calculate primary direction
    """

    def __init__(self, data) -> None:
        self.jday = swe.julday(data["year"],
                               data["month"],
                               data["day"],
                               data["hour"] + data["min"]/60
                               ) - data["tzone"] / 24

        self.__sin_e = sin(self.__inclination_ecliptic() / 180 * pi)
        self.__cos_e = cos(self.__inclination_ecliptic() / 180 * pi)
        self.__sin_phi = sin(data['lat'] / 180 * pi)
        self.__cos_phi = cos(data['lat'] / 180 * pi)
        self.__tan_phi = tan(data['lat'] / 180 * pi)

        # A list of house cusps longitudes
        cusps = swe.houses(
            self.jday,
            data['lat'],
            data['lon'],
            bytes("P", "ascii")
        )[0]

        # Real ascension MC
        self.RAMC = self.__ra_dec_from_lon_lat(
            lon_ecl=cusps[9],
            lat_ecl=0.0
        ).horz_angle

        # Real ascension IC
        self.RAIC = self.normalize_360(self.RAMC + 180)

        # Local sidereal time
        self.LST = self.RAMC / 15

        self.planet = []
        self.__planet_names = [
            'SUN', 'MON', 'MER', 'VEN',
            'MAR', 'JUP', 'SAT', 'URN',
            'NEP', 'PLT', 'NND',
        ]
        for planet_num, planet_name in enumerate(self.__planet_names):
            # Ecliptical lon/lat
            lon, lat = swe.calc_ut(self.jday, planet_num)[0][:2]
            lat = lat if planet_num else 0.0
            self.planet.append(self.set_object(planet_name, lon, lat))

        self.house = []
        cusp_names = [f'H{house}' for house in range(1, 13)]
        for cusp_num, cusp_degree in enumerate(cusps):
            self.house.append(self.set_object(
                cusp_names[cusp_num],
                cusp_degree,
                lat_ecl=0.0
            ))

        star_names = [
            "Abhijit", "Aboras", "Abyssus Aqueus", "Acamar", "Achernar",
            "Achird", "Acrab", "Acrab", "Acrux", "Acubens", "Aculeus",
            "Acumen", "Adara", "Adhab", "Adhafera", "Adhara", "Adhil",
            "Agena", "Ahadi", "Ain", "Ain", "Ain al Rami", "Ainalhai",
            "Akrab", "Al Anz", "Al Athfar", "Al Dhanab", "Al Durajah",
            "Al Fahd", "Al Fawaris", "Al Hamalain", "Al Haud", "Al Hecka",
            "Al Jabhah", "Al Jathiyah", "Al Kalb al Rai", "Al Khabdhilinan",
            "Al Krikab", "Al Kurud", "Al Minliar al Asad", "Al Minliar al Shuja",
            "Al Mizan", "Al Nitham", "Al Pherg", "Al Rescha", "Al Rihla",
            "Al Sadr al Ketus", "Al Tarf", "Al Thalimaim Anterior",
            "Al Thalimaim Posterior", "Aladfar", "Alagemin", "Alaraph",
            "Alathfar", "Alathfar", "Albaldah", "Albali", "Albireo",
            "Albulaan", "Alchiba", "Alchita", "Alcor", "Alcyone",
            "Aldafirah", "Aldebaran", "Aldebaran", "Alderamin", "Aldhibah",
            "Aldhibain", "Alfecca Meridiana", "Alfirk", "Algedi",
            "Algedi", "Algenib", "Algieba", "Algieba", "Algol", "Algol",
            "Algorab", "Algorab", "Alhakim", "Alhena", "Alherem",
            "Alifa Al Farkadain", "Alioth", "Alkaid", "Alkalurops",
            "Alkes", "Alkibash", "Alkidr", "Alkurhah", "Almaak",
            "Almac", "Almach", "Almak", "Almeisan", "Alminhar", "Alnair",
            "Alnasl", "Alnath", "Alnilam", "Alnitak", "Alniyat",
            "Alphard", "Alphecca", "Alphekka", "Alpheratz", "Alphirk",
            "Alradif", "Alrai", "Alredif", "Alrischa", "Alsafi",
            "Alsafi", "Alsciaukat", "Alshain", "Alsharasif", "Alshat",
            "Alsuhail", "Altager", "Altair", "Altais", "Altawk",
            "Alterf", "Althaur", "Aludra", "Alula Australis", "Alula Borealis",
            "Alvahet", "Alvashak", "Alwaid", "Alya", "Alzirr", "Anazitisi",
            "Ancha", "Andromeda Galaxy", "Andromeda Galaxy", "Angetenar",
            "Ankaa", "Anser", "Antares", "Antares", "Anunitum", "Anuradha",
            "Anwar al Farkadain", "Apex", "Ara", "Arcturus", "Arkab Posterior",
            "Arkab Prior", "Armus", "Arneb", "Arrakis", "Ascella",
            "Asellus Australis", "Asellus Borealis", "Asellus Primus",
            "Asellus Secundus", "Asellus Tertius", "Ashlesha", "Ashvini",
            "Aspidiske", "Asterion", "Asterope", "Athafi", "Atik",
            "Atiks", "Atirsagne", "Atlas", "Atria", "Auva", "Avior",
            "Avis Satyra", "Azelfafage", "Azha", "Azmidiske", "Baham",
            "Barnard's star", "Baten Algiedi", "Baten Kaitos", "Batentaban Australis",
            "Batentaban Borealis", "Bazak", "Beid", "Bellatrix",
            "Benetnash", "Bered", "Beteigeuse", "Betelgeuse", "Bharani",
            "Biham", "Birdun", "Bogardus", "Bos", "Botein", "Brachium",
            "Bunda", "Bungula", "Canopus", "Capella", "Caph", "Capulus",
            "Capulus", "Castor", "Castra", "Cauda Hydrae", "Ceginus",
            "Celbalrai", "Celeano", "Cestan", "Chara", "Chertan",
            "Chin", "Chort", "Chow", "Cih", "Citra", "Cor Caroli",
            "Cor Hydrae", "Cor Serpentis", "Coxa", "Cujam", "Cursa",
            "Dabih", "Decapoda", "Decrux", "Deli", "Deneb", "Deneb",
            "Deneb Adige", "Deneb Algedi", "Deneb Algenubi", "Deneb Dulphim",
            "Deneb Kaitos", "Deneb el Okab Australis", "Deneb el Okab Borealis",
            "Denebola", "Dhanab al Shuja", "Dhanishtha", "Dheneb",
            "Dhur", "Diadem", "Difda", "Diphda", "Dorsum", "Drus",
            "Drys", "Dschubba", "Dubhe", "Dziban", "Ed Asich", "Edasich",
            "Ekkhysis", "El Kophrah", "El Nath", "Electra", "Elnath",
            "Eltanin", "Enif", "Erakis", "Errai", "Etamin", "Facies",
            "Farkadain", "Fomalhaut", "Foramen", "Fornacis", "Fudail",
            "Fum Alsamakah", "Furibundus", "Furud", "GCRS00", "Gacrux",
            "Gal. Center", "Gal. Center", "Gemma", "Ghusn al Zaitun",
            "Gianfar", "Giansar", "Giedi Prima", "Giedi Secunda",
            "Gienah Corvi", "Gienah Cygni", "Gienah Ghurab", "Girtab",
            "Gomeisa", "Gorgona Quatra", "Gorgona Secunda", "Gorgona Tertia",
            "Graffias", "Grafias", "Great Attractor", "Gruid", "Grumium",
            "Hadar", "Haedi", "Haedus", "Hamal", "Han", "Haris",
            "Hasseleh", "Hasta", "Hatsya", "Hecatebolus", "Heka",
            "Helkath", "Hemelein Prima", "Hemelein Secunda", "Heze",
            "Hilasmus", "Hoedus I", "Hoedus II", "Homam", "Hyadum I",
            "Hyadum II", "Hydor", "Hydria", "Hydrobius", "Imad",
            "Ira Furoris", "Isis", "Izar", "Jabbah", "Jabhat al Akrab",
            "Jabhat al Akrab", "Jih", "Juxta Crucem", "Jyeshtha",
            "Kabkent Secunda", "Kabkent Tertia", "Kaffaljidhma",
            "Kaht", "Kaimana", "Kajam epsHer", "Kajam omeHer", "Kakkab",
            "Kattupothu", "Kaus Australis", "Kaus Borealis", "Kaus Medis",
            "Kaus Meridionalis", "Ke Kwan", "Kebash", "Keid", "Kekouan",
            "Kelb Alrai", "Kerb", "Ketu", "Khambalia", "Kissin",
            "Kitalpha", "Kochab", "Koo She", "Kornephoros", "Kraz",
            "Krttika", "Kuma", "Kuma", "Kurdah", "Kurhah", "Labrum",
            "Leiolepidotus", "Leiolepis", "Lesath", "Linteum", "Ma Ti",
            "Maasym", "Maaz", "Mabsuthat", "Mabsuthat", "Maculata",
            "Maculosa", "Magha", "Maia", "Manubrium", "Manus", "Marakk",
            "Marfak", "Marfik", "Marfik", "Markab", "Markeb", "Marsik",
            "Masym", "Matar", "Mati", "Mautinah", "Mebsuta", "Megrez",
            "Meissa", "Mekbuda", "Melkarth", "Men", "Menkalinan",
            "Menkar", "Menkar", "Menkent", "Menkib", "Merak", "Merga",
            "Merope", "Mesarthim", "Metallah", "Miaplacidus", "Mimosa",
            "Minchir", "Minkar", "Mintaka", "Mira", "Mira", "Mirach",
            "Mirak", "Miram", "Mirfak", "Mirphak", "Mirzam", "Misam",
            "Mizar", "Mrgashirsha", "Mufrid", "Muhlifain", "Mula",
            "Muliphein", "Muphrid", "Murzim", "Murzims", "Muscida",
            "Nageba", "Nair al Saif", "Nanto", "Naos", "Nash", "Nashira",
            "Nehushtan", "Nekkar", "Nihal", "Nodus I", "Nodus II",
            "Nulla Pambu", "Nunki", "Nusakan", "Oculus", "Pazan",
            "Pazhan", "Peacock", "Peregrini", "Phact", "Phaeo", "Phaesula",
            "Pharkadain", "Phecda", "Pherkad", "Pherkad Minor", "Phicares",
            "Phicareus", "Phycochroma", "Pleione", "Pleura", "Polaris",
            "Polaris", "Polaris Australis", "Polis", "Pollux", "Porrima",
            "Praecipua", "Praesepe Cluster", "Praesepe Cluster",
            "Prijipati", "Prima Hyadum", "Princeps", "Procyon", "Propus etaGem",
            "Propus iotGem", "Proxima Centauri", "Pulcherrima", "Punarvasu",
            "Purvabhadra", "Purvaphalguni", "Purvashadha", "Pushya",
            "Pushya", "Qin", "Rana", "Ras Algethi", "Ras Alkurki",
            "Ras Elased Australis", "Ras Elased Borealis", "Ras Mutallah",
            "Rasalas", "Rasalgethi", "Rasalhague", "Rastaban", "Regor",
            "Regulus", "Regulus", "Rehla", "Revati", "Rigel", "Rigel",
            "Rigel Kentaurus", "Rigil Kent", "Rijl al Awwa", "Ril Alauva",
            "Rohini", "Rotanev", "Ruc", "Rucha", "Ruchbah", "Ruchbah I",
            "Ruchbah II", "Rukbalgethi Genubi", "Rukbalgethi Shemali",
            "Rukbat", "Rukh", "Rutilicus", "Sabik", "Sadalachbia",
            "Sadalbari", "Sadalmelek", "Sadalmelik", "Sadalsuud",
            "Sadaltager", "Sadatoni", "Sadir", "Sador", "Sadr", "Saidak",
            "Saiph", "Salm", "Samakah", "Sanduleak", "Sargas", "Sarin",
            "Sasin", "Sataghni", "Sceptrum", "Scheat", "Schedar",
            "Schedir", "Scutulum", "Seat", "Secunda Hyadum", "Segin",
            "Segin", "Seginus", "Sephdar", "Sham", "Shatabhishaj",
            "Shaula", "Shedir", "Sheliak", "Shemali", "Sheratan",
            "Shir", "Shishimai", "Shravana", "Shravishtha", "Simiram",
            "Simmah", "Sinistra", "Sirius", "Sirius", "Sirrah", "Situla",
            "Skat", "Sofian", "Spica", "Spica", "Spiculum", "Sterope I",
            "Sterope II", "Sualocin", "Subra", "Suhail", "Suhail Hadar",
            "Suhail al Muhlif", "Sulafat", "Sulaphat", "Svati", "Syrma",
            "Tabit", "Tabit", "Talitha Australis", "Talitha Borealis",
            "Tang", "Tania Australis", "Tania Borealis", "Tarazed",
            "Taygeta", "Tegmen", "Tegmine", "Tejat", "Terebellium",
            "Thabit", "The Blaze Star", "The Garnet Star", "Theemin",
            "Thuban", "Thusia", "Tien Kang", "Toliman", "Torcularis Septentrionalis",
            "Trapezium", "Tse Tseng", "Tseen Foo", "Tseen Ke", "Tsih",
            "Tsze", "Tsze Tseang", "Turais", "Tureis", "Tyl", "Ukdah",
            "Ukdah", "Unukalhai", "Urakhga", "Urodelus", "Ushakaron",
            "Uttarabhadra", "Uttaraphalguni", "Uttarashadha", "Vathorz Posterior",
            "Vathorz Prior", "Vega", "Vernalis", "Vindemiatrix",
            "Virgo Cluster", "Vishakha", "Wasat", "Wazn", "Wei",
            "Wezen", "Xestus", "Yed Posterior", "Yed Prior", "Yildun",
            "Zaniah", "Zaurak", "Zavijava", "Zhou", "Zibal", "Zosma",
            "Zuben Elakrab", "Zuben Elakribi", "Zuben Elgenubi",
            "Zuben Eschamali", "Zuben Hakrabi", "Zubenelakrab", "Zubenelakribi",
        ]
        self.stars = []
        for name in star_names:
            lon, lat = swe.fixstar_ut(
                name,
                self.jday
            )[0][:2]
            self.stars.append(self.set_object(
                name, lon, lat
            ))

        self.__geo_lat = data['lat']
        self.__geo_lon = data['lon']

    Point_Cartasian = namedtuple(
        'Cartasian', ['x', 'y', 'z']
    )

    Point_Spherical = namedtuple(
        'Spherical',
        ['horz_angle', 'vert_angle']
    )

    Celestial_Object = namedtuple(
        'Object',
        [
            'id',
            'x_ecl', 'y_ecl', 'z_ecl',  # ecliptical xyz
            'x_eqt', 'y_eqt', 'z_eqt',  # equatorial xyz
            'x_loc', 'y_loc', 'z_loc',  # local horz xyz
            'RA', 'D',  # real ascension and declination
            'lon', 'lat',  # ecliptical lon & lat
            'AD', 'AD2',  # ascension diff. type 1 and 2
            'UMD', 'LMD',  # Upper/Lower meridian dist.
            'DSA', 'NSA',  # Day/Night semi-arc
            'quadrant',  # 1,2 -under horison, 4,1 - west
            'OA',  # Oblique ascention (time to ASC)
        ]
    )

    @staticmethod
    def normalize_360(degree: float) -> float:
        """
        Normalizes degree to 360 scale
        """
        if degree < 0:
            return degree + 360
        if degree > 360:
            return degree - 360
        return degree

    @staticmethod
    def shortest_distance(deg1, deg2):
        """
        Returns the shortest arc between two angles
        """
        arc1 = max(deg1, deg2) - min(deg1, deg2)
        arc2 = 360 - arc1
        return min(arc1, arc2)

    def __inclination_ecliptic(self) -> float:
        """
        Returns the angle of ecliptic inclanation
        in degrees
        """
        julian_era = (self.jday - 2415020.0) / 36525
        return (
            23 + 27 / 60 + 8.26 / 3600 -
            46.845 / 3600 * julian_era +
            0.0059 / 3600 * julian_era ** 2 +
            0.001811 / 3600 * julian_era ** 3
        )

    def __xyz_eqt_from_xyz_ecl(self,
                               point: Point_Cartasian
                               ) -> Point_Cartasian:
        """
        Transfers XYZ from ecliptic to XYZ equator plane
        x_eqt directed to the vernal equinox
        y_eqt directed to the projection of 0° Cancer on eqt
        z_eqt directed to the North Pole
        """
        x_eqt = point.x
        y_eqt = self.__cos_e * point.y - self.__sin_e * point.z
        z_eqt = self.__sin_e * point.y + self.__cos_e * point.z

        return self.Point_Cartasian(x_eqt, y_eqt, z_eqt)

    def __xyz_ecl_from_xyz_eqt(self,
                               point: Point_Cartasian
                               ) -> Point_Cartasian:
        """
        Transfers XYZ from equator to XYZ ecliptic plane.
        x_ecl directed to 0° Aries
        y_ecl directed to 0° Cancer
        z_ecl directed to northen eqt. hemisphere
        """
        x_ect = point.x
        y_ect = self.__cos_e * point.y + self.__sin_e * point.z
        z_ect = - self.__sin_e * point.y + self.__cos_e * point.z

        return self.Point_Cartasian(x_ect, y_ect, z_ect)

    def __xyz_loc_from_xyz_eqt(self,
                               point: Point_Cartasian
                               ) -> Point_Cartasian:
        """
        Transfers XYZ from equator to XYZ local horizon.
        x_loc directed to East Point
        y_loc directed to North Point
        z_loc directed to Azimuth
        """
        # time since Aries raising (tsa, in radians)
        tsa = (self.LST * 15 + 90) / 180 * pi

        sin_t = sin(tsa)
        cos_t = cos(tsa)

        sin_p = self.__sin_phi
        cos_p = self.__cos_phi

        x_loc = cos_t * point.x + sin_t * point.y
        y_loc = (- sin_p * sin_t * point.x
                 + sin_p * cos_t * point.y
                 + cos_p * point.z)
        z_loc = (cos_p * sin_t * point.x -
                 cos_p * cos_t * point.y +
                 sin_p * point.z)

        return self.Point_Cartasian(x_loc, y_loc, z_loc)

    def __xyz_from_lat_lon(self,
                           lon_ecl: float,
                           lat_ecl: float) -> Point_Cartasian:
        """
        Transfers XYZ from a given plane to
        celestial lon & lat of the same plane.
        Ecliptical latitudes are positive for
        northern hemisphere
        """
        cos_lat = cos(lat_ecl / 180 * pi)
        sin_lat = sin(lat_ecl / 180 * pi)
        cos_lon = cos(lon_ecl / 180 * pi)
        sin_lon = sin(lon_ecl / 180 * pi)

        x_ecl = cos_lat * cos_lon
        y_ecl = cos_lat * sin_lon
        z_ecl = sin_lat

        return self.Point_Cartasian(x_ecl, y_ecl, z_ecl)

    def __lon_lat_from_xyz(self,
                           point: Point_Cartasian
                           ) -> Point_Spherical:
        """
        Transfers XYZ from a given plane to lon
        and lat of the same plane (in degrees)
        """
        lon = atan(point.y / point.x) / pi * 180
        if point.x < 0 and point.y > 0:
            lon += 180
        elif point.x < 0 and point.y < 0:
            lon += 180
        elif point.x > 0 and point.y < 0:
            lon += 360
        lat = asin(point.z) / pi * 180

        return self.Point_Spherical(horz_angle=lon,
                                    vert_angle=lat)

    def __ra_dec_from_lon_lat(self,
                              lon_ecl: float,
                              lat_ecl: float) -> Point_Spherical:
        """
        Transfers lon & lat from ecliptic plane to real
        ascension & decl in equator plane (in degrees)
        """
        return self.__lon_lat_from_xyz(
            self.__xyz_eqt_from_xyz_ecl(
                self.__xyz_from_lat_lon(lon_ecl, lat_ecl)
            )
        )

    def __lon_lat_from_ra_dec(self,
                              r_asc: float,
                              decl: float) -> Point_Spherical:
        """
        Transfers real ascension & dec from eqt, plane
        lon & lat in ecliptic plane (in degrees)
        """
        return self.__lon_lat_from_xyz(
            self.__xyz_ecl_from_xyz_eqt(
                self.__xyz_from_lat_lon(r_asc, decl)
            )
        )

    def __ad_from_dec(self, decl: float) -> float:
        """
        Returns ascension difference (in degrees)
        of horizon raising of celestial body from
        its declination. AD = on geo latitude 0
        """
        multiplication = self.__tan_phi * tan(decl/180*pi)

        # The object never raises over horizon
        if abs(multiplication) > 1:
            return None

        return asin(self.__tan_phi *
                    tan(decl/180*pi)) / pi * 180

    def __ad2_from_dec(self, decl: float) -> float:
        """
        Returns ascension difference (in degrees)
        of celestial body crossing prime vertical
        from the body's declination. AD2 = 0 on
        geo latitude 0 (earth equator)
        """
        ratio = tan(decl/180*pi) / self.__tan_phi

        # Mathematically abs(ratio) cannot exceed 1,
        # but computationally it can by small value.
        # We need to avoid such an unsertainty:
        if ratio < -1:
            ratio = -1
        elif ratio > 1:
            ratio = 1

        return asin(ratio) / pi * 180

    def set_object_eqt(self,
                       name: str,
                       r_asc: float,
                       decl: float) -> Celestial_Object:
        """
        Create Celestial Object from
        equatorial coordinates
        """
        lon_ecl, lat_ecl = self.__lon_lat_from_ra_dec(r_asc, decl)

        return self.set_object(name,
                               self.normalize_360(lon_ecl),
                               self.normalize_360(lat_ecl))

    def set_object(self,
                   name: str,
                   lon_ecl: float,
                   lat_ecl: float) -> Celestial_Object:
        """
        Create Celestial Object from
        ecliptical coordinates
        """
        # See Celestial_Object named tuple
        # definition above for reference on
        # astronomical data prepared here
        xyz_ecl = self.__xyz_from_lat_lon(lon_ecl, lat_ecl)
        xyz_eqt = self.__xyz_eqt_from_xyz_ecl(xyz_ecl)
        xyz_loc = self.__xyz_loc_from_xyz_eqt(xyz_eqt)
        r_asc, decl = self.__lon_lat_from_xyz(xyz_eqt)
        asc_diff = self.__ad_from_dec(decl)

        upper_md = self.shortest_distance(self.RAMC, r_asc)
        lower_md = self.shortest_distance(self.RAIC, r_asc)
        upper_md = upper_md if upper_md < 180 else 360 - upper_md
        lower_md = lower_md if lower_md < 180 else 360 - upper_md

        day_sa = 90 + asc_diff if isinstance(asc_diff, float) else None
        night_sa = 90 - asc_diff if isinstance(asc_diff, float) else None

        # if body can raise over horizon
        if asc_diff:

            # If celestial object above or below the horizon
            if upper_md > day_sa or lower_md < night_sa:
                above_horizon = False
            else:
                above_horizon = True

            # To the east or west
            east_point_ra = self.normalize_360(self.RAMC + 90)
            if self.shortest_distance(east_point_ra, r_asc) < 90:
                to_the_east = True
            else:
                to_the_east = False

        # if body never raise or descend
        else:
            above_horizon = True if xyz_loc.z > 0 else False
            to_the_east = True if xyz_loc.x > 0 else False

        if above_horizon:
            quadrant = 4 if to_the_east else 3
        elif not above_horizon:
            quadrant = 1 if to_the_east else 2

        return self.Celestial_Object(
            id=name,
            x_ecl=xyz_ecl.x,
            y_ecl=xyz_ecl.y,
            z_ecl=xyz_ecl.z,
            x_eqt=xyz_ecl.x,
            y_eqt=xyz_eqt.y,
            z_eqt=xyz_eqt.z,
            x_loc=xyz_loc.x,
            y_loc=xyz_loc.y,
            z_loc=xyz_loc.z,
            RA=r_asc,
            D=decl,
            lon=lon_ecl,
            lat=lat_ecl,
            AD=asc_diff,
            AD2=self.__ad2_from_dec(decl),
            UMD=upper_md,
            LMD=lower_md,
            DSA=day_sa,
            NSA=night_sa,
            quadrant=quadrant,
            OA=r_asc - asc_diff if asc_diff else None,
        )

    def create_vertex(self,
                      obj: Celestial_Object) -> Celestial_Object:
        """
        Creates Celestal object for vertex point
        (cross of day circle of a given celestial
        object with western part of prime meridian)
        """
        decl = obj.D
        r_asc = self.RAMC - 90 + obj.AD2
        return self.set_object_eqt('VTX', r_asc, decl)

    def create_antivertex(self,
                          obj: Celestial_Object) -> Celestial_Object:
        """
        Creates Celestal object for a-vertex point
        (cross of day circle of a given celestial
        object with eastern part of prime meridian)
        """

        decl = obj.D
        r_asc = self.RAMC + 90 - obj.AD2
        return self.set_object_eqt('VTX', r_asc, decl)

    def create_parallel(self,
                        obj: Celestial_Object) -> Celestial_Object:
        """
        Creates Celestal object paralleled
        to a given one in Placidua system
        """
        decl = obj.D
        if obj.quadrant in [4, 1]:
            r_asc = self.normalize_360(self.RAMC - obj.UMD)
        else:
            r_asc = self.normalize_360(self.RAMC + obj.UMD)
        return self.set_object_eqt(f'||{obj.id}', r_asc, decl)

    def create_cont_parallel(self,
                             obj: Celestial_Object) -> Celestial_Object:
        """
        Creates Celestal object which is in contra
        parallel with a given one in Placidua system
        """
        decl = obj.D
        if obj.quadrant == 4:
            r_asc = self.RAMC + obj.LMD
        elif obj.quadrant == 3:
            r_asc = self.RAMC - obj.LMD
        elif obj.quadrant == 1:
            r_asc = self.RAIC - obj.UMD
        else:
            r_asc = self.RAIC + obj.LMD

        return self.set_object_eqt(f'-||{obj.id}', r_asc, decl)

    def conjunction_placidus(self,
                             obj1: Celestial_Object,
                             obj2: Celestial_Object) -> float:
        """
        Calculates the distance in degrees that
        the 1st object needs to pass to conjunct
        the 2nd object according Placidus method
        """
        # If object never raises or descends this
        # method is not applicable
        if (
            not isinstance(obj1.AD, float) or
            not isinstance(obj2.AD, float)
        ):
            return None

        if abs(obj1.lon - obj2.lon) < 1e-10:
            return 0

        if obj2.quadrant == 4:
            ratio = obj2.UMD / obj2.DSA
            target_ra = self.RAMC + obj1.DSA * ratio
        elif obj2.quadrant == 1:
            ratio = obj2.LMD / obj2.NSA
            target_ra = self.RAIC - obj1.NSA * ratio
        elif obj2.quadrant == 2:
            ratio = obj2.LMD / obj2.NSA
            target_ra = self.RAIC + obj1.NSA * ratio
        else:
            ratio = obj2.UMD / obj2.DSA
            target_ra = self.RAMC - obj1.DSA * ratio

        # It can happen that conjunctions to 1st or 7th
        # house may randomly be calculated by formulas
        # for object 2, located above/beyond horizon. It
        # happens due to computational uncertainty where
        # horizon object with MD almost = SA located -
        # beyond or above horizon. As a result the longest
        # path of direction may be applied. We need to take
        # it into account.
        arc = obj1.RA - target_ra
        if abs(arc) > 180:
            arc = arc - 360 if arc > 0 else arc + 360

        return arc

    def __get_aspect_point_zodiaco(self,
                                   promissor: Celestial_Object,
                                   aspect: int,
                                   lat_ecl: float) -> float:
        """
        Returns aspect point for zodiacal directions
        with the given celestial (ecliptic) latitude
        """
        return self.set_object(
            'Aspect point',
            lon_ecl=promissor.lon + aspect,
            lat_ecl=lat_ecl
        )

    def aspect_zodiaco(self,
                       promissor: Celestial_Object,
                       significator: Celestial_Object,
                       aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor's
        aspect point should pass to meets the significator
        in zodiac direction algorythm (with latitude 0)
        """
        # if object never raises or descends
        # this method is not applicable
        if not isinstance(significator.AD, float):
            return None

        aspect_point = self.__get_aspect_point_zodiaco(
            promissor,
            aspect,
            lat_ecl=0.0
        )

        return self.conjunction_placidus(aspect_point, significator)

    def aspect_zodiaco_lat(self,
                           promissor: Celestial_Object,
                           significator: Celestial_Object,
                           aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor's
        aspect point should pass to meets the significator
        in zodiac direction algorythm (with significator's
        latitude)
        """
        # if object never raises or descends
        # this method is not applicable
        if not isinstance(significator.AD, float):
            return None

        aspect_point = self.__get_aspect_point_zodiaco(
            promissor,
            aspect,
            lat_ecl=significator.lat
        )

        return self.conjunction_placidus(aspect_point, significator)

    def aspect_zodiaco_byankini(self,
                                promissor: Celestial_Object,
                                significator: Celestial_Object,
                                aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor's
        aspect point should pass to meets the significator
        in zodiac direction algorythm (with significator's
        latitude calculated according to Byankini formula)
        """
        # if object never raises or descends
        # this method is not applicable
        if not isinstance(significator.AD, float):
            return None

        lat = asin(sin(significator.lat/180*pi) *
                   cos(aspect/180*pi)) / pi * 180

        aspect_point = self.__get_aspect_point_zodiaco(
            promissor,
            aspect,
            lat_ecl=lat
        )

        return self.conjunction_placidus(aspect_point, significator)

    def aspect_mundi(self,
                     promissor: Celestial_Object,
                     significator: Celestial_Object,
                     aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor
        should pass to meet the significator's aspect
        point in proportional semi-arcs algorythm. Aspect
        is the point on the significator's hourly circle.
        """
        # if object never raises or descends
        # this method is not applicable
        if not isinstance(significator.AD, float):
            return None

        aspect_point = self.set_object_eqt(
            'Aspect point',
            self.normalize_360(significator.RA - aspect),
            significator.D
        )

        # Significator conjuncts this point when
        # it deviates from horizon to the same ratio
        # of his day/night semi-arc, as the aspect
        # point does on the significator's hour circle
        return self.conjunction_placidus(
            promissor,
            aspect_point
        )

    def aspect_mundi2(self,
                      promissor: Celestial_Object,
                      significator: Celestial_Object,
                      aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor
        should pass to meet the significator's aspect
        point in proportional semi-arcs algorythm. Aspect
        is the point on the equatorial's hourly circle in
        conjucnction with the significator
        """
        # if object never raises or descends
        # this method is not applicable
        if not isinstance(significator.AD, float):
            return None

        ratio = (significator.UMD / significator.DSA
                 if significator.quadrant in [4, 3] else
                 significator.LMD / significator.NSA)

        # The RA of the point on the equator
        # which "conjuncts" significator by
        # semi-arc proportions
        if significator.quadrant == 1:
            r_asc = self.RAIC - ratio * 90
        elif significator.quadrant == 2:
            r_asc = self.RAIC + ratio * 90
        elif significator.quadrant == 3:
            r_asc = self.RAMC - ratio * 90
        else:
            r_asc = self.RAMC + ratio * 90

        aspect_point = self.set_object_eqt(
            'Aspect point',
            r_asc=self.normalize_360(r_asc - aspect),
            decl=0
        )

        # Significator conjuncts this point when
        # it deviates from horizon to the same ratio
        # of his day/night semi-arc, as the aspect
        # point does on the equator
        return self.conjunction_placidus(
            promissor,
            aspect_point
        )

    def aspect_alexandrian(self,
                           promissor: Celestial_Object,
                           significator: Celestial_Object,
                           aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor's
        oblique ascention pass to meet the significator's
        oblique ascention (Alexandrian algorithm)
        """
        # If significator never raises ot descends
        # this method is not applicable
        if (
            not isinstance(promissor.AD, float) or
            not isinstance(significator.AD, float)
        ):
            return None

        # The point which will raise over horizon
        # after significator's ascension in N hours.
        # 8 hrs = Trine
        aspect_point = self.set_object_eqt(
            'Oblique Asc - aspect',
            # OA, time to asc. minus aspect
            significator.OA - aspect,
            significator.D
        )

        # The promissor will be in exact
        # mundane aspect with significator
        # when (significator's OA - promissot's
        # OA) = aspect, i.e. aspect point's RA
        # equalt ro promissor's OA. The distance
        # to that point is the mundane direction

        return promissor.OA - aspect_point.RA

    def aspect_field_plane(self,
                           promissor: Celestial_Object,
                           significator: Celestial_Object,
                           aspect: int) -> float:
        """
        Calculates the distance in degrees the promissor's
        aspect point should pass to meet the significator
        in fileld plane algorythm
        """
        try:
            planet = self.__planet_names.index(promissor.id)
        except ValueError:
            lon = self.normalize_360(promissor.lon + aspect)
            lat = promissor.lat
            aspect_point = self.set_object('Aspect point', lon, lat)
            return self.conjunction_placidus(aspect_point, significator)

        target_lon = self.normalize_360(self.planet[planet].lon + aspect)
        if aspect == 0:
            day = 0
        else:
            day = self.__time_to_travel(
                planet,
                target_lon=target_lon
            )

        lon, lat = swe.calc_ut(self.jday + day, planet)[0][:2]
        aspect_point = self.set_object('Aspect point', lon, lat)

        return self.conjunction_placidus(aspect_point, significator)

    def __time_to_travel(self,
                         planet: int,
                         target_lon: float) -> float:
        """
        Returns a number of julian days nesessary
        for planet (0..9 - Sun..Pluto) to pass to
        a certain zodiac degree. May return negative
        days in case the target longitude is less
        than starting point
        """
        # A very approximate number of days needed for
        # a planet in average to pass 30 degree distance.
        step = [30, 2.5, 30, 30, 80, 400, 1200, 3000, 6000, 7000][planet]
        _target_lon = self.normalize_360(target_lon)

        if _target_lon < swe.calc_ut(self.jday, planet)[0][0]:
            step = -step

        day = 0
        while True:
            while True:
                lon = swe.calc_ut(self.jday + day, planet)[0][0]
                next_lon = swe.calc_ut(self.jday + day + step, planet)[0][0]
                # For dev tracking
                # print('{:.2f} {:.2f} {:.2f}'.format(next_lon, target_lon, lon))
                if abs(lon - next_lon) > 180:
                    break
                if (_target_lon - lon) * (next_lon - _target_lon) > 0:
                    break
                day += step
            step /= 2
            if abs(_target_lon - lon) < 1e-3:
                break

        return day

    def __time_to_travel_sun(self, target_lon: float) -> float:
        """
        Returns the number of days the Sun needs to
        reach a certain longitute on ecliptic plane
        Always returns positive number
        """
        sun_lon = self.planet[0].lon

        # Shortest time may be negative
        # if shortest direction from the Sun
        # to target degree is a backword dir-n
        shortest_time = self.__time_to_travel(
            planet=0,
            target_lon=target_lon)

        # If normal direction
        if target_lon > self.planet[0].lon:
            return shortest_time

        # Retro movement
        time_to_360 = self.__time_to_travel(
            planet=0, target_lon=359.99)
        time_0_target = self.__time_to_travel(
            planet=0,
            target_lon=sun_lon + target_lon)

        return time_to_360 + abs(time_0_target)

    @staticmethod
    def get_years_ptolemey(arc: float) -> float:
        """
        According Prolemey each degree of the arc
        corresponds to 1 year of life. Defined for
        reference.
        """
        return arc

    @staticmethod
    def get_years_naibod(arc: float) -> float:
        """
        Transfers the arc of direction into the year
        according to Naibod method
        """
        return arc * 1.0146

    def get_years_simmonite(self, arc: float) -> float:
        """
        Transfers the arc of direction into the year
        according to Simmonite method
        """
        next_day_sun_ra = self.__ra_dec_from_lon_lat(
            lon_ecl=(swe.calc_ut(self.jday + 1, 0)[0][0]),
            lat_ecl=0.0
        ).horz_angle

        correction = abs(self.planet[0].RA - next_day_sun_ra)

        return arc / correction

    def get_years_placidus(self, arc: float) -> float:
        """
        Transfers the arc of direction into the year
        according to to True Solar Arc algorithm
        """

        r_asc = self.planet[0].RA + arc
        lon = self.__lon_lat_from_ra_dec(r_asc, decl=0.0).horz_angle
        lon = self.normalize_360(lon)

        return self.__time_to_travel_sun(target_lon=lon)

    def get_years_ascendant_arc(self, arc: float) -> float:
        """
        Transfers the arc of direction into the year
        according to ascendant arc algorithm
        """
        progressed_asc = swe.houses_armc(
            self.normalize_360(self.RAMC + arc),
            self.__geo_lat,
            self.__inclination_ecliptic()
        )[0][0]
        curr_asc = self.house[0].lon

        sun_lon = self.planet[0].lon
        target_lon = progressed_asc - curr_asc + sun_lon
        target_lon = self.normalize_360(target_lon)

        return self.__time_to_travel_sun(target_lon)

    def get_years_vertical_arc(self, arc: float) -> float:
        """
        Transfers the arc of direction into the year
        according to vertical arc measure
        """
        # Longitude of zodiac which crosses
        # prime vertical at the moment (y_loc
        # = 0 in local system)
        current_v = swe.houses_armc(
            self.RAIC,
            90 - self.__geo_lat,
            self.__inclination_ecliptic()
        )[0][0]

        progressed_v = swe.houses_armc(
            self.normalize_360(self.RAIC + arc),
            90 - self.__geo_lat,
            self.__inclination_ecliptic()
        )[0][0]

        sun_lon = self.planet[0].lon
        target_lon = progressed_v - current_v + sun_lon
        target_lon = self.normalize_360(target_lon)

        return self.__time_to_travel_sun(target_lon)

    def progressed_asc(self, degrees: float) -> float:
        """
        Returns progressed asc - an ecliptic longitude
        which raises after the celestial sphere rotates
        at a given number of degrees
        """
        return swe.houses_armc(
            self.normalize_360(self.RAMC + degrees),
            self.__geo_lat,
            self.__inclination_ecliptic()
        )[0][0]

    def progressed_mc(self, degrees: float) -> float:
        """
        Returns progressed MC - an ecliptic longitude
        which culminates after the celestial sphere rotates
        at a given number of degrees
        """
        progressed_jday = self.jday + (degrees / 15 / 24) * 365.25/366.25

        # A list of progressed house cusps longitudes
        cusps = swe.houses(
            progressed_jday,
            self.__geo_lat,
            self.__geo_lon,
            bytes("P", "ascii")
        )[0]

        return cusps[9]

    def show_date(self, age: float) -> str:
        """
        Transfers the number of years derived
        from into directions into a future date
        and prints result into a string
        """
        date = swe.jdut1_to_utc(self.jday + age * 365.24, 1)

        return datetime(
            date[0], date[1], date[2], date[3], date[4], int(date[5])
        ).strftime("%Y %B, %d")
