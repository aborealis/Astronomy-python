"""
Домашний калькулятор первичных дирекций
"""

from math import tan, sin, cos, asin, atan, pi
from collections import namedtuple
from datetime import datetime
import swisseph as swe


class Directions:
    """
    Создает объект промиссора и сигнификатора с атрибутами,
    содержащими с эклиптические данные для расчета дирекций.
    Включает методы расчета 1) мунданных, 2) зодиакальных и
    3) орбитальных дирекций в системе домов Плацидуса.
    """

    # Астрономические данные для отдельных планет и куспидов домов
    Coords = namedtuple('Coords', [
        'lon',  # эклиптическая долгота
        'lat',  # небесная широта (в эклиптической системе)
        'max_ecl_speed',  # средняя эклиптическая скорость (для планет)
        'RA',   # Прямое восхождение
        'dec',  # Склонение
        'AD',  # Разница восхождений до горизонта (на экваторе и на широте)
        'OA',  # Косое восхождение = RA - AD
        'OD',  # Косое захождение = RA + AD
        'AD2',  # Разница восхождений до вертекса (на экваторe и на широте)
        'OA2',  # Косое восхождение тип 2 = RA + AD2
        'OD2',  # Косое захождение тип 2 = RA - AD2
        'DSA',  # Дневная полудуга 90 + AD
        'NSA',  # Ночная полудуга, 90 - AD
        'UMD',  # Отклонение от MC, Upper MD, |RAMC - RA|
        'LMD',  # Отклонение от IC, Lower MD, |RAIC - RA|
        'MD',  # Отклонение (расстояние) от меридиана = min(UMD,LMD)
        'quadrant',  # Квадрант 1 - ночний восток, 2 ночной запад и т.д.
        'SA',  # Длина текущей (ночной или дневной) полудуги
        'O_AD',  # Косое восхождение/зах. для восточн./западных планет
        'O_AD2',  # Косое восхождение/зах. тип 2 для восточн./западных планет
        # Для зодиакальных дирекций
        'RA_zod',  # Прямое восхождение точки проекции планеты на эклиптику
        'dec_zod',  # Зодиакальное склонение (точка проекции на эклиптику)
        'AD_zod',  # Зодиакальное разница восхождений
        'OA_zod',  # Зодиакальное косое восхождение = RA_zod - AD_zod
        'OD_zod',  # Зодиакальное косое захождение = RA_zod + AD_zod
        'AD2_zod',  # Разница вертексных восхождений точки проекции на эклиптику
        'OA2_zod',  # Зодиакальное косое восхождение тип 2 = RA_zod + AD2_zod
        'OD2_zod',  # Зодиакальное косое захождение тип 2 = RA_zod - AD2_zod
    ])

    # Астрономические данные - свойства небесной сферы
    Sphere = namedtuple('Sphere', [
        'RAMC',  # Прямое восхождение MC
        'RAIC',  # Прямое восхождение IC
        'OA_asc',  # Косое восхождение ASC
        'OD_dsc',  # Косое захождение DSC
    ])

    def __init__(self, data) -> None:

        self.geo_lat = data["lat"]
        self.jday = self.__get_jday(data)

        # Эклиптические координаты куспидов
        _cusps = swe.houses(
            self.jday,
            self.geo_lat,
            data["lon"],
            bytes("P", "ascii")
        )[0]

        # Константы для дальнейших преобразований
        self.earth_angle = self.__get_earth_angle()
        _ramc = self.__get_ra(ecl_lon=_cusps[9], ecl_lat=0.0)
        _raic = self.__get_ra(ecl_lon=_cusps[3], ecl_lat=0.0)
        self.sphere = self.Sphere(
            RAMC=_ramc,
            RAIC=_raic,
            OA_asc=self.__smart_add(_ramc, 90),
            OD_dsc=self.__smart_subst(_ramc, 90)
        )

        # Задаем экваториальные данные для планет
        self.planet = []
        for _p, _max_speed in enumerate([
            1.0197905235424334,  # Sun
            15.391811060863061,  # Moon
            2.202957566153169,  # Merc
            1.25890306853398,  # Venus
            0.7913933948763741,  # Mars
            0.24246019880425115,  # Jupiter
            0.13266330060580733,  # Saturn
            0.06470421363975136,
            0.04225855279702155,
            0.040509615350373694
        ]):
            _lon, _lat = swe.calc_ut(self.jday, _p)[0][:2]
            self.planet.append(
                self.get_equatorial_data(_lon, _lat, max_ecl_speed=_max_speed)
            )

        # Задаем экваториальные данные для домов
        self.house = []
        for _cusp in _cusps:
            self.house.append(
                self.get_equatorial_data(_cusp, 0.0)
            )

    @staticmethod
    def __smart_add(deg_a, deg_b):
        """
        Складывает два градуса на круге 0..360
        """
        if deg_a + deg_b > 360:
            return deg_a + deg_b - 360
        return deg_a + deg_b

    @staticmethod
    def __smart_subst(deg_a, deg_b):
        """
        Вычитает два градуса на круге 0..360
        """
        if deg_a - deg_b < 0:
            return deg_a - deg_b + 360
        return deg_a - deg_b

    @staticmethod
    def __get_jday(data: dict) -> float:
        """
        Возвращает юлианскую дату момента рождения
        """
        return swe.julday(
            data["year"], data["month"], data["day"],
            data["hour"] + data["min"]/60
        ) - data["tzone"] / 24

    def __get_earth_angle(self) -> float:
        """
        Возвращает угол наклона земной оси в момент рождения
        """
        _t = (self.jday - 2415020.0) / 36525
        _delta = (46.845 * _t + 0.0059 * (_t ** 2) -
                  0.001811 * (_t ** 3)) / 3600
        return 23.452294 - _delta

    def __get_ra(self,
                 ecl_lon: float,
                 ecl_lat: float) -> float:
        """
        Возвращает прямое восхождение по эклиптическим координатам
        """
        _ra = atan(
            (
                sin(ecl_lon/180*pi) * cos(self.earth_angle/180 * pi) -
                tan(ecl_lat/180*pi) * sin(self.earth_angle/180*pi)
            ) / cos(ecl_lon/180*pi)
        ) / pi * 180

        # Корректировка неопределенности функции arctan
        if ecl_lon == 90:
            _ra = 90
        elif ecl_lon == 270:
            _ra = 270
        elif 90 < ecl_lon < 270:
            _ra += 180
        elif 270 < ecl_lon < 360:
            _ra += 360

        return _ra

    def __get_decl(self,
                   ecl_lon: float,
                   ecl_lat: float) -> float:
        """
        Возвращает склонение по эклиптическим координатам
        """
        return asin(sin(ecl_lat/180*pi) * cos(self.earth_angle/180*pi) +
                    cos(ecl_lat/180*pi) * sin(self.earth_angle/180*pi) *
                    sin(ecl_lon/180*pi)) / pi * 180

    def __get_ad(self, dec: float) -> float:
        """
        Возвращает разницу восхождений объекта на экваторе и
        на текущей широте по заданному склонению объекта
        """
        return asin(tan(self.geo_lat/180*pi) * tan(dec/180*pi)) / pi * 180

    def __get_ad2(self, dec: float) -> float:
        """
        Возвращает разницу вертексных восхождений на экваторе и
        на текущей широте по заданному склонению объекта
        """
        _a = tan(dec/180*pi)  # max (tan 23.44)
        _b = tan(self.geo_lat/180*pi)
        _c = _a/_b

        # Убираем отклонения в далеких знаках
        # после запятой в вычислениях
        if _c < -1:
            _c = -1
        elif _c > 1:
            _c = 1
        return asin(_c) / pi * 180

    def __get_eclipt_lon(self,
                         ra: float,
                         dec: float) -> float:
        """
        Возвращает эклиптическую долготу по экваториальным координатам
        """
        _lon = atan(
            (
                sin(ra/180*pi) * cos(self.earth_angle/180 * pi) +
                tan(dec/180*pi) * sin(self.earth_angle/180*pi)
            ) / cos(ra/180*pi)
        ) / pi * 180

        # Корректировка неопределенности функции arctan
        if ra == 90:
            _lon = 90
        elif ra == 270:
            _lon = 270
        elif 90 < ra < 270:
            _lon += 180
        elif 270 < ra < 360:
            _lon += 360

        if _lon < 0:
            _lon += 360
        elif _lon > 360:
            _lon -= 360

        return _lon

    def __get_asc_or_vtx_from_ramc(self,
                                   ramc: float,
                                   vertex=False) -> float:
        """
        Возсращает эклиптическую долготу асцендента
        или вертекса по прямому восхождению MC
        """

        if vertex:
            _tan_phi = tan((90 - self.geo_lat)/180*pi)
        else:
            _tan_phi = tan(self.geo_lat/180*pi)
        _cos_ramc = cos(ramc/180*pi)
        _sin_ramc = sin(ramc/180*pi)
        _sin_e = sin(self.earth_angle/180*pi)
        _cos_e = cos(self.earth_angle/180*pi)
        _tan_e = tan(self.earth_angle/180*pi)

        _y = asin(_tan_phi * _tan_e)/pi*180
        if abs(ramc - (180 + _y)) < 1e-10:
            _asc_lon = 270
        elif abs(ramc - (360 - _y)) < 1e-10:
            _asc_lon = 90
        else:
            _x = atan(
                -_cos_ramc / (_tan_phi * _sin_e + _sin_ramc * _cos_e)
            ) / pi * 180

            if ramc < 180 + _y:
                _asc_lon = _x + 180
            elif 180 + _y < ramc < 270:
                _asc_lon = _x + 360
            elif 270 <= ramc < 360 - _y:
                _asc_lon = _x
            else:
                _asc_lon = _x + 180

            return _asc_lon

    def __get_vertex_lon(self, promissor: Coords) -> float:
        """
        Находит эклиптическую долготу вертекса/антивертекса
        (что ближе к промиссору) по экваториальным координатам
        промиссора
        """

        # Если соединяем с антивертексом
        # (промиссор на вотоке)
        if promissor.quadrant in [1, 4]:
            _ra = self.sphere.RAMC + 90 - promissor.AD2
        else:
            _ra = self.sphere.RAMC - 90 - promissor.AD2
        _dec = promissor.dec
        return self.__get_eclipt_lon(_ra, _dec)

    def get_equatorial_data(self,
                            lon: float,
                            lat: float,
                            max_ecl_speed: float = 0) -> Coords:
        """
        Вычисляет экваториальные данные по заданным
        эклиптическим координатам
        """

        # Данные для мунданных дирекций
        _ra = self.__get_ra(lon, lat)  # Right Ascention (RA)
        _dec = self.__get_decl(lon, lat)  # Declination (D)
        _ad = self.__get_ad(_dec)  # Ascention Diff (AD)
        _ad2 = self.__get_ad2(_dec)  # Ascention Diff Type 2 (AD2)

        # Данные для зодиакальных дирекций
        _ra_zod = self.__get_ra(lon, ecl_lat=0.0)  # RA of Zodiac degree
        _dec_zod = self.__get_decl(lon, ecl_lat=0.0)  # D of Zodiac degree
        _ad_zod = self.__get_ad(_dec_zod)  # AD of Zodiac degree
        _ad2_zod = self.__get_ad2(_dec_zod)  # AD22 of Zodiac degree

        # Формирование объекта Coords
        _umd = abs(self.sphere.RAMC - _ra)  # Upper Median Distance
        _lmd = abs(self.sphere.RAIC - _ra)  # Lower Median Distance
        if _umd > 180:
            _umd = 180 - _lmd
        if _lmd > 180:
            _lmd = 180 - _umd
        _dsa = self.__smart_add(90, _ad)  # Day Semi-Arch (DSA)
        _nsa = self.__smart_subst(90, _ad)  # Night Semi-Arch (NSA)
        _oa = self.__smart_subst(_ra, _ad)  # Oblique Ascention (OA)
        _od = self.__smart_add(_ra, _ad)  # Oblique Discention (OD)
        _oa2 = self.__smart_add(_ra, _ad2)  # OA type 2
        _od2 = self.__smart_subst(_ra, _ad2)  # OD type 2

        if self.sphere.RAMC < _ra < self.sphere.RAIC:
            _west = False
            _east = True
        else:
            _west = True
            _east = False

        # Если 1 или 7 дома, относим их к подгоризонтным
        if abs(_umd - _dsa) < 0.005 or abs(_lmd - _nsa) < 0.005:
            _above = False
            _below = True
        elif _umd < _dsa or _lmd > _nsa:
            _above = True
            _below = False
        else:
            _above = False
            _below = True

        if _east and _below:
            _quadrant = 1
        elif _west and _below:
            _quadrant = 2
        elif _west and _above:
            _quadrant = 3
        else:
            _quadrant = 4

        return self.Coords(
            lon=lon,
            lat=lat,
            max_ecl_speed=max_ecl_speed,
            dec=_dec,
            RA=_ra,
            AD=_ad,
            OA=_oa,
            OD=_od,
            AD2=_ad2,
            OA2=_oa2,
            OD2=_od2,
            DSA=_dsa,
            NSA=_nsa,
            UMD=_umd,
            LMD=_lmd,
            MD=min(_umd, _lmd),
            quadrant=_quadrant,
            SA=_dsa if _above else _nsa,
            O_AD=_oa if _east else _od,
            O_AD2=_oa2 if _east else _od2,
            RA_zod=_ra_zod,
            dec_zod=_dec_zod,
            AD_zod=_ad_zod,
            OA_zod=self.__smart_subst(_ra_zod, _ad_zod),
            OD_zod=self.__smart_add(_ra_zod, _ad_zod),
            AD2_zod=_ad2_zod,
            OA2_zod=self.__smart_add(_ra_zod, _ad2_zod),
            OD2_zod=self.__smart_subst(_ra_zod, _ad2_zod),
        )

    def print_coordinates(self) -> None:
        """
        Выводит на экран таблицу экваториальных данных планет
        """
        planets = ['SUN', 'MON', 'MER', 'VEN',
                   'MAR', 'JUP', 'SAT', 'URN', 'NEP', 'PLT']

        header = f"\n \t{'LONG':^6}\t{'LAT':^5}\t{'DEC':^7}\t{'RA':^6}\t{'MD':^6}"
        print(f"{header}\n{'-'*45}")

        for _p, name in enumerate(planets):
            row = f"{name}\t{self.planet[_p].lon:6.2f}\t{self.planet[_p].lat:5.2f}\t"
            row += f"{self.planet[_p].dec:6.2f}\t{self.planet[_p].RA:6.2f}\t"
            row += f"{self.planet[_p].MD:5.2f}"
            print(row)

        header = f"\n \t{'SA':^6}\t{'AD':^6}\t{'OA/OD':^6}\t{'AD2':^7}\t{'OA/OD2':^6}\n{'-'*46}"
        print(header)

        for _p, name in enumerate(planets):
            row = f"{name}\t{self.planet[_p].SA:6.2f}\t"
            row += f"{self.planet[_p].AD:6.2f}\t{self.planet[_p].O_AD:6.2f}\t"
            row += f"{self.planet[_p].AD2:6.2f}\t{self.planet[_p].O_AD2:6.2f}"
            print(row)
        print('\n')

    def aspect_mundi(self,
                     promissor: Coords,
                     significator: Coords,
                     aspect: int = 0,
                     option: int = 0,
                     ) -> float:
        """
        Рассчитывает мунданное соединение/аспект (дирекцию) в
        домах Плацида между двумя объектами с заданными эклип-
        тическими параметрами. Параметр option принимает значения:
        1 - рассчитывать соединение по параллели
        2 - рассчитывать соединение по контрпараллели
        3 - рассчитывать соединение с примарной вертикалью
        По-умолчанию рассчитывает соединение между объектами,
        но если задан параметр aspect (30, 60, 90...), то рассчи-
        тывает мунданный аспект
        """
        # Если ищем соединение с вертексом/контрвертексом
        if option == 3:
            if promissor.AD2 < 0:
                _ratio = (90 + promissor.AD2)/90
            else:
                _ratio = (90 - promissor.AD2)/90

            if promissor.quadrant in [4, 1]:
                if promissor.dec > 0:
                    _quadrant = 4
                else:
                    _quadrant = 1
            else:
                if promissor.dec > 0:
                    _quadrant = 3
                else:
                    _quadrant = 2

        # Если задан мунданный аспект (кроме соединения)
        elif aspect:
            _ratio = significator.MD / significator.SA
            _quadrant = significator.quadrant

            # Определяем Placidus Mundane Position (PMP)
            # сигнификатора
            if _quadrant == 1:
                _pmp = 90 * (1 - _ratio)
            elif _quadrant == 2:
                _pmp = 90 * (1 + _ratio)
            elif _quadrant == 3:
                _pmp = 90 * (3 - _ratio)
            else:
                _pmp = 90 * (3 + _ratio)

            # Определяем Placidus Mundane Position (PMP)
            # точки аспекта
            _pmp = self.__smart_subst(_pmp, abs(aspect))

            # Уточняем пропорции точки аспекта к длине полудуги
            # и положение точки аспекта в квадранте
            if _pmp < 90:
                _ratio = 1 - _pmp / 90
                _quadrant = 1
            elif 90 <= _pmp < 180:
                _ratio = _pmp / 90 - 1
                _quadrant = 2
            elif 180 <= _pmp < 270:
                _ratio = 3 - _pmp / 90
                _quadrant = 3
            else:
                _ratio = _pmp / 90 - 3
                _quadrant = 4

        # Ищем мунданные соединения
        else:
            _ratio = significator.MD/significator.SA
            _quadrant = significator.quadrant

            # Устраняем неопределенность дневных/ночных полудуг
            # для объектов близких к горизонту
            if (abs(significator.lon - self.house[0].lon) < 0.01 or
                    abs(significator.lon - self.house[6].lon) < 0.01):
                _ratio = 1

        _t = 1 if _quadrant in [1, 3] else -1
        if _quadrant in [1, 2]:
            _v, _r, _contr_r = -1, self.sphere.RAIC, self.sphere.RAMC
        else:
            _v, _r, _contr_r = 1, self.sphere.RAMC, self.sphere.RAIC

        # Соединение с мунданной параллелью сигнификатора
        if option == 1:
            _t = -_t

        # Соединение с мунданной контрпараллелью сигнификатора
        if option == 2:
            _t = -_t
            _v = -_v
            _r = _contr_r

        # print("{} - {} + ({}) * (90 + ({})*{}) * {}".format(
        #     promissor.RA, _r, _t, _v, promissor.AD, _ratio
        # ))
        return promissor.RA - _r + _t * (90 + _v * promissor.AD) * _ratio

    def aspect_zodiaco(self,
                       promissor: Coords,
                       significator: Coords,
                       aspect: int = 0,
                       option: int = 0) -> float:
        """
        Возвращает зодиакальный аспект (по заданным
        градусам - 30, 60, 90 и т.д.). Options
        принимает те же значения, что и в методе
        aspect_mundi()
        """
        if option == 3:
            _new_lon = self.__get_vertex_lon(promissor)
        else:
            _new_lon = self.__smart_add(promissor.lon, abs(aspect))
        _new_promissor = self.get_equatorial_data(lon=_new_lon, lat=0.0)
        return self.aspect_mundi(_new_promissor,
                                 significator,
                                 option=option)

    def aspect_field_plane(self,
                           promissor: Coords,
                           significator: Coords,
                           planetary_orbit: int,
                           aspect: int = 0,
                           option: int = 0):
        """
        Возвращает орбитальный аспект (по заданным
        градусам - 30, 60, +90 и т.д.) по орбите
        заданной планеты (0..9 - Солнце-Плутон).
        Options принимает те же значения, что и в
        методе aspect_mundi()
        """

        if option == 3:
            _target_lon = self.__get_vertex_lon(promissor)
            # return _target_lon
        else:
            _target_lon = self.__smart_add(promissor.lon, abs(aspect))

        # Находим момент времени, когда промиссор
        # окажется в новом положении
        if promissor.max_ecl_speed:
            _step = 1 / promissor.max_ecl_speed
        else:
            _step = 1

        _delta = self.__travel_time_to_ecl_lat(planetary_orbit,
                                               _target_lon,
                                               _step)

        _target_lat = swe.calc_ut(self.jday + _delta, planetary_orbit)[0][1]
        _new_promissor = self.get_equatorial_data(
            lon=_target_lon, lat=_target_lat
        )

        return self.aspect_mundi(
            promissor=_new_promissor,
            significator=significator,
            option=option
        )

    def __travel_time_to_ecl_lat(self,
                                 planet: int,
                                 target_lon: float,
                                 step: float) -> float:
        """
        Вычисляет кол-во юлианских дней, за которые планета
        дойдет до заданной точки эклиптики. Начальный поиск
        происходит с шагом _step (кол-во дней, за который
        планета движется примерно на 1 градус эклиптической
        дуги.
        """
        _day = 0
        while True:
            _lon1 = swe.calc_ut(self.jday + _day, planet)[0][0]
            _lon2 = swe.calc_ut(self.jday + _day + step,
                                planet)[0][0]
            if _lon1 < target_lon < _lon2 and abs(_lon2 - _lon1) < 30:
                break
            _day += step

        # Метод дихотомии
        _a = _day
        _b = _day + step
        while True:
            _c = _a + (_b - _a)/2
            _lon_a = swe.calc_ut(
                self.jday + _a, planet)[0][0] - target_lon
            _lon_c = swe.calc_ut(
                self.jday + _c, planet)[0][0] - target_lon
            if _lon_a * _lon_c > 0:
                _a = _c
            else:
                _b = _c

            if abs(_lon_c) < 0.1:
                break
        return _c

    def years_to_date(self, years: float) -> str:
        """
        Переводит год жизни в дату с момента рождения
        """
        date = swe.jdut1_to_utc(self.jday + years * 365.24, 1)
        return datetime(
            date[0], date[1], date[2], date[3], date[4], int(date[5])
        ).strftime("%Y, %B")

    @staticmethod
    def get_years_ptolemey(arc: float) -> float:
        return arc

    @staticmethod
    def get_years_naibod(arc: float) -> float:
        """
        Переводит дугу в годы жизни по методу Найбода
        """
        return arc * 1.0146

    def get_years_simmonite(self, arc: float) -> float:
        """
        Переводит дугу в годы жизни по методу Симмонита
        """
        _next_day_sun = self.get_equatorial_data(
            lon=(swe.calc_ut(self.jday + 1, 0)[0][0]),
            lat=0.0
        )

        _correction = abs(self.planet[0].RA - _next_day_sun.RA)

        return arc / _correction

    def get_years_placidus(self, arc: float) -> float:
        """
        Переводит дугу в годы жизни по методу Плацида
        """

        _ra = self.planet[0].RA + arc
        _lon = self.__get_eclipt_lon(_ra, dec=0.0)
        if _lon < 0:
            _lon += 360
        _time_to_travel = self.__travel_time_to_ecl_lat(planet=1,
                                                        target_lon=_lon,
                                                        step=1)
        return _time_to_travel

    def get_years_ascendant_arc(self, arc: float) -> float:
        """
        Переводит дугу в годы жизни по методу восходящей дуги
        """
        _progressed_ramc = self.sphere.RAMC + arc
        _progressed_asc = self.__get_asc_or_vtx_from_ramc(_progressed_ramc)
        _target_lon = _progressed_asc - self.house[0].lon + self.planet[0].lon
        if _target_lon > 360:
            _target_lon -= 360
        elif _target_lon < 0:
            _target_lon += 360

        _time_to_travel = self.__travel_time_to_ecl_lat(
            planet=1,
            target_lon=_target_lon,
            step=1
        )

        return _time_to_travel

    def get_years_vertical_arc(self, arc: float) -> float:
        """
        Переводит дугу в годы жизни по методу вертикальной дуги
        """
        _progressed_ramc = self.sphere.RAIC + arc
        _progressed_vtx = self.__get_asc_or_vtx_from_ramc(
            _progressed_ramc,
            vertex=True
        )
        _natal_vtx = self.__get_asc_or_vtx_from_ramc(
            self.sphere.RAIC,
            vertex=True
        )

        _target_lon = _progressed_vtx - _natal_vtx + self.planet[0].lon
        if _target_lon > 360:
            _target_lon -= 360
        elif _target_lon < 0:
            _target_lon += 360

        _time_to_travel = self.__travel_time_to_ecl_lat(
            planet=1,
            target_lon=_target_lon,
            step=1
        )

        return _time_to_travel
