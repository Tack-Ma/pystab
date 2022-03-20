import numpy
import sympy


class Circle:
    def __init__(self, do, di):
        self.do = do
        self.di = di

    @property
    def area(self):
        return sympy.pi / 4 * (self.do ** 2 - self.di ** 2)

    @property
    def area_o(self):
        return sympy.pi / 4 * self.do ** 2

    @property
    def area_i(self):
        return sympy.pi / 4 * self.di ** 2

    @property
    def I(self):
        return sympy.pi / 64 * (self.do ** 4 - self.di ** 4)


class RoundBar:
    def __init__(self, do, di, l, rho, E=21000):
        self.section = Circle(do, di)
        self.l = l
        self.material = Material(rho, E)

    @property
    def volume(self):
        return self.section.area * self.l

    @property
    def volume_o(self):
        return self.section.area_o * self.l

    @property
    def volume_i(self):
        return self.section.area_i * self.l

    @property
    def mass(self):
        return self.volume * self.material.rho

    @property
    def mass_o(self):
        return self.volume_o * self.material.rho

    @property
    def mass_i(self):
        return self.volume_i * self.material.rho

    @property
    def Ip(self):
        return (self.section.do ** 2 + self.section.di ** 2) * self.mass / 8

    @property
    def Id(self):
        return (self.mass_o - self.mass_i) * self.l ** 2 / 3


class Bearing(RoundBar):

    def __init__(self, do, di, k, c):
        super().__init__(do, di, 0, 0)
        self.k = k
        self.c = c

    @property
    def kxx(self):
        return self.k[0]

    @property
    def kxy(self):
        return self.k[1]

    @property
    def kyx(self):
        return self.k[2]

    @property
    def kyy(self):
        return self.k[3]


class Material:
    def __init__(self, rho, E):
        self.rho = rho
        self.E = E


Steel = Material(7.85*10**-6, 21000)


class LateralVibration:

    l, m, E, I, Id, Ip, omg, lmd = sympy.symbols(r'l m E I I_d I_p \omega \lambda')
    kxx, kxy, kyx, kyy, cxx, cxy, cyx, cyy = sympy.symbols(r'k_{xx} k_{xy} k_{yx} k_{yy} c_{xx} c_{xy} c_{yx} c_{yy}')

    stiffness_matrix = sympy.Matrix([
        [1, l, l ** 2 / (2 * E * I), l ** 3 / (6 * E * I), 0, 0, 0, 0],
        [0, 1, l / (E * I), l ** 2 / (2 * E * I), 0, 0, 0, 0],
        [0, 0, 1, l, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, l, l ** 2 / (2 * E * I), l ** 3 / (6 * E * I)],
        [0, 0, 0, 0, 0, 1, l / (E * I), l ** 2 / (2 * E * I)],
        [0, 0, 0, 0, 0, 0, 1, l],
        [0, 0, 0, 0, 0, 0, 0, 1]
    ])

    # 質点マトリクス
    mass_matrix = sympy.Matrix([
        [1,              0,         0, 0, 0,          0, 0, 0],
        [0,              1,         0, 0, 0,          0, 0, 0],
        [0,              Id*lmd**2, 1, 0, 0, omg*Ip*lmd, 0, 0],
        [-m*lmd**2,      0,         0, 1, 0,          0, 0, 0],
        [0,              0,         0, 0, 1,          0, 0, 0],
        [0,              0,         0, 0, 0,          1, 0, 0],
        [0,           -omg*Ip*lmd, 0,         0, 0,   Id*lmd**2, 1, 0],
        [0,           0, 0, 0, -m*lmd**2, 0,         0, 1]
    ])

    # 軸受マトリクス
    bearing_matrix = sympy.Matrix([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [-(kxx+lmd*cxx), 0, 0, 1, kxy+lmd*cxy, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [kyx+lmd*cyx, 0, 0, 0, -(kyy+lmd*cyy), 0, 0, 1]
    ])

# TODO: 初期化方法直したい
#     def __init__(self, do, di, l, rho, E=21000)
    def __init__(self, model=None):
        self.shape = model
        # self.shape = RoundBar(do, di, l, rho, E)

    # TODO: シャフトの場合と軸受の場合の変換関数が必要
    @staticmethod
    def f_matrix_pattern(param):
        if type(param) is Bearing:
            pass
        if type(param) is RoundBar:
            return LateralVibration.stiffness_matrix.subs([
                    (LateralVibration.l, param.l),
                    (LateralVibration.E, param.material.E),
                    (LateralVibration.I, param.section.I)
                ])

    @property
    def f_matrix(self):
        if self.shape is None:
            return None
        return [self.f_matrix_pattern(i) for i in self.shape]

    @property
    def mass_list(self):
        if self.shape is None:
            return None
        return [[i.mass, i.Id, i.Ip] for i in self.shape]

    @property
    def mass_array(self):
        if self.shape is None:
            return None
        return numpy.array(self.mass_list)

    @property
    def p_matrix(self):
        if self.shape is None:
            return None
        return [LateralVibration.stiffness_matrix.subs([
            (LateralVibration.l, i.shape.l),
            (LateralVibration.E, i.shape.material.E),
            (LateralVibration.I, i.shape.section.I)
        ]) for i in self.shape]