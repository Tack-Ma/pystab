import numpy
import sympy


class Circle:
    """
    円断面の材料力学上の属性を求めるクラス

    Attributes
    ----------
    do : float
        円の外径
    di : float
        円の内径
    """
    def __init__(self, do, di):
        """

        Parameters
        ----------
        do : float
            円の外径
        di : float
            円の内径
        """
        self.do = do
        self.di = di

    @property
    def area(self):
        """

        Returns
        -------
        area : float
            断面積
        """
        return sympy.pi / 4 * (self.do ** 2 - self.di ** 2)

    @property
    def area_o(self):
        """

        Returns
        -------
        area_of_outer_circle : float
            外径円の面積
        """
        return sympy.pi / 4 * self.do ** 2

    @property
    def area_i(self):
        """

        Returns
        -------
        area_of_internal_circle : float
            内径円の面積

        """
        return sympy.pi / 4 * self.di ** 2

    @property
    def I(self):
        """

        Returns
        -------
        area_moment_of_inertia : float
            断面2次モーメント

        """
        return sympy.pi / 64 * (self.do ** 4 - self.di ** 4)


class RoundBar:
    def __init__(self, do, di, length, rho, E=21000.0):
        self.section = Circle(do, di)
        self.length = length
        self.material = Material(rho, E)

    @property
    def volume(self):
        return self.section.area * self.length

    @property
    def volume_o(self):
        return self.section.area_o * self.length

    @property
    def volume_i(self):
        return self.section.area_i * self.length

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
        return (self.mass_o - self.mass_i) * self.length ** 2 / 3


class Shaft(RoundBar):

    def __init__(self, do, di, length, rho, E=21000.0, add_weight=0):
        super().__init__(do, di, length, rho, E)
        self.add_weight = add_weight

    @property
    def mass(self):
        return self.volume * self.material.rho + self.add_weight


class Bearing(RoundBar):
    """
    軸受の定義

    """
    def __init__(self, do, di, k: list, c: list):
        """

        Parameters
        ----------
        do : float
            Journal outer diameter
        di : float
            Journal inter diameter
        k : list
            spring parameter list (xx, xy, yx, yy)
        c : list
            damping parameter list (xx, xy, yx, yy)
        """
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

    @property
    def cxx(self):
        return self.c[0]

    @property
    def cxy(self):
        return self.c[1]

    @property
    def cyx(self):
        return self.c[2]

    @property
    def cyy(self):
        return self.c[3]


class Material:
    """
    材料ライブラリクラス

    Attributes
    ----------
    rho : float
        材料の密度
    E : float
        材料の縦弾性係数

    """
    def __init__(self, rho: float, E: float):
        """

        Parameters
        ----------
        rho : float
            材料の密度
        E : float
            材料の縦弾性係数
        """
        self.rho = rho
        self.E = E


class MaterialLibrary(Material):
    """
    Materialクラスの材料データをJSONファイルから取得する


    Attributes
    ----------
    """

    def __init__(self, material_name, library_database):
        """
        Parameters
        ----------
        material_name : str
            材料名称
        library_database : str
            材料ライブラリのパス
        """
        # super().__init__()
        pass


# 材料のライブラリ
# TODO: JSONデータ化する
Steel = Material(7.85 * 10 ** -6, 21000)


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
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, Id * lmd ** 2, 1, 0, 0, omg * Ip * lmd, 0, 0],
        [-m * lmd ** 2, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, -omg * Ip * lmd, 0, 0, 0, Id * lmd ** 2, 1, 0],
        [0, 0, 0, 0, -m * lmd ** 2, 0, 0, 1]
    ])

    # 軸受マトリクス
    bearing_matrix = sympy.Matrix([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [-(kxx + lmd * cxx), 0, 0, 1, kxy + lmd * cxy, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [kyx + lmd * cyx, 0, 0, 0, -(kyy + lmd * cyy), 0, 0, 1]
    ])

    # TODO: 初期化方法直したい
    #     def __init__(self, do, di, l, rho, E=21000)
    def __init__(self):
        self.shape = ()
        # self.shape = RoundBar(do, di, l, rho, E)

    @staticmethod
    def f_matrix_pattern(param):
        """
        シャフトの場合と軸受の場合の変換関数

        Parameters
        ----------
        param :
            Shaft or Bearing class
        Returns
        -------
        Matrix : sympy
            f matrix

        """
        if type(param) is Bearing:
            return LateralVibration.bearing_matrix.subs([
                (LateralVibration.kxx, param.kxx),
                (LateralVibration.kxy, param.kxy),
                (LateralVibration.kyx, param.kyx),
                (LateralVibration.kyy, param.kyy),
                (LateralVibration.cxx, param.cxx),
                (LateralVibration.cxy, param.cxy),
                (LateralVibration.cyx, param.cyx),
                (LateralVibration.cyy, param.cyy)
            ])
        if type(param) is Shaft:
            return LateralVibration.stiffness_matrix.subs([
                (LateralVibration.l, param.l),
                (LateralVibration.E, param.material.E),
                (LateralVibration.I, param.section.I)
            ])

    @staticmethod
    def p_matrix_param(param):
        return LateralVibration.mass_matrix.subs([
            (LateralVibration.m, param[0]),
            (LateralVibration.Ip, param[1]),
            (LateralVibration.Id, param[2])
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

    # 質点マトリクスの定義
    @property
    def p_matrix(self):
        if self.shape is None:
            return None
        tmp_array = numpy.vstack(
            numpy.zeros(len(self.mass_array[0])),
            self.mass_array,
            numpy.zeros(len(self.mass_array[0]))
        )
        new_tmp_array = tmp_array[:-1]/2+tmp_array[1:]/2
        return [self.p_matrix_param(i) for i in new_tmp_array.tolist()]

    @property
    def t_matrix(self):
        return [f*p for f, p in zip(self.f_matrix, self.p_matrix[1:])]

    # TODO: 自由振動の固有振動数を求める
    def solve_free_vibration(self):
        pass

