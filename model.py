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


class Material:
    def __init__(self, rho, E):
        self.rho = rho
        self.E = E


Steel = Material(7.85*10**-6, 21000)


class TMatrix:
    def __init__(self, do, di, l, rho, E=21000):
        self.shape = RoundBar(do, di, l, rho, E)

    @property
    def f_matrix(self):
        pass

    @property
    def p_in(self):
        pass
