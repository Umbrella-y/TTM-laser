class Grid:
    # Grid类，里面保存了密度，热导率，热容，长，宽，电子温度，晶格温度信息
    def __init__(self, density, conductivity, heat_capacity, length, width,temperature, lattice_temp):
        self.density = density
        self.conductivity = conductivity
        self.heat_capacity = heat_capacity
        self.length = length
        self.width = width
        self.temperature = temperature
        self.lattice_temp = lattice_temp

