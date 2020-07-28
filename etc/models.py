from astropy import units as u

class Site(object):
    """Model for a site location and the atmosphere above it"""

    def __init__(self, name=None, altitude=None, latitude= None, longitude=None, **kwargs):
        self.name = name if name is not None else "Undefined"
        self.altitude = altitude * u.m
        self.latitude = latitude
        self.longitude = longitude

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

    def __str__(self):
        return "{}: lon={}, lat={}, altitude={}".format(self.name, self.longitude, self.latitude, self.altitude)


class Telescope(object):
    def __init__(self, name=None, size=0, area=0, num_mirrors=2, **kwargs):
        self.name = name if name is not None else "Undefined"
        self.size = size * u.m
        self.area_unit = u.m * u.m
        self.area = area * self.area_unit
        self.num_mirrors = num_mirrors

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

    def __str__(self):
        return "{} (M1: {} diameter, {} area; {} mirrors)".format(self.name, self.size.to(u.m), self.area.to(self.area_unit), self.num_mirrors)


class Instrument(object):
    def __init__(self, name=None, **kwargs):
        self.name = name if name is not None else "Undefined"

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

