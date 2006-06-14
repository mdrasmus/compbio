import util
import algorithms


class Region:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def overlap(self, region):
        return util.overlap(self.start, self.end, region.start, region.end)


class EndPoint:
    def __init__(self, region, boundary):
        self.region = region
        self.boundary = boundary


