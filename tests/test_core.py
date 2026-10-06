"""geohash.core 的验收测试。

断言编码结果、格子几何与邻域不变量，覆盖正常编解码、边界归属、解码误差、
格子大小、八个方向的邻域、跨 180 度回绕、极点与非法输入。
"""

import unittest

from geohash import (
    GeoError,
    adjacent,
    bbox,
    cell_size,
    contains,
    decode,
    decode_exact,
    encode,
    neighbors,
)

DIRECTION_STEPS = {
    "n": (1, 0),
    "s": (-1, 0),
    "e": (0, 1),
    "w": (0, -1),
    "ne": (1, 1),
    "nw": (1, -1),
    "se": (-1, 1),
    "sw": (-1, -1),
}


class EncodingTests(unittest.TestCase):

    def test_known_cells_round_trip(self):
        self.assertEqual(encode(57.64911, 10.40744, 11), "u4pruydqqvj")
        self.assertEqual(encode(-90.0, -180.0, 1), "0")
        self.assertEqual(decode("ezs42"), (42.60498046875, -5.60302734375))
        self.assertEqual(bbox("ezs42"), (42.5830078125, -5.625,
                                         42.626953125, -5.5810546875))
        lat, lon = decode("ezs42")
        self.assertEqual(encode(lat, lon, 5), "ezs42")
        self.assertTrue(contains("ezs42", lat, lon))
        self.assertEqual(len(encode(0.0, 0.0, 7)), 7)


class BoundaryTests(unittest.TestCase):

    def test_boundary_points_belong_to_the_north_east_cell(self):
        self.assertEqual(encode(0.0, 0.0, 1), "s")
        self.assertEqual(encode(0.0, 0.0, 3), "s00")
        self.assertEqual(encode(0.0, 0.0, 5), "s0000")
        self.assertEqual(encode(45.0, 135.0, 1), "z")
        self.assertEqual(encode(45.0, 135.0, 5), "z0000")
        self.assertEqual(encode(-45.0, -90.0, 1), "6")
        lat_min, lon_min, _, _ = bbox(encode(0.0, 0.0, 3))
        self.assertEqual((lat_min, lon_min), (0.0, 0.0))


class CellGeometryTests(unittest.TestCase):

    def test_decode_error_radius_is_half_the_cell(self):
        self.assertEqual(decode_exact("ezs42")[2:],
                         (0.02197265625, 0.02197265625))
        for precision in (1, 2, 3, 4, 5):
            text = encode(10.0, 20.0, precision)
            lat, lon, lat_err, lon_err = decode_exact(text)
            lat_min, lon_min, lat_max, lon_max = bbox(text)
            self.assertAlmostEqual(lat, (lat_min + lat_max) / 2.0, places=9)
            self.assertAlmostEqual(lon, (lon_min + lon_max) / 2.0, places=9)
            self.assertAlmostEqual(lat_max - lat_min, 2.0 * lat_err, places=9)
            self.assertAlmostEqual(cell_size(precision)[0], 2.0 * lat_err,
                                   places=9)
            self.assertGreater(lon_err, 0.0)

    def test_longitude_error_matches_half_the_cell_width(self):
        for precision in (1, 3, 5):
            _, _, _, lon_err = decode_exact(encode(10.0, 20.0, precision))
            self.assertAlmostEqual(cell_size(precision)[1] / 2.0, lon_err,
                                   places=9)

    def test_cell_size_shrinks_with_precision(self):
        self.assertEqual(cell_size(1), (45.0, 45.0))
        self.assertEqual(cell_size(2), (5.625, 11.25))
        self.assertEqual(cell_size(3), (1.40625, 1.40625))
        heights = [cell_size(precision)[0] for precision in range(1, 7)]
        widths = [cell_size(precision)[1] for precision in range(1, 7)]
        self.assertEqual(heights, sorted(heights, reverse=True))
        self.assertEqual(widths, sorted(widths, reverse=True))
        self.assertEqual(len(set(heights)), len(heights))


class NeighborTests(unittest.TestCase):

    def test_adjacent_cells_share_an_edge(self):
        text = encode(10.0, 20.0, 4)
        lat_min, lon_min, lat_max, lon_max = bbox(text)
        height = lat_max - lat_min
        width = lon_max - lon_min
        around = neighbors(text)
        self.assertEqual(sorted(around), sorted(DIRECTION_STEPS))
        self.assertEqual(len(set(around.values())), 8)
        for name, (dlat, dlon) in DIRECTION_STEPS.items():
            self.assertEqual(around[name], adjacent(text, name))
            cell_lat_min, cell_lon_min, cell_lat_max, cell_lon_max = \
                bbox(around[name])
            self.assertAlmostEqual(cell_lat_min, lat_min + dlat * height,
                                   places=9)
            self.assertAlmostEqual(cell_lon_min, lon_min + dlon * width,
                                   places=9)
            self.assertAlmostEqual(cell_lat_max - cell_lat_min, height,
                                   places=9)
            self.assertAlmostEqual(cell_lon_max - cell_lon_min, width,
                                   places=9)
        self.assertEqual(neighbors("s"), {"n": "u", "s": "k", "e": "t",
                                          "w": "e", "ne": "v", "nw": "g",
                                          "se": "m", "sw": "7"})

    def test_antimeridian_wraps_both_ways(self):
        west = encode(10.0, -170.0, 2)
        east = encode(10.0, 170.0, 2)
        self.assertEqual((west, east), ("81", "xc"))
        self.assertEqual((bbox(west)[1], bbox(east)[3]), (-180.0, 180.0))
        self.assertEqual(adjacent(west, "w"), east)
        self.assertEqual(adjacent(east, "e"), west)
        self.assertEqual(neighbors(west)["w"], east)
        self.assertEqual(neighbors(east)["e"], west)
        self.assertEqual(adjacent(adjacent(west, "w"), "e"), west)
        self.assertEqual(bbox(adjacent(east, "e")), bbox(west))

    def test_pole_directions_have_no_neighbour(self):
        north = encode(89.0, 20.0, 2)
        south = encode(-89.0, 20.0, 2)
        self.assertEqual(bbox(north)[2], 90.0)
        self.assertEqual(bbox(south)[0], -90.0)
        self.assertIsNone(adjacent(north, "n"))
        self.assertIsNone(neighbors(north)["n"])
        self.assertIsNone(neighbors(north)["ne"])
        self.assertIsNone(neighbors(north)["nw"])
        self.assertIsNone(adjacent(south, "s"))
        self.assertIsNone(neighbors(south)["s"])
        self.assertIsNotNone(adjacent(north, "s"))
        self.assertIsNotNone(adjacent(north, "e"))
        self.assertEqual(bbox(adjacent(north, "s"))[2], bbox(north)[0])


class InputValidationTests(unittest.TestCase):

    def test_invalid_input_is_rejected(self):
        self.assertRaises(GeoError, decode, "")
        self.assertRaises(GeoError, decode, "u4pruydqqvjzz")
        self.assertRaises(GeoError, decode, "ezs42a")
        self.assertRaises(GeoError, decode, "ezs4o2")
        self.assertRaises(GeoError, decode, "EZs42")
        self.assertRaises(GeoError, bbox, "ezs4i2")
        self.assertRaises(GeoError, adjacent, "ezs42", "np")
        self.assertRaises(GeoError, encode, 90.5, 0.0, 4)
        self.assertRaises(GeoError, encode, 0.0, -180.5, 4)
        self.assertRaises(GeoError, encode, 0.0, 0.0, 0)
        self.assertRaises(GeoError, encode, 0.0, 0.0, 13)
        self.assertRaises(TypeError, encode, "10.0", 0.0, 4)
        self.assertRaises(TypeError, decode, 42)
        self.assertRaises(TypeError, contains, "ezs42", None, 0.0)


class InvariantTests(unittest.TestCase):

    def test_every_point_lies_in_its_own_cell(self):
        points = [(0.0, 0.0), (0.0, 90.0), (45.0, 135.0), (45.0, -90.0),
                  (-45.0, -90.0), (10.0, 20.0), (-12.345, 67.89),
                  (89.9, 179.9)]
        for lat, lon in points:
            for precision in (1, 2, 3):
                text = encode(lat, lon, precision)
                lat_min, lon_min, lat_max, lon_max = bbox(text)
                self.assertLessEqual(lat_min, lat)
                self.assertLessEqual(lat, lat_max)
                self.assertLessEqual(lon_min, lon)
                self.assertLessEqual(lon, lon_max)
                self.assertTrue(contains(text, lat, lon))


if __name__ == "__main__":
    unittest.main()
