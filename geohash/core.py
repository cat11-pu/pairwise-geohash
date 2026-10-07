"""geohash —— 地理编码与邻域查询内核（纯标准库，行为完全确定）。

把经纬度按交错位编成 base32 文本，也能从文本还原格子中心、格子边界、
解码误差与八个方向的相邻格子。输入全部由调用方给出：模块不读写文件、
不联网、不使用真实时钟或随机数，同一输入永远得到同样的结果。

约定：

* 精度以字符计，取值 1..12；每个字符 5 位，其中经度先出位，
  所以总位数为奇数时经度比纬度多一位；
* 格子是闭合区间：落在格子边界上的点同时属于相邻两侧的格子，
  编码时判给北侧或东侧的那一格；
* 南北方向没有回绕，极点方向的邻居不存在；东西方向跨过 180 度时绕到对面。
"""

BASE32 = "0123456789bcdefghjkmnpqrstuvwxyz"

LAT_MIN = -90.0
LAT_MAX = 90.0
LON_MIN = -180.0
LON_MAX = 180.0

LAT_SPAN = LAT_MAX - LAT_MIN
LON_SPAN = LON_MAX - LON_MIN

MIN_PRECISION = 1
MAX_PRECISION = 12

__all__ = [
    "BASE32",
    "GeoError",
    "adjacent",
    "bbox",
    "bit_lengths",
    "cell_size",
    "contains",
    "decode",
    "decode_exact",
    "encode",
    "neighbors",
]

_INDEX = {char: index for index, char in enumerate(BASE32)}

_DIRECTIONS = {
    "n": (1, 0),
    "s": (-1, 0),
    "e": (0, 1),
    "w": (0, -1),
    "ne": (1, 1),
    "nw": (1, -1),
    "se": (-1, 1),
    "sw": (-1, -1),
}


class GeoError(ValueError):
    """坐标越界、精度非法、文本非法或方向未知时抛出。"""


def _check_coords(lat, lon):
    """确认纬度与经度是范围内的数值。"""
    for value, label, low, high in ((lat, "纬度", LAT_MIN, LAT_MAX),
                                    (lon, "经度", LON_MIN, LON_MAX)):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError("%s 必须是数值: %r" % (label, value))
        if not low <= value <= high:
            raise GeoError("%s 超出范围 [%r, %r]: %r" % (label, low, high, value))


def _check_precision(precision):
    """确认精度是 1..12 的整数。"""
    if isinstance(precision, bool) or not isinstance(precision, int):
        raise TypeError("精度必须是整数: %r" % (precision,))
    if not MIN_PRECISION <= precision <= MAX_PRECISION:
        raise GeoError("精度必须在 %d..%d 之间: %r"
                       % (MIN_PRECISION, MAX_PRECISION, precision))


def _check_hash(text):
    """确认文本是 1..12 个字符的字符串。"""
    if not isinstance(text, str):
        raise TypeError("geohash 文本必须是字符串: %r" % (text,))
    if not MIN_PRECISION <= len(text) <= MAX_PRECISION:
        raise GeoError("geohash 长度必须在 %d..%d 之间: %r"
                       % (MIN_PRECISION, MAX_PRECISION, text))
    for char in text:
        if char not in _INDEX:
            raise GeoError("geohash 含字母表外字符 %r: %r" % (char, text))
    return text


def _char_value(char):
    """取一个字符对应的五位值。"""
    try:
        return _INDEX[char]
    except KeyError:
        raise GeoError("geohash 含字母表外字符: %r" % (char,))


def bit_lengths(precision):
    """返回该精度下 (纬度位数, 经度位数)；经度先出位。"""
    _check_precision(precision)
    total = precision * 5
    return total // 2, total - total // 2


def cell_size(precision):
    """返回该精度下格子的大小 (纬度高度, 经度宽度)。"""
    _check_precision(precision)
    lat_bits, lon_bits = bit_lengths(precision)
    return LAT_SPAN / (1 << lat_bits), LON_SPAN / (1 << lon_bits)


def encode(lat, lon, precision=12):
    """把经纬度编成指定精度的 geohash 文本。"""
    _check_coords(lat, lon)
    _check_precision(precision)
    lat_min, lat_max = LAT_MIN, LAT_MAX
    lon_min, lon_max = LON_MIN, LON_MAX
    chars = []
    chunk = 0
    filled = 0
    take_lon = True
    while len(chars) < precision:
        if take_lon:
            mid = (lon_min + lon_max) / 2.0
            if lon >= mid:
                chunk = (chunk << 1) | 1
                lon_min = mid
            else:
                chunk = chunk << 1
                lon_max = mid
        else:
            mid = (lat_min + lat_max) / 2.0
            if lat >= mid:
                chunk = (chunk << 1) | 1
                lat_min = mid
            else:
                chunk = chunk << 1
                lat_max = mid
        take_lon = not take_lon
        filled += 1
        if filled == 5:
            chars.append(BASE32[chunk])
            chunk = 0
            filled = 0
    return "".join(chars)


def _bounds(text):
    """返回格子的 (纬度下界, 纬度上界, 经度下界, 经度上界)。"""
    lat_min, lat_max = LAT_MIN, LAT_MAX
    lon_min, lon_max = LON_MIN, LON_MAX
    take_lon = True
    for char in text:
        value = _char_value(char)
        for shift in range(4, -1, -1):
            bit = (value >> shift) & 1
            if take_lon:
                mid = (lon_min + lon_max) / 2.0
                if bit:
                    lon_min = mid
                else:
                    lon_max = mid
            else:
                mid = (lat_min + lat_max) / 2.0
                if bit:
                    lat_min = mid
                else:
                    lat_max = mid
            take_lon = not take_lon
    return lat_min, lat_max, lon_min, lon_max


def decode_exact(text):
    """解码：返回 (纬度中心, 经度中心, 纬度误差半径, 经度误差半径)。"""
    _check_hash(text)
    lat_min, lat_max, lon_min, lon_max = _bounds(text)
    lat = (lat_min + lat_max) / 2.0
    lon = (lon_min + lon_max) / 2.0
    return (lat, lon,
            (lat_max - lat_min) / 2.0,
            (lon_max - lon_min) / 2.0)


def decode(text):
    """解码：返回格子中心的 (纬度, 经度)。"""
    lat, lon, _, _ = decode_exact(text)
    return lat, lon


def bbox(text):
    """返回格子的外接框 (纬度下界, 经度下界, 纬度上界, 经度上界)。"""
    _check_hash(text)
    lat_min, lat_max, lon_min, lon_max = _bounds(text)
    return lat_min, lon_min, lat_max, lon_max


def contains(text, lat, lon):
    """判断经纬度是否落在该 geohash 的格子里；边界上的点算在格子里。"""
    _check_hash(text)
    _check_coords(lat, lon)
    lat_min, lat_max, lon_min, lon_max = _bounds(text)
    return (lat_min <= lat <= lat_max) and (lon_min <= lon <= lon_max)


def _to_values(text):
    """把文本展开成 (纬度位值, 经度位值, 纬度位数, 经度位数)。"""
    lat_bits, lon_bits = bit_lengths(len(text))
    lat_value = 0
    lon_value = 0
    index = 0
    for char in text:
        value = _char_value(char)
        for shift in range(4, -1, -1):
            bit = (value >> shift) & 1
            if index % 2 == 0:
                lon_value = (lon_value << 1) | bit
            else:
                lat_value = (lat_value << 1) | bit
            index += 1
    return lat_value, lon_value, lat_bits, lon_bits


def _to_hash(lat_value, lon_value, lat_bits, lon_bits):
    """把两个方向的位值按经度先出的顺序拼回 geohash 文本。"""
    bits = []
    for index in range(lat_bits + lon_bits):
        if index % 2 == 0:
            bits.append((lon_value >> (lon_bits - 1 - index // 2)) & 1)
        else:
            bits.append((lat_value >> (lat_bits - 1 - index // 2)) & 1)
    chars = []
    for start in range(0, len(bits), 5):
        chunk = 0
        for bit in bits[start:start + 5]:
            chunk = (chunk << 1) | bit
        chars.append(BASE32[chunk])
    return "".join(chars)


def adjacent(text, direction):
    """返回相邻格子的 geohash 文本。

    direction 取 n / s / e / w / ne / nw / se / sw 之一。南北方向没有邻居时
    返回 None；东西方向跨过 180 度时从对面绕回来。
    """
    _check_hash(text)
    if direction not in _DIRECTIONS:
        raise GeoError("未知方向: %r" % (direction,))
    lat_value, lon_value, lat_bits, lon_bits = _to_values(text)
    dlat, dlon = _DIRECTIONS[direction]
    lat_value += dlat
    if lat_value < 0 or lat_value >= (1 << lat_bits):
        return None
    lon_limit = 1 << lon_bits
    lon_value += dlon
    lon_value %= lon_limit
    return _to_hash(lat_value, lon_value, lat_bits, lon_bits)


def neighbors(text):
    """返回八个方向的相邻格子；该方向没有邻居时取 None。"""
    _check_hash(text)
    return {name: adjacent(text, name) for name in _DIRECTIONS}
