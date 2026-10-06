# geohash

只依赖 Python 标准库的地理编码内核：经纬度交错位编码与解码、精度对应的位长、
格子大小与解码误差、格子边界框、点是否落在格子里，以及八个方向的相邻格子。
所有输入由调用方给出，不使用真实时钟、网络或随机数，同一输入永远得到同样的结果。

- `geohash/core.py` — 编码内核：base32 文本、格子几何、邻域与回绕。
- `tests/test_core.py` — 验收用例。

## 运行测试

在项目根目录执行：

```
python3 -m unittest discover -s tests -v
```

Windows 上如果没有 `python3`，可用：

```
python -m unittest discover -s tests -v
```
