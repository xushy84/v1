# 代码检查报告

检查时间: 自动生成

## 执行的检查
- `ruff check ...`（Python 静态检查）
- `python -m compileall -q .`（Python 语法编译检查）
- `bash -n <script>`（Shell 语法检查）

## 按文件结果
- `README.md`: ⚠️ 未做语义校验（文档文件）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/README.md`: ⚠️ 未做语义校验（文档文件）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/10_prepare_relax0.py`: ❌ 有问题（7）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/11_submit_relax0.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/12_collect_relax0.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/20_prepare_mag_relax.py`: ❌ 有问题（4）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/21_submit_mag_relax.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/22_collect_mag_gs.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/40_prepare_static_gap.py`: ❌ 有问题（4）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/41_submit_static_gap.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/42_collect_gap.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/50_symmetry_analyze.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/51_prepare_polar_kick.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/52_submit_polar_kick.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/53_update_polar_flag.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/60_prepare_relax_tight.py`: ❌ 有问题（3）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/61_submit_relax_tight.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/62_collect_relax_tight.py`: ✅ 通过（0）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/70_prepare_phonon.py`: ❌ 有问题（6）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/71_submit_phonon.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/72_analyze_phonon.py`: ❌ 有问题（8）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/80_prepare_berryP.py`: ❌ 有问题（1）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/81_submit_berryP.sh`: ✅ 通过（bash -n）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/82_collect_berryP.py`: ❌ 有问题（2）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/10_prepare_relax0.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/12_collect_relax0.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/20_prepare_mag_relax.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/22_collect_mag_gs.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/40_prepare_static_gap.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/42_collect_gap.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/50_symmetry_analyze.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/51_prepare_polar_kick.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/53_update_polar_flag.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/60_prepare_relax_tight.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/62_collect_relax_tight.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/70_prepare_phonon.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/72_analyze_phonon.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/80_prepare_berryP.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/82_collect_berryP.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/_common.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/__pycache__/_vasp_inputs.cpython-311.pyc`: ⚠️ 未检查（非源码或二进制）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/_common.py`: ❌ 有问题（3）
- `vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/_vasp_inputs.py`: ✅ 通过（0）

## Python问题汇总（ruff）
### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/10_prepare_relax0.py
- L13 `F401`: `_common.log` imported but unused
- L13 `F401`: `_common.vlog` imported but unused
- L71 `E702`: Multiple statements on one line (semicolon)
- L112 `E701`: Multiple statements on one line (colon)
- L113 `E701`: Multiple statements on one line (colon)
- L114 `E701`: Multiple statements on one line (colon)
- L161 `E701`: Multiple statements on one line (colon)

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/12_collect_relax0.py
- L17 `F401`: `_common.log` imported but unused
- L17 `F401`: `_common.vlog` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/20_prepare_mag_relax.py
- L9 `F401`: `_common.log` imported but unused
- L9 `F401`: `_common.vlog` imported but unused
- L31 `E702`: Multiple statements on one line (semicolon)
- L103 `E701`: Multiple statements on one line (colon)

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/22_collect_mag_gs.py
- L18 `F401`: `_common.log` imported but unused
- L18 `F401`: `_common.vlog` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/40_prepare_static_gap.py
- L6 `F401`: `_common.log` imported but unused
- L6 `F401`: `_common.vlog` imported but unused
- L68 `E702`: Multiple statements on one line (semicolon)
- L79 `E701`: Multiple statements on one line (colon)

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/42_collect_gap.py
- L25 `F401`: `_common.log` imported but unused
- L25 `F401`: `_common.vlog` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/50_symmetry_analyze.py
- L22 `F401`: `_common.log` imported but unused
- L22 `F401`: `_common.vlog` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/51_prepare_polar_kick.py
- L24 `F401`: `_common.log` imported but unused
- L24 `F401`: `_common.vlog` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/53_update_polar_flag.py
- L19 `F401`: `_common.log` imported but unused
- L19 `F401`: `_common.vlog` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/60_prepare_relax_tight.py
- L46 `E701`: Multiple statements on one line (colon)
- L48 `E701`: Multiple statements on one line (colon)
- L170 `E701`: Multiple statements on one line (colon)

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/70_prepare_phonon.py
- L50 `F401`: `_common.vlog` imported but unused
- L99 `E702`: Multiple statements on one line (semicolon)
- L114 `E701`: Multiple statements on one line (colon)
- L119 `E702`: Multiple statements on one line (semicolon)
- L122 `E701`: Multiple statements on one line (colon)
- L196 `F401`: `numpy` imported but unused

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/72_analyze_phonon.py
- L24 `F811`: Redefinition of unused `os` from line 22: `os` redefined here
- L26 `F401`: `typing.Tuple` imported but unused
- L28 `F401`: `_common.log` imported but unused
- L28 `F401`: `_common.vlog` imported but unused
- L75 `F401`: `phonopy` imported but unused
- L76 `F401`: `phonopy.Phonopy` imported but unused
- L78 `F401`: `phonopy.file_IO.parse_FORCE_SETS` imported but unused
- L131 `F841`: Local variable `pa` is assigned to but never used

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/80_prepare_berryP.py
- L127 `E701`: Multiple statements on one line (colon)

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/82_collect_berryP.py
- L23 `F401`: `_common.vlog` imported but unused
- L63 `E702`: Multiple statements on one line (semicolon)

### vasp_multiferroic_project_full_params_FINAL/vasp_multiferroic_project_full_params/scripts/_common.py
- L7 `E401`: Multiple imports on one line
- L81 `E402`: Module level import not at top of file
- L82 `E402`: Module level import not at top of file
