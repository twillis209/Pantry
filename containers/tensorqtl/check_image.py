"""Check the tensorQTL image: baseline versions, torch<->numpy interop, R qvalue.

Run at build time (Dockerfile) and on the HPC; prints CUDA availability.
"""

import os
import sys

import numpy
import pandas
import pandas_plink  # noqa: F401
import scipy
import tensorqtl
import torch
from rpy2.robjects.packages import importr

got = {
    "python": sys.version.split()[0],
    "tensorqtl": tensorqtl.__version__,
    "torch": torch.__version__,
    "pandas": pandas.__version__,
    "numpy": numpy.__version__,
    "scipy": scipy.__version__,
}
want = {
    "python": "3.12.16",
    "tensorqtl": "1.0.10",
    "torch": "2.3.1+cu121",
    "pandas": "2.3.3",
    "numpy": "2.5.4",
    "scipy": "1.18.1",
}
print(got)
assert got == want, f"versions differ from the baseline: {got} != {want}"
# numpy/scipy must be the pip wheels (bundled OpenBLAS), as in the baseline.
for mod in (numpy, scipy):
    libs = os.path.join(os.path.dirname(os.path.dirname(mod.__file__)), mod.__name__ + ".libs")
    assert os.path.isdir(libs), f"{mod.__name__} is not the pip wheel ({libs} missing)"
assert torch.from_numpy(numpy.ones(3)).numpy().sum() == 3
importr("qvalue")
print("rpy2 + qvalue: ok")
print("cuda available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device:", torch.cuda.get_device_name(0))
