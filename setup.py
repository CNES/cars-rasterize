#!/usr/bin/env python
# coding: utf8
#
# Copyright (C) 2023 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of cars-rasterize
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
Packaging setup.py for compatibility
All packaging in setup.cfg
"""

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

try:
    ext_modules = [
        Pybind11Extension("rasterize", ["cars_rasterize/lib.cpp"]),
    ]

    setup(
        ext_modules=ext_modules,
        # Currently, build_ext only provides an optional "highest supported C++
        # level" feature, but in the future it may provide more features.
        cmdclass={"build_ext": build_ext},
    )

except Exception:
    print(
        "\n\nAn error occurred while building the project, "
        "please ensure you have the most updated version of pip, setuptools, "
        "setuptools_scm and wheel with:\n"
        "   pip install -U pip setuptools setuptools_scm wheel\n\n"
    )
    raise
