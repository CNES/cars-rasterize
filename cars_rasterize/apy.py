#!/usr/bin/env python
# coding: utf8
#
# Copyright (C) 2023 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of cars_rasterize
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
Python api for cars_rasterize.
"""

import numpy as np
import rasterio as rio
from rasterio.profiles import DefaultGTiffProfile
from rasterio.transform import Affine

import rasterize


def main(cloud, infos):
    """
    Convert point cloud to dsm
    """
    layers = [cloud[key].to_numpy() for key in cloud]
    resolution = infos["resolution"]
    xstart = infos["xstart"]
    ystart = infos["ystart"]
    xsize = infos["xsize"]
    ysize = infos["ysize"]
    radius = infos["radius"]
    sigma = -1
    if infos["sigma"]:
        sigma = infos["sigma"]

    pointcloud = np.vstack(layers)

    # pylint: disable=c-extension-no-member
    dsm = rasterize.pc_to_dsm(
        pointcloud,
        xstart,
        ystart,
        xsize,
        ysize,
        resolution,
        int(radius),
        sigma,
        False,
    )

    transform = Affine.translation(xstart, ystart)
    transform = transform * Affine.scale(resolution, -resolution)

    profile = DefaultGTiffProfile(
        count=1, dtype=dsm.dtype, width=xsize, height=ysize, transform=transform
    )

    with rio.open("dsm.tif", "w", **profile) as dst:
        dst.write(dsm, 1)
