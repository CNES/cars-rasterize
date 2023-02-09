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


def generate_synthetic_point_cloud(
    x_offset: float, y_offset: float, nb_rows: int, nb_cols: int, reso: float
):
    """
    Convention bord pixel pour x_offset et y_offset
    """

    x = np.arange(nb_cols) * reso + x_offset + reso / 2
    y = np.arange(nb_rows) * -1 * reso - y_offset - reso / 2
    z = np.arange(nb_rows * nb_cols) * 100
    xgrid, ygrid = np.meshgrid(x, y)
    pointcloud = np.vstack((np.ravel(xgrid), np.ravel(ygrid), z))

    # // All bands
    #   // data_valid  0
    #   // x           1
    #   // y           2
    #   // z           3
    #   // msk         4
    #   // clr0        5
    #   // clr1        6
    #   // clr2        7
    #   // clr3        8
    #   // coord_epi_geom_i   9
    #   // coord_epi_geom_j   10
    #   // idx_im_epi         11
    #   // ambiguity_confidence 12

    pointcloud = np.vstack(
        (
            np.zeros(nb_rows * nb_cols),
            pointcloud,
            *np.zeros((8, nb_rows * nb_cols)),
            np.ones(nb_rows * nb_cols, dtype=np.float32),
        )
    )

    return pointcloud


if __name__ == "__main__":
    XSTART = 0.0
    YSTART = 0.0
    XSIZE = 2
    YSIZE = 2
    RESOLUTION = 1.0
    RADIUS = 3
    SIGMA = 0.5

    POINTCLOUD = generate_synthetic_point_cloud(
        XSTART, YSTART, XSIZE, YSIZE, RESOLUTION
    )

    # pylint: disable=c-extension-no-member
    dsm = rasterize.pc_to_dsm(
        POINTCLOUD,
        XSTART,
        YSTART,
        XSIZE,
        YSIZE,
        RESOLUTION,
        RADIUS,
        SIGMA,
        True,
    )

    transform = Affine.translation(XSTART, YSTART)
    transform = transform * Affine.scale(RESOLUTION, -RESOLUTION)

    profile = DefaultGTiffProfile(
        count=1, dtype=dsm.dtype, width=XSIZE, height=YSIZE, transform=transform
    )

    with rio.open("dsm.tif", "w", **profile) as dst:
        dst.write(dsm, 1)
