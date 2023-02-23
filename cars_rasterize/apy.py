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

import json
from pathlib import Path

import laspy
import numpy as np
import rasterio as rio
from rasterio.profiles import DefaultGTiffProfile
from rasterio.transform import Affine

import rasterize


def main(cloud_in, dsm_out, resolution=0.5, radius=1, sigma=None, roi=None):
    """
    Convert point cloud to dsm
    """
    with laspy.open(cloud_in) as creader:
        las = creader.read()
        pointcloud = np.vstack((las.x, las.y, las.z))

    if roi is None:
        attrs = str(Path(cloud_in).with_suffix("")) + "_attrs.json"
        if Path(attrs).exists():
            with open(attrs, "r", encoding="utf-8") as attrs_reader:
                roi = json.load(attrs_reader)["attributes"]
        else:
            roi = {
                "xmin": resolution
                * ((np.amin(las.x) - resolution / 2) // resolution),
                "ymax": resolution
                * ((np.amax(las.y) + resolution / 2) // resolution),
                "xmax": resolution
                * ((np.amax(las.x) + resolution / 2) // resolution),
                "ymin": resolution
                * ((np.amin(las.y) - resolution / 2) // resolution),
            }

        roi["xstart"] = roi["xmin"]
        roi["ystart"] = roi["ymax"]
        roi["xsize"] = (roi["xmax"] - roi["xmin"]) / resolution
        roi["ysize"] = (roi["ymax"] - roi["ymin"]) / resolution

    if sigma is None:
        sigma = resolution

    # pylint: disable=c-extension-no-member
    dsm = rasterize.pc_to_dsm(
        pointcloud,
        roi["xstart"],
        roi["ystart"],
        int(roi["xsize"]),
        int(roi["ysize"]),
        resolution,
        radius,
        sigma,
    )

    transform = Affine.translation(roi["xstart"], roi["ystart"])
    transform = transform * Affine.scale(resolution, -resolution)

    profile = DefaultGTiffProfile(
        count=2,
        dtype=dsm.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=np.nan,
    )

    with rio.open(dsm_out, "w", **profile) as dst:
        dst.write(dsm[..., 0], 1)
        dst.write(dsm[..., 1], 2)
