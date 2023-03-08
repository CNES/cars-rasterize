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


def main(
    cloud_in,
    dsm_out,
    clr_out=None,
    resolution=0.5,
    radius=1,
    sigma=None,
    roi=None,
):
    """
    Convert point cloud to dsm
    """
    with laspy.open(cloud_in) as creader:
        las = creader.read()
        points = np.vstack((las.x, las.y))
        if clr_out is None:
            values = las.z
            values = values[np.newaxis, ...]
        else:
            values = np.vstack((las.z, las.red, las.green, las.blue))
        valid = np.ones((1, points.shape[1]))

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
    out, mean, stdev, nb_pts_in_disc, nb_pts_in_cell = rasterize.pc_to_dsm(
        points,
        values,
        valid,
        roi["xstart"],
        roi["ystart"],
        int(roi["xsize"]),
        int(roi["ysize"]),
        resolution,
        radius,
        sigma,
    )

    # reshape data as a 2d grid.
    shape_out = (int(roi["ysize"]), int(roi["xsize"]))
    out = out.reshape(shape_out + (-1,))
    mean = mean.reshape(shape_out + (-1,))
    stdev = stdev.reshape(shape_out + (-1,))
    nb_pts_in_disc = nb_pts_in_disc.reshape(shape_out)
    nb_pts_in_cell = nb_pts_in_cell.reshape(shape_out)

    # save dsm
    # out: gaussian interpolation
    transform = Affine.translation(roi["xstart"], roi["ystart"])
    transform = transform * Affine.scale(resolution, -resolution)

    profile = DefaultGTiffProfile(
        count=1,
        dtype=out.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=np.nan,
    )

    with rio.open(dsm_out, "w", **profile) as dst:
        dst.write(out[..., 0], 1)

    if clr_out is not None:
        # clr: color r, g, b
        transform = Affine.translation(roi["xstart"], roi["ystart"])
        transform = transform * Affine.scale(resolution, -resolution)

        profile = DefaultGTiffProfile(
            count=3,
            dtype=out.dtype,
            width=roi["xsize"],
            height=roi["ysize"],
            transform=transform,
            nodata=np.nan,
        )

        with rio.open(clr_out, "w", **profile) as dst:
            for band in range(3):
                dst.write(out[..., band + 1], band + 1)
