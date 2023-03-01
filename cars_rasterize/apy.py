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
        points = np.vstack((las.x, las.y))
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

    # save all images
    # out: gaussian interpolation
    transform = Affine.translation(roi["xstart"], roi["ystart"])
    transform = transform * Affine.scale(resolution, -resolution)

    profile = DefaultGTiffProfile(
        count=out.shape[-1],
        dtype=out.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=np.nan,
    )

    with rio.open(dsm_out, "w", **profile) as dst:
        for band in range(out.shape[-1]):
            dst.write(out[..., band], band + 1)

    # mean: simple mean
    profile = DefaultGTiffProfile(
        count=mean.shape[-1],
        dtype=mean.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=np.nan,
    )

    with rio.open(dsm_out + "_mean.tif", "w", **profile) as dst:
        for band in range(mean.shape[-1]):
            dst.write(mean[..., band], band + 1)

    # stdev: standard deviation
    profile = DefaultGTiffProfile(
        count=stdev.shape[-1],
        dtype=stdev.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=np.nan,
    )

    with rio.open(dsm_out + "_stdev.tif", "w", **profile) as dst:
        for band in range(stdev.shape[-1]):
            dst.write(stdev[..., band], band + 1)

    # nb_pts_in_disc: nb points used for interpolation
    profile = DefaultGTiffProfile(
        count=1,
        dtype=nb_pts_in_disc.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=0,
    )

    with rio.open(dsm_out + "_nb_disc.tif", "w", **profile) as dst:
        dst.write(nb_pts_in_disc, 1)

    # nb_pts_in_cell: nb points in interpolated cell
    profile = DefaultGTiffProfile(
        count=1,
        dtype=nb_pts_in_cell.dtype,
        width=roi["xsize"],
        height=roi["ysize"],
        transform=transform,
        nodata=0,
    )

    with rio.open(dsm_out + "_nb_cell.tif", "w", **profile) as dst:
        dst.write(nb_pts_in_cell, 1)
