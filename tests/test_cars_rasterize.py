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
"""Tests for `cars-rasterize` package."""

from typing import Any, Dict

import numpy as np

# Third party imports
import pytest

# rasterize imports
import rasterize


@pytest.fixture
def synthetic_cloud() -> Dict[str, Any]:
    """Generate synthetic point cloud."""
    cloud = {
        "xstart": 0.0,
        "ystart": 0.0,
        "xsize": 2,
        "ysize": 2,
        "resolution": 1.0,
    }  # type: Dict[str, Any]

    x = np.arange(cloud["xsize"]) * cloud["resolution"]
    x += cloud["xstart"] + cloud["resolution"] / 2
    y = np.arange(cloud["ysize"]) * -cloud["resolution"]
    y -= cloud["ystart"] + cloud["resolution"] / 2
    xgrid, ygrid = np.meshgrid(x, y)

    cloud["points"] = np.vstack((xgrid.ravel(), ygrid.ravel()))
    cloud["values"] = np.arange(cloud["xsize"] * cloud["ysize"]) * 100
    cloud["values"] = cloud["values"][np.newaxis, ...]
    cloud["valid"] = np.ones((1, cloud["points"].shape[1]))
    return cloud


@pytest.mark.parametrize("radius, sigma", [[0, 1], [0, 2], [1, np.inf]])
def test_synthetic(
    synthetic_cloud, radius, sigma
):  # pylint: disable=redefined-outer-name
    """Test on synthetic points cloud."""
    out, mean, stdev, nb_pts_in_disc, nb_pts_in_cell = rasterize.pc_to_dsm(
        synthetic_cloud["points"],
        synthetic_cloud["values"],
        synthetic_cloud["valid"],
        synthetic_cloud["xstart"],
        synthetic_cloud["ystart"],
        synthetic_cloud["xsize"],
        synthetic_cloud["ysize"],
        synthetic_cloud["resolution"],
        radius,
        sigma,
    )

    if radius == 0:
        mean_ref = np.array([[[0], [100]], [[200], [300]]])
        out_ref = mean_ref
        stdev_ref = np.zeros((2, 2, 1))
        nb_pts_in_disc_ref = np.ones((2, 2))
        nb_pts_in_cell_ref = np.ones((2, 2))

    elif radius == 1 and sigma == np.inf:
        mean_ref = 150 * np.ones((2, 2, 1))
        out_ref = mean_ref
        stdev_ref = np.float32(50 * np.sqrt(5)) * np.ones((2, 2, 1))
        nb_pts_in_disc_ref = np.array([[4, 4], [4, 4]])
        nb_pts_in_cell_ref = np.array([[1, 1], [1, 1]])

    assert np.array_equal(out, out_ref)
    assert np.array_equal(mean, mean_ref)
    assert np.array_equal(stdev, stdev_ref)
    assert np.array_equal(nb_pts_in_disc, nb_pts_in_disc_ref)
    assert np.array_equal(nb_pts_in_cell, nb_pts_in_cell_ref)
