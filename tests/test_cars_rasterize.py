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
"""Tests for `cars_rasterize` package."""

# Third party imports
import pytest

# cars_rasterize imports
import cars_rasterize


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # Example to edit
    return "response"


def test_content(response):  # pylint: disable=redefined-outer-name
    """Sample pytest test function with the pytest fixture as an argument."""
    # Example to edit
    print(response)


def test_cars_rasterize():
    """Sample pytest cars_rasterize module test function"""
    assert cars_rasterize.__author__ == "CNES"
    assert cars_rasterize.__email__ == "cars@cnes.fr"
