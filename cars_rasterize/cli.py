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
Console script for cars_rasterize.
"""

import argparse
import sys

from cars_rasterize import apy


def main():
    """Console script for cars_rasterize."""
    parser = argparse.ArgumentParser()
    parser.add_argument("cloud_in")
    parser.add_argument("dsm_out")
    parser.add_argument("--clr_out", default=None)
    parser.add_argument("--resolution", default=0.5, type=float)
    parser.add_argument("--radius", default=1, type=int)
    parser.add_argument("--sigma", default=None, type=float)
    parser.add_argument(
        "--roi",
        nargs=4,
        metavar=("xstart", "ystart", "xsize", "ysize"),
        default=None,
        type=float,
    )

    args = parser.parse_args()
    user_roi = None
    if args.roi is not None:
        user_roi = {
            "xstart": args.roi[0],
            "ystart": args.roi[1],
            "xsize": int(args.roi[2]),
            "ysize": int(args.roi[3]),
        }

    apy.main(
        args.cloud_in,
        args.dsm_out,
        clr_out=args.clr_out,
        resolution=args.resolution,
        radius=args.radius,
        sigma=args.sigma,
        roi=user_roi,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
