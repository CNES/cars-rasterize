<div align="center">
  <a href="https://github.com/CNES/cars"><img src="https://raw.githubusercontent.com/CNES/cars-rasterize/master/docs/images/picto_transparent.png" alt="CARS" title="CARS"  width="20%"></a>

<h4>cars-rasterize</h4>

[![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)](CONTRIBUTING.md)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0/)

<p>
  <a href="#overview">Overview</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#how-it-works">How It Works</a> •
  <a href="#contribution">Contribution</a>
</p>
</div>

## Overview

**cars-rasterize** aims to convert a point cloud into a digital surface (or terrain) model with colors.

It is a part of the  photogrammetry tool [cars](https://github.com/cnes/cars) extracting Digital Surface Models from satellite images.

## Installation
**cars-rasterize** is available on Pypi and can be installed by:
```
pip install cars-rasterize
```

## Quick start

1. Download **subsampled_nimes.laz***:
```
wget https://raw.githubusercontent.com/CNES/cars-rasterize/master/data/subsampled_nimes.laz
```

subsampled_nimes.laz |
:-------------------------:|
<img src="https://raw.githubusercontent.com/CNES/cars-rasterize/master/docs/images/nimes.gif" alt="drawing" width="400"/> 

[subsampled_nimes.laz*](./data/subsampled_nimes.laz) is from https://geoservices.ign.fr/lidarhd. and has been downsampled (1 point every 50cm) to make the file smaller.

2. Run **las2tif** executable:
```
las2tif subsampled_nimes.laz dsm.tif --clr_out clr.tif
```

3. ✅ Done! The executable generates two files:
- **dsm.tif**: the elevation of the points (Z dimension) are projected into a regular grid to generate a raster file named Digital Surface Model.
- **clr.tif**: the red, the green and the blue dimensions can be also projected producing a color interpretation map superimposable on DSM

dsm.tif |  clr.tif
:-------------------------:|:-------------------------:
<img src="https://raw.githubusercontent.com/CNES/cars-rasterize/master/docs/images/nimes_elevation.png" alt="drawing" width="300"/>|   <img src="https://raw.githubusercontent.com/CNES/cars-rasterize/master/docs/images/nimes_colors.png" alt="drawing" width="300"/>


## How it works

A LAS file contains a set of points $P = \{(x, y, z, r, g, b)_k\}$ each having several dimensions:
- $x$ and $y$ correspond to planimetric information
- $z$ corresponds to the altitude
- $r$, $g$ and $b$ correspond to colorimetric information (respectively red, green, blue )


To create a raster digital surface model, we define a regular grid on a region of interest **roi** of origin $(x_{start}, y_{start})$, size $(x_{size}, y_{size})$ with a constant **resolution**.

For each cell of center $(c_x, c_y)$, we consider the subset of points contained in the disk $D$ (parameter **radius**) centered on this cell (see figure below):

Contributing points |
:-------------------------:|
<img src="https://raw.githubusercontent.com/CNES/cars-rasterize/master/docs/images/contributing_points.png" alt="drawing" width="600"/>

Then, the altitude assigned $z(c_x, c_y)$ to the cell is a Gaussian  weighted average (standard deviation **sigma** $\sigma$) of the distance $d$ to its center :

$$z(c_x, c_y) = \frac{\sum_{p_k \in D} z_k e^{-d_k^2/2\sigma^2}}{\sum_{p_k \in D} e^{-d_k^2/2\sigma^2}}$$

Finally, to have a superimposable color to this dsm, the colors are averaged in the same way.

## Contribution
**cars-rasterize** is a free software: Apache Software License 2.0. See [Contribution](./CONTRIBUTING.md) manual.