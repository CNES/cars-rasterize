#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <array>
#include <list>
#include <cmath>
#include <algorithm>
#include <limits>

// -------------
// pure C++ code
// -------------

using Coords = std::array<float, 4>;
using CoordsList = std::list<Coords>;
const float epsilon = std::numeric_limits<float>::epsilon();

float distance( const Coords& coords1, 
                const Coords& coords2 ) 
{
  // Euclidean square root distance
  return sqrt( (coords1[0] - coords2[0]) * (coords1[0] - coords2[0]) +
               (coords1[1] - coords2[1]) * (coords1[1] - coords2[1]) );
}

std::vector<size_t> getNeighboringCells(const size_t cellCol,
                                        const size_t cellRow,
                                        const size_t xSize,
                                        const size_t ySize,
                                        const size_t radius) 
{
    
  std::vector<size_t> neighboringCells;

  // window of 2 * radius + 1 centered on the current cell coordinates
  size_t minCol = std::max<size_t>(0, cellCol - radius);
  size_t minRow = std::max<size_t>(0, cellRow - radius);
  size_t maxCol = std::min<size_t>(xSize-1, cellCol + radius);
  size_t maxRow = std::min<size_t>(ySize-1, cellRow + radius);

  for(size_t r = minRow; r <= maxRow; r++){
    for(size_t c = minCol; c <= maxCol; c++){
      neighboringCells.push_back( r * xSize + c );
    }
  }

  return neighboringCells;
}

float getIdw(const std::vector< CoordsList >& gridToInterpol,
             const std::vector<size_t>& neighbors,
             const size_t cellCol,
             const size_t cellRow,
             const size_t xSize,
             const size_t ySize) 
{
  
  // We take the coordinates at the center of the pixel
  Coords cellCoordsCenter = Coords{ static_cast<float>(cellCol + 0.5f),
                                    static_cast<float>(cellRow + 0.5f),
                                    0.f, 0.f };

  float sumNumerator = 0.f;
  float sumDenominator = 0.f;
  float dist;
  bool atLeastOne = false;
  for(const auto neigh : neighbors){
    for(const auto& coords: gridToInterpol[neigh]){
      dist = distance(cellCoordsCenter, coords) + epsilon;
      sumNumerator += (coords[2] * coords[3]) / dist;
      sumDenominator += coords[3] / dist;
      atLeastOne = true;
    }
  }

  if(atLeastOne){
    return sumNumerator / sumDenominator;
  } else {
    return 0.f;
  }

}

std::vector<double> pc_to_dsm(const std::vector<double>& pos,
			      const size_t nb_bands, 
                              const size_t nb_points,
			      const float xstart, 
                              const float ystart, 
                              const size_t xsize,
			      const size_t ysize,
			      const float resolution,
			      const size_t radius,
			      const bool trace)
{
  size_t N = pos.size();

  if (trace) {
    std::cout << "vector size: " << N << std::endl;
    std::cout << "nb bands (excepted (x,y)): " << nb_bands << std::endl;
    std::cout << "nb points: " << nb_points << std::endl;
    std::cout << "dstwin (xoff yoff xsize ysize): (" 
	      << xstart << " " << ystart << " "
	      << xsize << " " << ysize << ")" << std::endl;
    std::cout << "resolution: " << resolution << std::endl;
  }
  const size_t outsize = xsize*ysize;
  std::vector<double> output(outsize);

  size_t cellDCol = xsize / 2;
  size_t cellDRow = ysize / 2;

  // All bands
  // data_valid  0            
  // x           1  
  // y           2
  // z           3 
  // msk         4
  // clr0        5
  // clr1        6
  // clr2        7
  // clr3        8
  // coord_epi_geom_i   9  
  // coord_epi_geom_j   10 
  // idx_im_epi         11
  // ambiguity_confidence 12

  // First attempt idw
  // Second attempt idw weighted with ambiguity
  // For each target pixel, we store a list of fractionnal Z coordinates
  std::vector< CoordsList > gridToInterpol(outsize);
  float x, y, z, confidence, col, row;
  size_t cellCol, cellRow, rowByCol;

  for ( size_t k = 0 ; k < nb_points ; ++k ) {
    
    x = pos[(1*nb_points)+k];
    y = pos[(2*nb_points)+k];
    z = pos[(3*nb_points)+k];
    confidence = pos[(12*nb_points)+k];

    col = (x - xstart) / resolution;
    row = (ystart - y) / resolution;

    cellCol = std::lround(col);
    cellRow = std::lround(row);

    if ((cellCol >= 0) && (cellCol < xsize) && (cellRow >= 0) && (cellRow < ysize)) {
      rowByCol= cellCol + cellRow * xsize;
      gridToInterpol[rowByCol].push_back(Coords{col, row, z, confidence});
    }

  }

  const size_t idwRadius = radius;

  // Loop over the grid to interpolate the z for each cell
  for ( size_t k = 0 ; k < outsize ; ++k ) {

    cellCol = k % xsize;
    cellRow = k / xsize;

    // Get neighboring cells with radius defined by the user
    auto neighbors = getNeighboringCells(cellCol,
                                         cellRow,
                                         xsize,
                                         ysize,
                                         idwRadius);
    
    // Get idw value
    output[k] = getIdw(gridToInterpol,
                       neighbors,
                       cellCol,
                       cellRow,
                       xsize,
                       ysize);
  }

  // for ( size_t k = 0 ; k < nb_points ; ++k ) {
  //   const float x = pos[(1*nb_points)+k];
  //   const float y = pos[(2*nb_points)+k];
  //   const float z = pos[(3*nb_points)+k];

  //   const ssize_t i = std::lround((x - xstart) / resolution);
  //   const ssize_t j = std::lround((ystart - y) / resolution);

  //   if (trace) {
  //     std::cout << "pt " << k << ": "
  // 	<< "x, y, z > " 
  // 	<< x << ", " << y << ", " << z << std::endl;
  //     std::cout << "i, j > " << i << ", " << j << std::endl;
  //   }

  //   if ((i >= 0) && (i < xsize) && (j >= 0) && (j < ysize)) {
  //     const size_t ij = i+j*xsize;
  //     output[ij] = z;
  //   }
  // }

  return output;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

// wrap C++ function with NumPy array IO
py::array py_pc_to_dsm(py::array_t<double, py::array::c_style | py::array::forcecast> array,
		       float xstart,
		       float ystart,
		       size_t xsize,
		       size_t ysize,
		       float resolution,
		       size_t radius,
		       bool trace)
{
  // check input dimensions
  if ( array.ndim()     != 2 )
    throw std::runtime_error("Input should be 2-D NumPy array");

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> pos(array.size());

  // copy py::array -> std::vector
  // Do we need to make a copy here ? Is there something similar
  // to cython memory view
  // https://stackoverflow.com/questions/54793539/pybind11-modify-numpy-array-from-c
  // It seems that passing the py::array by reference should do the job if
  // there is not conflicting type, however the function pc_to_dsm needs a double * ptr as
  // input instead of a std::vector.
  std::memcpy(pos.data(),array.data(),array.size()*sizeof(double));

  size_t nb_bands = array.shape()[0]-2;
  size_t nb_points = array.shape()[1];

  // call pure C++ function
  std::vector<double> result = pc_to_dsm(pos, nb_bands, nb_points,
					 xstart, ystart, xsize, ysize,
					 resolution,
					 radius,
					 trace);

  ssize_t             ndim    = 2;
  std::vector<size_t> shape   = { xsize, ysize };
  std::vector<size_t> strides = { sizeof(double)*ysize, sizeof(double) };

  // return 2-D NumPy array, I think here it is ok since the expected argument is
  // a pointer so there is no copy
  return py::array(py::buffer_info(
				   result.data(),                           /* data as contiguous array  */
				   sizeof(double),                          /* size of one scalar        */
				   py::format_descriptor<double>::format(), /* data type                 */
				   ndim,                                    /* number of dimensions      */
				   shape,                                   /* shape of the matrix       */
				   strides                                  /* strides for each axis     */
				   ));
}

// wrap as Python module
PYBIND11_MODULE(rasterize, m)
{
  m.doc() = "rasterize";

  m.def("pc_to_dsm", &py_pc_to_dsm, "Convert point cloud to digital surface model");
}
