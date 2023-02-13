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

using Coords = std::array<double, 4>;
using CoordsList = std::list<Coords>;
const float epsilon = std::numeric_limits<float>::epsilon();

float distance( const Coords& coords1, 
                const Coords& coords2 ) 
{
  // Euclidean square root distance
  return sqrt( (coords1[0] - coords2[0]) * (coords1[0] - coords2[0]) +
               (coords1[1] - coords2[1]) * (coords1[1] - coords2[1]) );
}

std::vector<size_t> getNeighboringCells(const long int cellCol,
                                        const long int cellRow,
                                        const long int xSize,
                                        const long int ySize,
                                        const long int radius) 
{
    
  std::vector<size_t> neighboringCells;

  // window of 2 * radius + 1 centered on the current cell coordinates
  long int minCol = std::max<long int>(0, cellCol - radius);
  long int minRow = std::max<long int>(0, cellRow - radius);
  long int maxCol = std::min<long int>(xSize-1, cellCol + radius);
  long int maxRow = std::min<long int>(ySize-1, cellRow + radius);

  for(long int r = minRow; r <= maxRow; r++){
    for(long int c = minCol; c <= maxCol; c++){
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

float getGaussian(const std::vector< CoordsList >& gridToInterpol,
                  const std::vector<size_t>& neighbors,
                  const size_t cellCol,
                  const size_t cellRow,
                  const size_t xSize,
                  const size_t ySize,
                  const float sigma,
		  const long int radius,
		  const float resolution)
{

  // We take the coordinates at the center of the pixel
  Coords cellCoordsCenter = Coords{ static_cast<float>(cellCol + 0.5f),
                                    static_cast<float>(cellRow + 0.5f),
                                    0.f, 0.f };
  
  float sumNumerator = 0.f;
  float sumDenominator = 0.f;
  bool atLeastOne = false;
  float weight;
  float dist;

  for(const auto neigh : neighbors){
    for(const auto& coords: gridToInterpol[neigh]){
      dist = distance(cellCoordsCenter, coords);
      if (dist < 0.5 + radius) {
	dist *= resolution;
	weight = exp( (- dist * dist) / (2 * sigma * sigma) );
	sumNumerator += (coords[2] * coords[3])  * weight;
	sumDenominator += coords[3] * weight;
	atLeastOne = true;
      }
    }
  }

  if(atLeastOne){
    return sumNumerator / sumDenominator;
  } else {
    return 0.f;
  }

}

std::vector<double> pointCloudToDSM(const std::vector<double>& pos,
				    const size_t nbBands, 
				    const size_t nbPoints,
				    const float xStart, 
				    const float yStart, 
				    const size_t xSize,
				    const size_t ySize,
				    const float resolution,
				    const size_t radius,
				    const float sigma,
				    const bool considerConfidence,
				    const bool trace)
{
  size_t N = pos.size();
  const size_t outSize = xSize*ySize;
  std::vector<double> output(outSize);

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
  std::vector< CoordsList > gridToInterpol(outSize);
  double x, y, z, confidence, col, row;
  size_t cellCol, cellRow, rowByCol;

  for ( size_t k = 0 ; k < nbPoints ; ++k ) {
    
    x = pos[(1*nbPoints)+k];
    y = pos[(2*nbPoints)+k];
    z = pos[(3*nbPoints)+k];

    if(considerConfidence){
      confidence = pos[(12*nbPoints)+k];
    } else {
      confidence = 1;
    }

    col = (x - xStart) / resolution;
    row = (yStart - y) / resolution;

    cellCol = floor(col);
    cellRow = floor(row);

    if ((cellCol >= 0) && (cellCol < xSize) && (cellRow >= 0) && (cellRow < ySize)) {
      rowByCol= cellCol + cellRow * xSize;
      gridToInterpol[rowByCol].push_back(Coords{col, row, z, confidence});
    }

  }

  const size_t idwRadius = radius;

  // Loop over the grid to interpolate the z for each cell
  for ( size_t k = 0 ; k < outSize ; ++k ) {

    cellCol = k % xSize;
    cellRow = k / xSize;

    // Get neighboring cells with radius defined by the user
    auto neighbors = getNeighboringCells(cellCol,
                                         cellRow,
                                         xSize,
                                         ySize,
                                         idwRadius);
    
    
    
    // Gaussian interpolation
    output[k] = getGaussian(gridToInterpol,
                            neighbors,
                            cellCol,
                            cellRow,
                            xSize,
                            ySize,
                            sigma,
			    idwRadius,
			    resolution);


    // Get idw value
    // output[k] = getIdw(gridToInterpol,
    //                    neighbors,
    //                    cellCol,
    //                    cellRow,
    //                    xSize,
    //                    ySize);
  }

  return output;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

// wrap C++ function with NumPy array IO
py::array pyPointCloudToDSM(py::array_t<double, py::array::c_style | py::array::forcecast> array,
			    float xStart,
			    float yStart,
			    size_t xSize,
			    size_t ySize,
			    float resolution,
			    size_t radius,
			    float sigma,
			    bool considerConfidence,
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
  // there is not conflicting type, however the function pointCloudToDSM needs a double * ptr as
  // input instead of a std::vector.
  std::memcpy(pos.data(),array.data(),array.size()*sizeof(double));

  size_t nbBands = array.shape()[0]-2;
  size_t nbPoints = array.shape()[1];

  // call pure C++ function
  std::vector<double> result = pointCloudToDSM(pos, nbBands, nbPoints,
					       xStart, yStart,
					       xSize, ySize,
					       resolution,
					       radius,
                                               sigma,
                                               considerConfidence,
					       trace);

  ssize_t             ndim    = 2;
  std::vector<size_t> shape   = { xSize, ySize };
  std::vector<size_t> strides = { sizeof(double)*ySize, sizeof(double) };

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

  m.def("pc_to_dsm", &pyPointCloudToDSM, "Convert point cloud to digital surface model");
}
