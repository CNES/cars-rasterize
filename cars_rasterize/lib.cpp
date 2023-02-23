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


float standard_deviation( std::vector<float> vect)
{
  double sum = std::accumulate(vect.begin(), vect.end(), 0.0);
  double mean = sum / vect.size();
  std::vector<double> diff(vect.size());
  std::transform(vect.begin(), vect.end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  return std::sqrt(sq_sum / vect.size());
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


std::pair<float, float> getGaussian(const std::vector< CoordsList >& gridToInterpol,
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
                                    0.f, 0.f};
  
  float sumNumerator = 0.f;
  float sumDenominator = 0.f;
  float weight;
  float dist;
  
  // std::vector<float> heights;
  std::list < std::pair<float, Coords> > coordsList;

  size_t nbPoints = 0; 

  for(const auto neigh : neighbors){
    for(const auto& coords: gridToInterpol[neigh]){
      dist = distance(cellCoordsCenter, coords);
      if (dist < 0.5 + radius) {
	weight = exp( (- dist * dist) / (2 * sigma * sigma) );
	sumNumerator += coords[2] * weight;
	sumDenominator += weight;
	nbPoints++;
      }
    }
  }

  if(nbPoints > 0){
    return std::make_pair(sumNumerator / sumDenominator,
    			  nbPoints);

  } else {
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(), 0.f);
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
				    const float sigma)
{
  const size_t outSize = xSize*ySize;
  std::vector<double> output(2*outSize);
  std::vector<bool> mask(nbPoints);

  // For each target pixel, we store a list of fractionnal Z coordinates
  std::vector< CoordsList > gridToInterpol(outSize);
  double x, y, z, col, row;
  size_t cellCol, cellRow, rowByCol;

  for ( size_t k = 0 ; k < nbPoints ; ++k ) {
    
    x = pos[k];
    y = pos[(1*nbPoints)+k];
    z = pos[(2*nbPoints)+k];

    col = (x - xStart) / resolution;
    row = (yStart - y) / resolution;

    cellCol = floor(col);
    cellRow = floor(row);

    if ((cellCol >= 0) && (cellCol < xSize) && (cellRow >= 0) && (cellRow < ySize)) {
      rowByCol= cellCol + cellRow * xSize;
      gridToInterpol[rowByCol].push_back(Coords{col, row, z, static_cast<float> (k)});
    }

  }

  // Loop over the grid to interpolate the z for each cell
  for ( size_t k = 0 ; k < outSize ; ++k ) {

    cellCol = k % xSize;
    cellRow = k / xSize;

    // Get neighboring cells with radius defined by the user
    auto neighbors = getNeighboringCells(cellCol,
                                         cellRow,
                                         xSize,
                                         ySize,
                                         radius);
    
    // Gaussian interpolation
    auto [mean, stdev] = getGaussian(gridToInterpol,
				     neighbors,
				     cellCol,
				     cellRow,
				     xSize,
				     ySize,
				     sigma,
				     radius,
				     resolution);
				     
    output[2*k+0] = mean;
    output[2*k+1] = stdev;
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
			    float sigma)
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
                                               sigma);

  ssize_t             ndim    = 3;
  std::vector<size_t> shape   = { xSize, ySize, 2 };
  std::vector<size_t> strides = { 2*ySize*sizeof(double), 2*sizeof(double), sizeof(double)};

  // return 2-D NumPy array, I think here it is ok since the expected argument is
  // a pointer so there is no copy
  return py::array(py::buffer_info(result.data(),                           /* data as contiguous array  */
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
