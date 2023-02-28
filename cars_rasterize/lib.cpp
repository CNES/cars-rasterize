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

using Coords = std::array<double, 3>;
using CoordsList = std::list<Coords>;
const float epsilon = std::numeric_limits<float>::epsilon();

float distance( const Coords& coords1, 
                const Coords& coords2 ) 
{
  // Euclidean square root distance
  return sqrt( (coords1[0] - coords2[0]) * (coords1[0] - coords2[0]) +
               (coords1[1] - coords2[1]) * (coords1[1] - coords2[1]) );
}


std::pair<float, float> vector_statistics(std::vector<float> vect)
{
  double sum = std::accumulate(vect.begin(), vect.end(), 0.0);
  double mean = sum / vect.size();
  std::vector<double> diff(vect.size());
  std::transform(vect.begin(), vect.end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  return std::make_pair(mean, std::sqrt(sq_sum / vect.size()));
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


typedef std::tuple<std::vector<float>, // gaussian interpolation
		   std::vector<float>, // mean
		   std::vector<float>, // standard deviation
		   uint16_t, // nb points in discus
		   uint16_t  // nb points in cell
		   > gaussianType;

gaussianType getGaussian(const std::vector<double>& pos,
			 const std::vector< CoordsList >& gridToInterpol,
			 const std::vector<size_t>& neighbors,
			 const size_t cellCol,
			 const size_t cellRow,
			 const size_t xSize,
			 const size_t ySize,
			 const float sigma,
			 const long int radius,
			 const float resolution,
			 const size_t nbBands,
			 const size_t nbPoints)
{

  // We take the coordinates at the center of the pixel
  Coords cellCoordsCenter = Coords{ static_cast<float>(cellCol + 0.5f),
                                    static_cast<float>(cellRow + 0.5f),
                                    -1.f };
  std::list < std::pair<float, Coords> > coordsList;
  uint16_t nbPointsInCell = 0;
  std::vector<size_t> indexes;
  std::vector<float> weights; 
  std::vector<float> gaussian_interp(nbBands, std::numeric_limits<double>::quiet_NaN());
  std::vector<float> mean(nbBands, std::numeric_limits<double>::quiet_NaN());
  std::vector<float> stdev(nbBands, 0);

  // allows retrocompatibility with cars: but it is a bug
  float dist_min = std::numeric_limits<float>::max();
  for(const auto neigh : neighbors){
    for(const auto& coords: gridToInterpol[neigh]){
      dist_min = std::min<float>(dist_min, distance(cellCoordsCenter, coords));
    }
  }

  for(const auto neigh : neighbors){
    for(const auto& coords: gridToInterpol[neigh]){
      float dist = distance(cellCoordsCenter, coords);
      if (dist < 0.5 + radius) {
	dist -= dist_min;
	dist *= resolution;
	weights.push_back(exp( (- (dist * dist) / (2 * sigma * sigma) )));
	indexes.push_back(coords[2]);

	if ((std::abs(coords[0] - (cellCol + 0.5f)) < 0.5f) &&
	    (std::abs(coords[1] - (cellRow + 0.5f)) < 0.5f)) {
	  nbPointsInCell++;
	}
      }
    }
  }

  uint16_t nbPointsInDisc = indexes.size(); 
  float weightsSum = std::accumulate(weights.begin(),
				     weights.end(),
				     0.0);

  if(nbPointsInDisc > 0){
    for( size_t band = 0; band < nbBands ; ++band) {
      std::vector<float> indexesValue(nbPointsInDisc);
      gaussian_interp[band] = 0;
      for( size_t point = 0; point < nbPointsInDisc ; ++point) {
	float weight = weights[point];
	float index = indexes[point];
	float value = pos[(band+2)*nbPoints+index];
	indexesValue[point] = value;
	gaussian_interp[band] += weight*value;
      }
      gaussian_interp[band] /= weightsSum;
     
      auto [mean_, stdev_] = vector_statistics(indexesValue);
      mean[band] = mean_;
      stdev[band] = stdev_;

    }
  }

  return std::make_tuple(gaussian_interp,
			 mean, stdev,
			 nbPointsInDisc,
			 nbPointsInCell);
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
				    std::vector<double> &meanVector,
				    std::vector<double> &stdevVector,
				    std::vector<uint16_t> &nbPointsInDiscVector,
				    std::vector<uint16_t> &nbPointsInCellVector)
{
  const size_t outSize = xSize*ySize;
  std::vector<double> outputVector(nbBands*outSize);
  meanVector.resize(nbBands*outSize);
  stdevVector.resize(nbBands*outSize);
  nbPointsInDiscVector.resize(nbBands*outSize);
  nbPointsInCellVector.resize(nbBands*outSize);

  // For each target pixel, we store a list of fractional Z coordinates
  std::vector< CoordsList > gridToInterpol(outSize);
  double x, y, col, row;
  size_t cellCol, cellRow, rowByCol;

  for ( size_t k = 0 ; k < nbPoints ; ++k ) {
    
    x = pos[(0*nbPoints)+k];
    y = pos[(1*nbPoints)+k];

    col = (x - xStart) / resolution;
    row = (yStart - y) / resolution;
    
    cellCol = floor(col);
    cellRow = floor(row);

    if ((cellCol >= 0) && (cellCol < xSize) && (cellRow >= 0) && (cellRow < ySize)) {
      rowByCol= cellCol + cellRow * xSize;
      gridToInterpol[rowByCol].push_back(Coords{col, row, static_cast<float> (k)});
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
    auto [out,
	  mean,
	  stdev,
	  nbPointsInDisc,
	  nbPointsInCell] = getGaussian(pos,
					gridToInterpol,
					neighbors,
					cellCol,
					cellRow,
					xSize,
					ySize,
					sigma,
					radius,
					resolution,
					nbBands,
					nbPoints);

    for( size_t band = 0; band < nbBands ; ++band) {
      outputVector[k*nbBands+band] = out[band];
      meanVector[k*nbBands+band] = mean[band];
      stdevVector[k*nbBands+band] = stdev[band];
    }

    nbPointsInDiscVector[k] = nbPointsInDisc;
    nbPointsInCellVector[k] = nbPointsInCell;
  }

  return outputVector;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

// wrap C++ function with NumPy array IO
typedef std::tuple<py::array, // gaussian interpolation
		   py::array, // mean
		   py::array, // standard deviation
		   py::array, // nb points in discus
		   py::array // nb points in cell
		   > pyPointCloudtoDSMType;

pyPointCloudtoDSMType pyPointCloudToDSM(py::array_t<double,
					py::array::c_style | py::array::forcecast> array,
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

  std::vector<double> meanVector;
  std::vector<double> stdevVector;
  std::vector<uint16_t> nbPointsInDiscVector;
  std::vector<uint16_t> nbPointsInCellVector;

  // call pure C++ function
  auto outputVector = pointCloudToDSM(pos,
				      nbBands,
				      nbPoints,
				      xStart,
				      yStart,
				      xSize,
				      ySize,
				      resolution,
				      radius,
				      sigma,
				      meanVector,
				      stdevVector,
				      nbPointsInDiscVector,
				      nbPointsInCellVector);

  // out, mean, stdev for each band + n_pts and n_in_cell
  ssize_t             ndim    = 3;
  std::vector<size_t> shape   = { xSize, ySize, nbBands };
  std::vector<size_t> strides = { nbBands*ySize*sizeof(double),
				  nbBands*sizeof(double),
				  sizeof(double)};
  
  // return 2-D NumPy array, I think here it is ok since the expected argument is
  // a pointer so there is no copy
  auto out = py::array(py::buffer_info(outputVector.data(),
				       sizeof(double), 
				       py::format_descriptor<double>::format(),
				       ndim,
				       shape,
				       strides
				       ));

  auto mean = py::array(py::buffer_info(meanVector.data(),
					sizeof(double), 
					py::format_descriptor<double>::format(),
					ndim,  
					shape,
					strides
					));

  auto stdev = py::array(py::buffer_info(stdevVector.data(),
					 sizeof(double), 
					 py::format_descriptor<double>::format(),
					 ndim,  
					 shape,
					 strides
					 ));

  ndim    = 2;
  shape   = { xSize, ySize };
  strides = { ySize*sizeof(uint16_t),
	      sizeof(uint16_t)};
  
  auto nbPointsInDisc = py::array(py::buffer_info(nbPointsInDiscVector.data(),
						    sizeof(size_t), 
						    py::format_descriptor<uint16_t>::format(),
						    ndim,  
						    shape,
						    strides
						    ));

  auto nbPointsInCell = py::array(py::buffer_info(nbPointsInCellVector.data(),
						  sizeof(size_t), 
						  py::format_descriptor<uint16_t>::format(),
						  ndim,  
						  shape,
						  strides
						  ));

  return std::make_tuple(out, mean, stdev, nbPointsInDisc, nbPointsInCell);

}

// wrap as Python module
PYBIND11_MODULE(rasterize, m)
{
  m.doc() = "rasterize";

  m.def("pc_to_dsm", &pyPointCloudToDSM, "Convert point cloud to digital surface model");
}
