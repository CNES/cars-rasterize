#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>

// -------------
// pure C++ code
// -------------

using Coords = std::array<double, 3>;
using CoordsVector = std::vector<Coords>;

double distance( const Coords& coords1,
		 const Coords& coords2 )
{
  // Euclidean square root distance
  return sqrt( (coords1[0] - coords2[0]) * (coords1[0] - coords2[0]) +
               (coords1[1] - coords2[1]) * (coords1[1] - coords2[1]) );
}

std::vector<float> pointCloudToDSM(const std::vector<double>& pointsVector,
				   const std::vector<double>& valuesVector,
				   const std::vector<int>& validVector,
				   const long int nbBands,
				   const long int nbPoints,
				   const double xStart,
				   const double yStart,
				   const long int xSize,
				   const long int ySize,
				   const double resolution,
				   const double radius,
				   const double sigma,
				   std::vector<float>& weightsSumVector,
				   std::vector<float>& meanVector,
				   std::vector<float>& stdevVector,
				   std::vector<uint16_t>& nbPointsInDiscVector,
				   std::vector<uint16_t>& nbPointsInCellVector)
{
  const long int outSize = xSize*ySize;
  std::vector<float> outputVector(nbBands*outSize);
  weightsSumVector.resize(outSize);
  meanVector.resize(nbBands*outSize);
  stdevVector.resize(nbBands*outSize);
  nbPointsInDiscVector.resize(outSize);
  nbPointsInCellVector.resize(outSize);

  // Initialize output vectors
  std::fill(weightsSumVector.begin(), weightsSumVector.end(), 0.0f);
  std::fill(meanVector.begin(), meanVector.end(), 0.0f);
  std::fill(stdevVector.begin(), stdevVector.end(), 0.0f);
  std::fill(nbPointsInDiscVector.begin(), nbPointsInDiscVector.end(), 0);
  std::fill(nbPointsInCellVector.begin(), nbPointsInCellVector.end(), 0);

  // Internal computation vectors for weighted sum and variance (Welford's algorithm)
  // This allow to perform only one pass over the points
  std::vector<double> bandSumWeighted(outSize*nbBands); 
  std::vector<double> bandM2(outSize*nbBands);          

  // For each target pixel, we update the stats data in all neighboorg that it apply.
  for ( long int pointIdx = 0 ; pointIdx < nbPoints ; ++pointIdx ) {
    if(validVector[pointIdx] == 0) {
      continue;
    }

    double x = pointsVector[(0*nbPoints)+pointIdx];
    double y = pointsVector[(1*nbPoints)+pointIdx];

    double pointCol = (x - xStart) / resolution;
    double pointRow = (yStart - y) / resolution;

    long int pointCellCol = floor(pointCol);
    long int pointCellRow = floor(pointRow);

    // compute min/max indexes of neighboring cells to avoid borders
    long int minNeighCol = std::max<long int>(0, pointCellCol - static_cast<long int>(std::ceil(radius)));
    long int maxNeighCol = std::min<long int>(xSize - 1, pointCellCol + static_cast<long int>(std::ceil(radius)));
    long int minNeighRow = std::max<long int>(0, pointCellRow - static_cast<long int>(std::ceil(radius)));
    long int maxNeighRow = std::min<long int>(ySize - 1, pointCellRow + static_cast<long int>(std::ceil(radius)));

    // Iterate over all neightbor
    for (long int neighRow = minNeighRow; neighRow <= maxNeighRow; ++neighRow) {
      for (long int neighCol = minNeighCol; neighCol <= maxNeighCol; ++neighCol) {
        long int neighIdx = neighRow * xSize + neighCol;
        double neighColCenter = static_cast<double>(neighCol) + 0.5;
        double neighRowCenter = static_cast<double>(neighRow) + 0.5;

        double dist = (std::sqrt( (pointCol - neighColCenter) * (pointCol - neighColCenter) +
                                  (pointRow - neighRowCenter) * (pointRow - neighRowCenter) ));
        if (dist < 0.5 + radius) {
          dist *= resolution;
          double weight = exp( (- (dist * dist) / (2 * sigma * sigma) ));
          weightsSumVector[neighIdx] += weight;
          if ((std::abs(neighColCenter - (pointCol)) < 0.5f) &&
              (std::abs(neighRowCenter - (pointRow)) < 0.5f)) {
            nbPointsInCellVector[neighIdx]++;
          }
          nbPointsInDiscVector[neighIdx]++;
          
          for( long int band = 0; band < nbBands ; ++band) {
            const long int bandIdx = neighIdx * nbBands + band;
            double value = valuesVector[band*nbPoints+pointIdx];
            double delta = value - meanVector[bandIdx];
            double newMean = meanVector[bandIdx] + delta/nbPointsInDiscVector[neighIdx];
            meanVector[bandIdx] = newMean;
            bandM2[bandIdx] += delta * (value - newMean);
            bandSumWeighted[bandIdx] += value*weight;
          }
        }
      }
    }
  }

  // Iterate over all output pixel to fill output
  for ( long int k = 0 ; k < outSize ; ++k ) {
    // Test if there is no point
    if (nbPointsInDiscVector[k] == 0) {
      for (long int band = 0; band < nbBands; ++band) {
        outputVector[k * nbBands + band] = std::numeric_limits<float>::quiet_NaN();
        meanVector[k * nbBands + band] = std::numeric_limits<float>::quiet_NaN();
        stdevVector[k * nbBands + band] = std::numeric_limits<float>::quiet_NaN();
      }
      continue;
    }

    // Fill expected data
    //   weightsSumVector[k] already set
    //   nbPointsInDiscVector[k] already set
    //   nbPointsInCellVector[k] already set
    for (long int band = 0; band < nbBands; ++band) {
      outputVector[k * nbBands + band] = bandSumWeighted[k * nbBands + band] / weightsSumVector[k];
      // meanVector[k * nbBands + band] already set
      stdevVector[k * nbBands + band] = std::sqrt(bandM2[k * nbBands + band] / nbPointsInDiscVector[k]);
    }
  }
  return outputVector;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

// wrap C++ function with NumPy array IO
typedef std::tuple<py::array, // gaussian interpolation
		   py::array, // sum of weights
		   py::array, // mean
		   py::array, // standard deviation
		   py::array, // nb points in discus
		   py::array // nb points in cell
		   > pyPointCloudtoDSMType;

pyPointCloudtoDSMType pyPointCloudToDSM(py::array_t<double,
					py::array::c_style | py::array::forcecast> points,
					py::array_t<double,
					py::array::c_style | py::array::forcecast> values,
					py::array_t<int,
					py::array::c_style | py::array::forcecast> valid,
					double xStart,
					double yStart,
					size_t xSize,
					size_t ySize,
					double resolution,
					double radius,
					double sigma)
{
  // check input dimensions
  if ( points.ndim()     != 2 )
    throw std::runtime_error("points should be 2-D NumPy array");

  if ( points.shape()[0] != 2 )
    throw std::runtime_error("points should have size [2, N]");

  if ( values.ndim()     != 2 )
    throw std::runtime_error("values should be 2-D NumPy array");

  if ( values.shape()[1] != points.shape()[1] ) {
    std::string message = "values and points should have the same second dimension";
    message = message + std::string("(") + std::to_string(values.shape()[1]);
    message = message + std::string(" ~= ") + std::to_string(points.shape()[1]) + std::string(")");
    throw std::runtime_error(message);
  }

  if ( valid.ndim()     != 2 )
    throw std::runtime_error("valid should be 2-D NumPy array");

  if ( valid.shape()[1] != points.shape()[1] ) {
    std::string message = "valid and points should have the same second dimension";
    message = message + std::string("(") + std::to_string(valid.shape()[1]);
    message = message + std::string(" ~= ") + std::to_string(points.shape()[1]) + std::string(")");
    throw std::runtime_error(message);
  }

  size_t nbBands = values.shape()[0];
  size_t nbPoints = points.shape()[1];

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> pointsVector(points.size());
  std::vector<double> valuesVector(values.size());
  std::vector<int> validVector(valid.size());

  // copy py::array -> std::vector
  // Do we need to make a copy here ? Is there something similar
  // to cython memory view
  // https://stackoverflow.com/questions/54793539/pybind11-modify-numpy-array-from-c
  // It seems that passing the py::array by reference should do the job if
  // there is not conflicting type, however the function pointCloudToDSM needs a double * ptr as
  // input instead of a std::vector.
  std::memcpy(pointsVector.data(),points.data(),points.size()*sizeof(double));
  std::memcpy(valuesVector.data(),values.data(),values.size()*sizeof(double));
  std::memcpy(validVector.data(),valid.data(),valid.size()*sizeof(int));

  std::vector<float> weightsSumVector;
  std::vector<float> meanVector;
  std::vector<float> stdevVector;
  std::vector<uint16_t> nbPointsInDiscVector;
  std::vector<uint16_t> nbPointsInCellVector;

  // call pure C++ function
  auto outputVector = pointCloudToDSM(pointsVector,
				      valuesVector,
				      validVector,
				      nbBands,
				      nbPoints,
				      xStart,
				      yStart,
				      xSize,
				      ySize,
				      resolution,
				      radius,
				      sigma,
				      weightsSumVector,
				      meanVector,
				      stdevVector,
				      nbPointsInDiscVector,
				      nbPointsInCellVector);

  // out, mean, stdev for each band + n_pts and n_in_cell
  size_t              ndim    = 3;
  std::vector<size_t> shape   = { xSize, ySize, nbBands };
  std::vector<size_t> strides = { nbBands*ySize*sizeof(float),
				  nbBands*sizeof(float),
				  sizeof(float)};

  // return 2-D NumPy array, I think here it is ok since the expected argument is
  // a pointer so there is no copy
  auto out = py::array(py::buffer_info(outputVector.data(),
				       sizeof(float),
				       py::format_descriptor<float>::format(),
				       ndim,
				       shape,
				       strides
				       ));

  auto mean = py::array(py::buffer_info(meanVector.data(),
					sizeof(float),
					py::format_descriptor<float>::format(),
					ndim,
					shape,
					strides
					));

  auto stdev = py::array(py::buffer_info(stdevVector.data(),
					 sizeof(float),
					 py::format_descriptor<float>::format(),
					 ndim,
					 shape,
					 strides
					 ));

  ndim    = 2;
  shape   = { xSize, ySize };
  strides = { ySize*sizeof(float),
	      sizeof(float)};

  auto weightsSum = py::array(py::buffer_info(weightsSumVector.data(),
					      sizeof(float),
					      py::format_descriptor<float>::format(),
					      ndim,
					      shape,
					      strides
					      ));

  ndim    = 2;
  shape   = { xSize, ySize };
  strides = { ySize*sizeof(uint16_t),
	      sizeof(uint16_t)};

  auto nbPointsInDisc = py::array(py::buffer_info(nbPointsInDiscVector.data(),
						  sizeof(uint16_t),
						  py::format_descriptor<uint16_t>::format(),
						  ndim,
						  shape,
						  strides
						  ));

  auto nbPointsInCell = py::array(py::buffer_info(nbPointsInCellVector.data(),
						  sizeof(uint16_t),
						  py::format_descriptor<uint16_t>::format(),
						  ndim,
						  shape,
						  strides
						  ));

  return std::make_tuple(out, weightsSum, mean, stdev, nbPointsInDisc, nbPointsInCell);

}

// wrap as Python module
PYBIND11_MODULE(rasterize, m)
{
  m.doc() = "rasterize";

  m.def("pc_to_dsm", &pyPointCloudToDSM, "Convert point cloud to digital surface model",
	py::arg("points"),
	py::arg("values"),
	py::arg("valid"),
	py::arg("xstart"),
	py::arg("ystart"),
	py::arg("xsize"),
	py::arg("ysize"),
	py::arg("resolution"),
	py::arg("radius"),
	py::arg("sigma")
	);
}
