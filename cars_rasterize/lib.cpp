#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <immintrin.h>  // AVX2 intrinsics

// -------------
// pure C++ code
// -------------

std::vector<float> pointCloudToDSM(const double* pointsVector,
				   const double* valuesVector,
				   const int* validVector,
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
  
  // Pre-compute constants
  const double radiusSq =  (0.5 + radius) * (0.5 + radius);
  const double resolutionSq = resolution * resolution;
  const double sigmaSq2 = - 2 * sigma * sigma; 
  const double resSigFactor = resolutionSq / sigmaSq2;

  // For each target pixel, we update the internal vectors in all neighboring pixels that it applies.
  for ( long int pointIdx = 0 ; pointIdx < nbPoints ; ++pointIdx ) {
    // Verify if the point is valid
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

        double distSq = ( (pointCol - neighColCenter) * (pointCol - neighColCenter) +
                          (pointRow - neighRowCenter) * (pointRow - neighRowCenter) );
        if (distSq < radiusSq) {
          // Compute weight
          double weight = exp(distSq * resSigFactor);
          weightsSumVector[neighIdx] += weight;
          // Count points in cell and in disc
          nbPointsInDiscVector[neighIdx]++;
          if ((std::abs(neighColCenter - (pointCol)) < 0.5) &&
              (std::abs(neighRowCenter - (pointRow)) < 0.5)) {
            nbPointsInCellVector[neighIdx]++;
          }
          
          // Compute mean and variance (Welford's algorithm)
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
    // Test if there is no point, if so fill with NaN
    if (nbPointsInDiscVector[k] == 0) {
      for (long int band = 0; band < nbBands; ++band) {
        outputVector[k * nbBands + band] = std::numeric_limits<float>::quiet_NaN();
        meanVector[k * nbBands + band] = std::numeric_limits<float>::quiet_NaN();
        stdevVector[k * nbBands + band] = std::numeric_limits<float>::quiet_NaN();
      }
      continue;
    }

    // Fill expected data using AVX2 vectorization when possible
    const float weightsSum = weightsSumVector[k];
    const float nbPointsInDisc = static_cast<float>(nbPointsInDiscVector[k]);
    const __m128 weightsSum_vec = _mm_set1_ps(weightsSum);
    const __m128 nbPointsInDisc_vec = _mm_set1_ps(nbPointsInDisc);
    
    long int band = 0;
    // Process 4 bands at a time with SSE/AVX
    for (; band + 3 < nbBands; band += 4) {
      const long int baseIdx = k * nbBands + band;
      
      // Load 4 bandSumWeighted values and convert from double to float
      __m256d bandSumWeighted_d = _mm256_loadu_pd(&bandSumWeighted[baseIdx]);
      __m128 bandSumWeighted_vec = _mm256_cvtpd_ps(bandSumWeighted_d);
      
      // Load 4 bandM2 values and convert from double to float  
      __m256d bandM2_d = _mm256_loadu_pd(&bandM2[baseIdx]);
      __m128 bandM2_vec = _mm256_cvtpd_ps(bandM2_d);
      
      // Compute outputVector: bandSumWeighted / weightsSum
      __m128 output_vec = _mm_div_ps(bandSumWeighted_vec, weightsSum_vec);
      _mm_storeu_ps(&outputVector[baseIdx], output_vec);
      
      // Compute stdevVector: sqrt(bandM2 / nbPointsInDisc)
      __m128 variance_vec = _mm_div_ps(bandM2_vec, nbPointsInDisc_vec);
      __m128 stdev_vec = _mm_sqrt_ps(variance_vec);
      _mm_storeu_ps(&stdevVector[baseIdx], stdev_vec);
    }
    
    // Handle remaining bands (scalar fallback)
    for (; band < nbBands; ++band) {
      const long int bandIdx = k * nbBands + band;
      outputVector[bandIdx] = static_cast<float>(bandSumWeighted[bandIdx]) / weightsSum;
      stdevVector[bandIdx] = std::sqrt(static_cast<float>(bandM2[bandIdx]) / nbPointsInDisc);
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
  const double* pointsPtr = points.data();
  const double* valuesPtr = values.data();
  const int* validPtr = valid.data();

  std::vector<float> weightsSumVector;
  std::vector<float> meanVector;
  std::vector<float> stdevVector;
  std::vector<uint16_t> nbPointsInDiscVector;
  std::vector<uint16_t> nbPointsInCellVector;

  // call pure C++ function
  auto outputVector = pointCloudToDSM(pointsPtr,
				      valuesPtr,
				      validPtr,
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
