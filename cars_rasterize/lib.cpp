#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>

// -------------
// pure C++ code
// -------------

std::vector<double> pc_to_dsm(const std::vector<double>& pos,
			      const size_t nb_bands, const size_t nb_points,
			      const float xstart, const float ystart, const size_t xsize, const size_t ysize,
			      const float resolution,
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

  for ( size_t k = 0 ; k < nb_points ; ++k ) {
    const float x = pos[(1*nb_points)+k];
    const float y = pos[(2*nb_points)+k];
    const float z = pos[(3*nb_points)+k];

    const ssize_t i = std::lround((x - xstart) / resolution);
    const ssize_t j = std::lround((ystart - y) / resolution);

    if (trace) {
      std::cout << "pt " << k << ": "
		<< "x, y, z > " 
		<< x << ", " << y << ", " << z << std::endl;
      std::cout << "i, j > " << i << ", " << j << std::endl;
    }

    if ((i >= 0) && (i < xsize) && (j >= 0) && (j < ysize)) {
      const size_t ij = i+j*xsize;
      output[ij] = z;
    }
  }

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
		       bool trace)
{
  // check input dimensions
  if ( array.ndim()     != 2 )
    throw std::runtime_error("Input should be 2-D NumPy array");

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> pos(array.size());

  // copy py::array -> std::vector
  std::memcpy(pos.data(),array.data(),array.size()*sizeof(double));

  size_t nb_bands = array.shape()[0]-2;
  size_t nb_points = array.shape()[1];

  // call pure C++ function
  std::vector<double> result = pc_to_dsm(pos, nb_bands, nb_points,
					 xstart, ystart, xsize, ysize, resolution,
					 trace);

  ssize_t             ndim    = 2;
  std::vector<size_t> shape   = { xsize, ysize };
  std::vector<size_t> strides = { sizeof(double)*ysize, sizeof(double) };

  // return 2-D NumPy array
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
