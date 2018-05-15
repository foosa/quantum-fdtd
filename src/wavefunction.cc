//
// Created by Merrill, True on 5/14/18.
//

#include "wavefunction.h"


qfdtd::Wavefunction::Wavefunction(std::vector<unsigned> &shape,
                                  std::vector<double> &lower,
                                  std::vector<double> &upper) :
  m_stride (std::vector<double>(shape.size()))
{
  m_shape = shape;
  m_lower_bounds = lower;
  m_upper_bounds = upper;

  /* Compute the stride */
  unsigned dim = shape.size();
  unsigned num_points = 1;
  for (unsigned i = 0; i < dim; i++)
  {
    if (m_lower_bounds[i] == m_upper_bounds[i])
    {
      /* Lower and upper bounds coincide, so therefore this dimension is
       * removed from the problem */
      m_stride[i] = 0;
      m_shape[i] = 1;
    }

    else if (m_shape[i] == 1) {
      /* Only one sample point in this direction, so this dimension is
       * removed from the problem */
      m_upper_bounds[i] = m_lower_bounds[i];
    }

    else {
      m_stride[i] = (m_upper_bounds[i] - m_lower_bounds[i]) /
              (m_shape[i] - 1);
    }

    num_points *= m_shape[i];
  }

  /* Initialize the data for the wavefunction */
  m_psi = std::vector<std::complex<double>>(num_points);
}

const std::vector<unsigned>& qfdtd::Wavefunction::getShape()
{
  return m_shape;
}

const std::vector<double>& qfdtd::Wavefunction::getLowerBounds()
{
  return m_lower_bounds;
}

const std::vector<double>& qfdtd::Wavefunction::getUpperBounds()
{
  return m_upper_bounds;
}

const std::vector<double>& qfdtd::Wavefunction::getStride()
{
  return m_stride;
}

unsigned wrap_index(int index, unsigned bound)
{
  unsigned wrapped = index;
  if (index > bound)
    wrapped = index % bound;

  if (index < 0)
  {
    unsigned reversed = wrap_index(-index, bound);
    wrapped = bound - reversed;
  }

  return wrapped;
}

/**
 * @brief Computes the row-major ordered index
 *
 * @param indices
 * @param shape
 * @return
 */
unsigned row_major_index(const std::vector<int> &indices,
                         const std::vector<unsigned> &shape)
{
  unsigned index = 0;
  unsigned term;
  unsigned wrapped;

  for (unsigned k = 0; k < indices.size(); k++)
  {
    term = 1;
    for (unsigned l = k + 1; l < indices.size(); l++)
      term *= shape[l];

    wrapped = wrap_index(indices[k], shape[k]);
    index += term * wrapped;
  }

  return index;
}

std::complex<double> &qfdtd::Wavefunction::operator()(
        const std::vector<int> &indices)
{
  unsigned index = row_major_index(indices, m_shape);
  return m_psi[index];
}

const std::complex<double>& qfdtd::Wavefunction::operator()(
        const std::vector<int> &indices) const
{
  unsigned index = row_major_index(indices, m_shape);
  return m_psi[index];
}

std::vector<double> qfdtd::Wavefunction::position(
        const std::vector<int> &indices)
{
  std::vector<double> x(m_shape.size());
  for (unsigned k = 0; k < m_shape.size(); k++)
  {
    unsigned idx = wrap_index(indices[k], m_shape[k]);
    x[k] = m_lower_bounds[k] + (m_stride[k] * idx);
  }

  return x;
}
