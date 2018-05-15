//
// Created by Merrill, True on 5/14/18.
//

#ifndef QUANTUM_FDTD_WAVEFUNCTION_H
#define QUANTUM_FDTD_WAVEFUNCTION_H

#include <complex>
#include <vector>
#include <cstdlib>


namespace qfdtd {

/* Fundamental constants that define the atomic units system */
static const double hbar = 1;      //! Reduced Planck's constant
static const double me = 1;        //! Electron mass
static const double qe = 1;        //! Electron charge
static const double ke = 1;        //! Coulomb's constant

/**
 * @brief Class representing a wavefunction in the space domain.  The
 *  wavefunction is sampled on a hypergrid of points.  On each grid point we
 *  store the complex amplitude.
 */
class Wavefunction {
 private:
  std::vector<unsigned> m_shape;          //! Shape of the hypergrid
  std::vector<double> m_lower_bounds;   //! Hypercube lower bounds
  std::vector<double> m_upper_bounds;   //! Hypercube upper bounds
  std::vector<double> m_stride;         //! Hypercube stride
  std::vector<std::complex<double>> m_psi;  //! Wavefunction data

 public:

  /**
   * @brief Create a Wavefunction object
   *
   * @param shape Hypercube shape
   * @param lower Lower hypercube bounds
   * @param upper Upper hypercube bounds
   */
  Wavefunction(std::vector<unsigned> &shape, std::vector<double> &lower,
               std::vector<double> &upper);


  /* Getters */
  const std::vector<unsigned> &getShape(void);
  const std::vector<double> &getLowerBounds(void);
  const std::vector<double> &getUpperBounds(void);
  const std::vector<double> &getStride(void);

  /**
   * @brief Multi-array like indexing for wavefunction values
   *
   * @param indices vector of array indices
   * @return value of the wavefunction
   */
  std::complex<double> &operator()(const std::vector<int> &indices);
  const std::complex<double> &operator()(const std::vector<int> &indices) const;

  /**
   * @brief Returns position coordinates corresponding to `indices`
   *
   * @param indices vector array of indices
   * @return vector of position coordinates
   */
  std::vector<double> position(const std::vector<int> &indices);

  /**
   * @brief Normalize the wavefunction
   */
  void normalize(void);
};

}

#endif //QUANTUM_FDTD_WAVEFUNCTION_H
