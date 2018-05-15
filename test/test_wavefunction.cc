#include "wavefunction.h"
#include "gtest/gtest.h"

TEST(wavefunction, constructor)
{
  std::vector<unsigned> shape = {101, 101, 11};
  std::vector<double> lower = {-100, -100, -10};
  std::vector<double> upper = {100, 100, 10};
  qfdtd::Wavefunction psi = qfdtd::Wavefunction(shape, lower, upper);

  const std::vector<double> &stride = psi.getStride();
  EXPECT_NEAR(stride[0], stride[1], 1e-4);
  EXPECT_NEAR(stride[0], stride[2], 1e-4);
}

TEST(wavefunction, index_access)
{
  std::vector<unsigned> shape = {101, 101};
  std::vector<double> lower = {-100, -100};
  std::vector<double> upper = {100, 100};
  qfdtd::Wavefunction psi = qfdtd::Wavefunction(shape, lower, upper);

  std::complex<double> v(0, -1);
  psi({0, 0}) = v;
  EXPECT_EQ(psi({0, 0}), v);
}

TEST(wavefunction, position)
{
  std::vector<unsigned> shape = {101, 101};
  std::vector<double> lower = {-100, -100};
  std::vector<double> upper = {100, 100};
  qfdtd::Wavefunction psi = qfdtd::Wavefunction(shape, lower, upper);

  std::vector<double> x = psi.position({0, 0});
  EXPECT_EQ(x, lower);
  std::vector<double> y = psi.position({-1, -1});
  EXPECT_EQ(y, upper);
}