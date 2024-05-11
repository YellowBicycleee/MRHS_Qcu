#pragma once

namespace qcu {
template <typename T, int N>
struct Matrix {
  T[N * N] = {};
  constexpr int rows() const { return N; }
  constexpr int cols() const { return N; }
  constexpr int size() const { return N * N; }

  Matrix() = default;
  Matrix(const Matrix<T, N> &) = default;
  Matrix(Matrix<T, N> &&) = default;
  Matrix &operator=(const Matrix<T, N> &) = default;
  Matrix &operator=(Matrix<T, N> &&) = default;

 private:
  __device__ __host__ inline int index(int i, int j) const { return i * N + j; }
};

}  // namespace qcu