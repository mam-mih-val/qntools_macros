#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cmath>

template<typename T, size_t N>
class Matrix{
public:
  Matrix() = default;
  ~Matrix() = default;
  Matrix(std::array< std::array< T, N>, N> matrix) : matrix_( std::move(matrix) ){ };
  Matrix( const Matrix<T, N>& ) = default;
  Matrix( Matrix<T, N>&& ) = default;
  auto operator=( const Matrix<T, N>& ) -> Matrix<T,N>& = default;
  auto operator=( Matrix<T, N>&& ) -> Matrix<T,N>& = default;

  auto operator[]( size_t i ) -> std::array<T, N>& { return matrix_[i]; }

  template<typename U, size_t M>
  friend auto Reduce( const Matrix<U, M>& A, size_t col, size_t row ) -> Matrix<U, M-1>;

  template<typename U, size_t M>
  friend auto Det( const Matrix<U, M>& A ) -> U;

  template<typename U>
  friend auto Det( const Matrix<U, 1>& A ) -> U;

  template<typename U, size_t M>
  friend auto Transpose( const Matrix<U, M>& A ) -> Matrix<U, M>;

  template<typename U, size_t M>
  friend auto Inverse( const Matrix<U, M>& A ) -> Matrix<U, M>;

private:
  std::array< std::array< T, N>, N> matrix_;
};

template<typename U, size_t M>
auto Reduce( const Matrix<U, M>& A, size_t col, size_t row ) ->  Matrix<U, M-1>{
  static_assert( M > 1 );
  Matrix<U, M-1> result;
  for( auto i=size_t{0}, res_i=size_t{0}; i<M; ++i ){
    if( i == col )
      continue;
    for( auto j=size_t{0}, res_j = size_t{0}; j<M; ++j ){
      if( j == row )
        continue;
      result.matrix_[res_i][res_j] = A.matrix_[i][j];
      ++res_j;
    }
    ++res_i;
  }
  return result;
}

template<typename U, size_t M>
auto Det( const Matrix<U, M>& A ) -> U{
  auto result =  A.matrix_[0][0] * Det( Reduce(A, 0, 0) );
  for( auto i=size_t{1}; i<M; ++i ){
    result += A.matrix_[i][0] * Det(Reduce( A, i, 0)) * pow( -1, i );
  }
  return result;
}

template<typename U>
auto Det( const Matrix<U, 1>& A ) -> U{
  return A.matrix_[0][0];
}

template<typename U, size_t M>
auto Transpose( const Matrix<U, M>& A ) -> Matrix<U, M>{
  auto res = Matrix<U, M>{};
  for( auto i=size_t{0}; i<M; ++i ){
    for( auto j = size_t{0}; j<M; ++j ){
      res.matrix_[i][j] = A.matrix_[j][i];
    }
  }
  return res;
}

template<typename U, size_t M>
auto Inverse( const Matrix<U, M>& A ) -> Matrix<U, M>{
  static_assert(M > 1);
  auto result = Matrix<U, M>{};
  auto AT = Transpose(A);
  auto detA = Det(A) ;
  for( auto i=size_t{0}; i<M; ++i ){
    for( auto j = size_t{0}; j<M; ++j ){
      auto Mij = pow(-1, i+j) * Det(Reduce(A, i, j));
      result.matrix_[i][j] = Mij / detA;
    }
  }
  return result;
}


#endif