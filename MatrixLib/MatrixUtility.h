/*!
	\file MatrixUtility.h
	\brief ������ ������ ��������� ��������� ��� ������ � ���������.

	 �� ������ ������ �� ������������ ��������� ���������:
	 Create,
	 ToMatrix,
	 ScalarProduct,
	 Transpose,
	 Tr,

*/

#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

#include "Matrix.h"
#include "MatrixTraits.h"

#define GetT(M, Trait) MatrixDecay<typename RowType<M>::T>::Trait

/*!
  \brief C������ �������� ������ ������� \f${rows x cols}\f$ ���� T (�������, ������-�������, ������-������).

  \details ������������ ��� ������������� ��������� ��� �������� ��������� ��������. 
  ��������� RowVector � ColVector ����� ����������� � ����� ����������, � �� ����� ��� Matrix ����� ����������� � ����� �����������,
  ��� ����� ������� �������� ��� ��������� ���������� ����. ������� Create ����� ��� ��������� � ����� ���������
  ����� ��������� �������: Matrix, RowVector, ColVector.

  \param[in] rows ����� ����� ���������� �������.
  \param[in] cols ����� �������� ���������� �������.

  \return ��������� ������ ������� \f${rows x cols}\f$ ���� T, ����������� ������.

  \throw std::exception ���� ��� �������-������ ���������� ����� ���� ������ 1 ��� ��� �������-������� ���������� �������� ���� ������ 1.
  �����, ���� ����� ����� ��� �������� ���� ������� �������.

*/

template<typename T, 
	typename = std::enable_if_t<(GetT(T, isMatrix)), int>
>
auto Create(size_t rows, size_t cols) {
	if (rows == 0 || cols == 0) throw std::exception(__FUNCSIG__" error: A matrix object should have at least one row and at least one column");
	if constexpr (GetT(T, isMatrixPolicy)) {
		return T(rows, cols);
	}
	else if constexpr (GetT(T, isRowVectorPolicy)) {
		if (rows != 1) throw std::exception(__FUNCSIG__" error: RowVector must have only one row.");
		return T(cols);
	}
	else if constexpr (GetT(T, isColVectorPolicy)) {
		if (cols != 1) throw std::exception(__FUNCSIG__" error: ColVector must have only one col.");
		return T(rows);
	}
}

/*!
  \brief ������ ��� ���������� ������� ������� ��� � �������. ���������� ������ �� ��������. ������������ ����� ������ ���� Matrix.

  \details �������� ���������� ������� - ��� ���������� ������� Matrix �� RowVector ��� ColVector.
  ������ ������� ������������� ����������� ��������� �������:
	Matrix -> Matrix
	RowVector -> Matrix
	ColVector -> Matrix
  ���������� ������������� ���������� ������� � ��� ������ ����� ������ ��, ��� � ������� ��������� �� ����.

  \param[in] v ��������� ������, �� �������� �������� Matrix.
  \return ������� �� ��������� ���� �� ����, ��� ��� � V.

*/
template<typename V, 
	typename = std::enable_if_t<(GetT(V, isMatrix)), int>
>
auto ToMatrix(const V& v) {
	typename ChangePolicy<MatrixIndexingPolicy, RowType<V>::T>::T ret(v.Rows(), v.Cols());
	for (size_t i = 0; i < v.Rows(); i++) for (size_t j = 0; j < v.Cols(); j++) ret.At(i, j) = v.At(i, j);
	return ret;
}


/*!
\brief ��������� ������������ ���� ��������.

\details ������������ ���� ������ ���� ���������. ���������� ������� (������-�������, ������-������) �� ������ ����.

\param[in] v1  ������ ������.
\param[in] v2  ����� ������.
\return ���������� ����� - ��������� ������������ ���� ��������. 
��� ������������� �������� - ��� ����� ��� ��� ��������� v1 � v2.

\throw std::exception ���� ������� �������� �� ���������.

*/

template<typename V1, typename V2,  
	typename = std::enable_if_t<(GetT(V1, isVectorPolicy), GetT(V2, isVectorPolicy)), int>
>
auto ScalarProduct(const V1& v1,const V2& v2) {
	if (v1.Size() != v2.Size()) throw std::exception("Vectors of different size cannot be multiplied");
	using C = std::common_type_t<GetT(V1, T), GetT(V2, T)>;
	C res = 0;
	for (size_t i = 0; i < v1.Size(); i++) {
		res += v1[i] * v2[i];
	}
	return res;
}


#define MatrixMult \
for (size_t i = 0; i < t1.Rows(); i++) for (size_t j = 0; j < t2.Cols(); j++) \
	ret.At(i, j) = ScalarProduct(t1.GetRow(i), t2.GetCol(j));\
return ret;\


template<typename T1, typename T2, 
	typename = std::enable_if_t<(GetT(T1, isMatrix) || GetT(T2, isMatrix)), int>
>
auto operator*(const T1& t1,const T2& t2) {
	if constexpr (
		GetT(T1, isMatrixPolicy) && GetT(T2, isMatrixPolicy) ||    //������� * ������� = �������
		GetT(T1, isColVectorPolicy) && GetT(T2, isRowVectorPolicy) //������-������� * ������-������ = �������
	) {
		if (t1.Cols() != t2.Rows()) throw std::exception("Matrices of this size cannot be multiplied");
		using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
		Matrix<C> ret(t1.Rows(), t2.Cols());

		MatrixMult
	}
	if constexpr (GetT(T1, isMatrixPolicy) && GetT(T2, isColVectorPolicy)) { //�������  * ������-������� = ������-�������
		if (t1.Cols() != t2.Rows()) throw std::exception("Matrices of this size cannot be multiplied");
		using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
		ColVector<C> ret(t1.Rows());

		MatrixMult
	}
	else if constexpr (GetT(T1, isRowVectorPolicy) && GetT(T2, isMatrixPolicy)) {//������-������ * ������� = ������-������.
		if (t1.Cols() != t2.Rows()) throw std::exception("Matrices of this size cannot be multiplied");
		using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
		RowVector<C> ret(t2.Cols());

		MatrixMult
	}
	else if constexpr (GetT(T1, isRowVectorPolicy) && GetT(T2, isColVectorPolicy)) {//������-������ * ������-������� = �����.
		return ScalarProduct(t1, t2);
	}
	else if constexpr (std::is_arithmetic_v<RowType<T1>::T> && GetT(T2, isMatrix)) { //����� * ������� = �������
		typedef std::common_type_t<typename RowType<T1>::T, GetT(T2, T)> C;
		typename ChangeT<C,typename RowType<T2>::T>::T ret = t2;
		for (size_t i = 0; i < t2.Rows(); i++) for (size_t j = 0; j < t2.Cols(); j++)
			ret.At(i, j) *= t1;
		return ret;
	}
	else if constexpr (GetT(T1, isMatrix) && std::is_arithmetic_v<RowType<T2>::T>) { //������� * ����� = �������.
		typedef std::common_type_t<GetT(T1, T), typename RowType<T2>::T> C;
		typename ChangeT<C,typename RowType<T1>::T>::T ret = t1;
		for (size_t i = 0; i < t1.Rows(); i++) for (size_t j = 0; j < t1.Cols(); j++)
			ret.At(i, j) *= t2;
		return ret;
	}
}


/*!
\brief ������������� ��������� ������.

\details ��������� �������� ������������� ��� ����������������.
	Matrix -> Matrix
	RowVector -> ColVector
	ColVector -> RowVector

\param[in] t ��������� ������.

\return ����������������� ��������� ������.

*/

template<typename T, 
	typename = std::enable_if_t<GetT(T, isMatrix), int>
>
auto Transpose(const T& t) {
	if constexpr (GetT(T, isMatrixPolicy)) {
		typename RowType<T>::T transposed(t.Cols(), t.Rows());
		for (size_t i = 0; i < t.Rows(); i++) for (size_t j = 0; j < t.Cols(); j++)
			transposed.At(j, i) = t.At(i, j);
		return transposed;
	}
	else if constexpr (GetT(T, isRowVectorPolicy)) {
		ColVector<GetT(T, T)> transposed(t.Size());
		for (size_t i = 0; i < t.Size(); i++)
			transposed[i] = t[i];
		return transposed;
	}
	else {
		RowVector<GetT(T, T)> transposed(t.Size());
		for (size_t i = 0; i < t.Size(); i++)
			transposed[i] = t[i];
		return transposed;
	}
}


template<typename T1, typename T2, 
	typename = std::enable_if_t<(GetT(T1, isMatrix) && GetT(T2, isMatrix)), int>
>
auto operator+(const T1& t1,const T2& t2) {
	if (t1.Rows() != t2.Rows() || t1.Cols() != t2.Cols()) throw std::exception(__FUNCSIG__" error: The matrices must be the same size.");
	using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
	using RetT = std::conditional_t<GetT(T1, isMatrixPolicy) || GetT(T2, isMatrixPolicy), 
		Matrix<C>,
		typename ChangeT<C, typename RowType<T1>::T>::T
	>;
	auto ret = Create<RetT>(t1.Rows(), t2.Cols());
	for (size_t i = 0; i < t1.Rows(); i++) for (size_t j = 0; j < t1.Cols(); j++) ret.At(i, j) = t1.At(i, j) + t2.At(i, j);
	return ret;
}


template<typename T1, typename T2,
	typename = std::enable_if_t<(GetT(T1, isMatrix) && GetT(T2, isMatrix)), int>
>
auto operator-(const T1& t1,const T2& t2) {
	if (t1.Rows() != t2.Rows() || t1.Cols() != t2.Cols()) throw std::exception(__FUNCSIG__" error: The matrices must be the same size.");
	using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
	using RetT = std::conditional_t<GetT(T1, isMatrixPolicy) || GetT(T2, isMatrixPolicy),
		Matrix<C>,
		typename ChangeT<C, typename RowType<T1>::T>::T
	>;
	auto ret = Create<RetT>(t1.Rows(), t2.Cols());
	for (size_t i = 0; i < t1.Rows(); i++) for (size_t j = 0; j < t1.Cols(); j++) ret.At(i, j) = t1.At(i, j) - t2.At(i, j);
	return ret;
}



/*!
\brief ������������ �������.

\details ����������� �������� ���������� �������.

\param[in] t1,t2 ��������� �������.

\throw std::exception ���� ������ ���������� ������ �� ���������.

\return �������, ���������� ��� ������������ ��������� ������ t1 � t2. 
��� ��������� ����������� ������� �������� ����� ����� ��� ��������� ���������� ������.

*/
template<typename T1,typename T2, 
	typename = std::enable_if_t<(GetT(T1, isMatrix) && GetT(T2, isMatrix)), int>
>
auto HadamardProduct(const T1& t1, const T2& t2) {
	//static_assert(GetT(T1, isMatrix) && GetT(T2, isMatrix));
	if (t1.Rows() != t2.Rows() || t1.Cols() != t2.Cols()) throw std::exception(__FUNCSIG__" error: Hadamard product cannot be applied to matrices of the differ size.");
	using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
	using RetT = std::conditional_t<GetT(T1, isMatrixPolicy) || GetT(T2, isMatrixPolicy),
		Matrix<C>,
		typename ChangeT<C, typename RowType<T1>::T>::T
	>;
	auto ret = Create<RetT>(t1.Rows(), t2.Cols());
	for (size_t i = 0; i < t1.Rows(); i++) for (size_t j = 0; j < t1.Cols(); j++) ret.At(i, j) = t1.At(i, j) * t2.At(i, j);
	return ret;

}

/*!
\brief ���� �������.

\details ���������� ����� ��������� ������� ���������.

\param[in] m ��������� ������.

\throw std::exception ���� ���������� ������ �� �������� ���������� ��������.

\return ���� �������. ��� ������������� �������� ��������� � ����� ��������� �������. 

*/
template<typename T,
	typename = std::enable_if_t<(GetT(T, isMatrix)), int>
>
auto Tr(const T& m) {
	if (m.Rows() != m.Cols()) throw std::exception(__FUNCSIG__" error: The matrix must be square.");
	using C = typename GetT(T, T);
	C tr = 0;
	for (size_t i = 0; i < m.Rows(); i++) {
		tr += m.At(i, i);
	}
	return tr;
}



namespace details {

	template<typename T, typename = std::enable_if_t<(GetT(T, isMatrix)), int> >
	void RowSum(T&& m, size_t row1, size_t row2, double k = 1) {
		if (row1 >= m.Rows() || row2 >= m.Rows()) throw std::exception(__FUNCSIG__" error: Index out of ragne.");

		for (size_t i = 0; i < m.Cols(); i++)
			m.At(row1, i) = m.At(row1, i) + k * m.At(row2, i);
	}


	template<typename T, typename = std::enable_if_t<(GetT(T, isMatrix)), int> >
	void SwapRows(T&& m, size_t row1, size_t row2) {
		if (row1 >= m.Rows() || row2 >= m.Rows()) throw std::exception(__FUNCSIG__" error: Index out of ragne.");
		if (row1 == row2) return;
		for (size_t i = 0; i < m.Cols(); i++)
			std::swap(m.At(row1, i), m.At(row2, i));
	}

	template<typename T, typename = std::enable_if_t<(GetT(T, isMatrix)), int> >
	void RowDivision(T&& m, size_t row, double k) {
		if (row >= m.Rows()) throw std::exception(__FUNCSIG__" error: Index out of ragne.");
		if (k == 0) throw std::exception(__FUNCSIG__" error: Division by zero");
		if (k == 1) return;
		for (size_t i = 0; i < m.Cols(); i++)
			m.At(row, i) /= k;
	}

}

/*!
\brief ����� ������.

\details �������� ������� ������� ������. ���������� ���������� �������.
����������� ����� ��������� ������������ ������� � ���� �������, ���� �������������� ��������� ���� ��������.

\param[in] m ��������� ������.
\param[out] det ������������ �������
\param[out] rank ���� �������.


\return �������, ���������� ������� ������.
*/
template<typename M, 
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
typename RowType<M>::T RowReduce(const M& m_, double* det = nullptr, size_t* rank = nullptr) {
	
	auto m = m_;

	if (m.Rows() != m.Cols()) det = nullptr;
	if (rank) *rank = 0;
	if (det) *det = 1;

	EntryPosition start_pos = { 0, 0 }; //�������, � ������� ���������� ����� �������� ��������.

	while (start_pos.col < m.Cols() && start_pos.row < m.Rows()) {
		EntryPosition lead_elem_i = start_pos; //������� �������
		size_t j = lead_elem_i.col;
		while (!m.At(lead_elem_i.row, lead_elem_i.col) && j < m.Cols())
		{
			for (size_t i = start_pos.row; i < m.Rows(); i++) {
				if (m.At(i, j) != 0) {
					lead_elem_i.row = i;
					lead_elem_i.col = j;
					break;
				}
			}
			j++;
		}
		if (!m.At(lead_elem_i.row, lead_elem_i.col)) {
			if (det) *det = 0;
			return m;
		}

		if (lead_elem_i.row != start_pos.row)
			if (det) *det *= -1;
		details::SwapRows(m, lead_elem_i.row, start_pos.row);
		lead_elem_i.row = start_pos.row;

		if (det) *det *= m.At(lead_elem_i.row, lead_elem_i.col);
		details::RowDivision(m, lead_elem_i.row, m.At(lead_elem_i.row, lead_elem_i.col));

		for (size_t i = 0; i < m.Rows(); i++)
			if (m.At(i, lead_elem_i.col) != 0 && i != lead_elem_i.row) {
				details::RowSum(m, i, lead_elem_i.row, -(double)m.At(i, lead_elem_i.col));
			}


		if (rank) *rank = start_pos.row + 1;
		start_pos.col = lead_elem_i.col + 1;
		start_pos.row++;
	}
	
	for (size_t i = 0; i < m.Cols() && i < m.Rows(); i++) {
		if (m.At(i, i) == 0) { 
			if (det) *det = 0;
			break;
		}
	}

	return m;

}

/*!
\brief ������������ �������.

\details ������� ������������ ������� � ������� ������ ������.

\param[in] m ��������� ������.

\return ������������ �������.
*/
template<typename M,
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
double Det(const M& m) {
	double det;
	RowReduce(m, &det);
	return det;
}

/*!
\brief ���� �������.

\details ������� ���� ������� � ������� ������ ������.

\param[in] m ��������� ������.

\return ���� �������.
*/
template<typename M, 
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
double Rg(const M& m) {
	size_t rank;
	RowReduce(m, nullptr, &rank);
	return rank;
}


/*!
\brief ��������� ����� �������.

\details ��������� ��������� ����� �������. 
 \f$ Norm(\vec{v}) = \sqrt{(\vec{v} , \vec{v})}\f$


\param[in] m ��������� ������.

\return ��������� ����� �������.
*/
template<typename V, 
	typename = std::enable_if_t<(GetT(V, isVectorPolicy)), int>
>
auto EuclidNorm(const V& m) {
	return sqrt(ScalarProduct(m, m));
}


/*!
\brief ������������ ����� �������.

\details ��������� ������������ ����� �������.
 \f$ Maxnorm(\vec{v}) = \max(|v_1|, |v_2|, \dots, |v_n|)\f$


\param[in] m ��������� ������.

\return ������������ ����� �������.
*/
template<typename V,
	typename = std::enable_if_t<(GetT(V, isVectorPolicy)), int>
>
auto MaxNorm(const V& m) {
	auto a = abs(m[0]);
	for (size_t i = 1; i < m.Size(); i++)
		if (a <= abs(m[i])) a = abs(m[i]);
	return a;
}


/*!
\brief ����� ����������.

\details ��������� ����� ����������.


\f$ FNorm(M) =  \sqrt{\sum\limits_{i = 0}^n \sum\limits_{j=0}^m M^2_{i,j}}\f$

\param[in] m �������.

\return ����� ����������.
*/
template<typename M,
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
auto FrobeniusNorm(const M& m) {
	double n = 0;
	for (size_t i = 0; i < m.Rows(); i++) for (size_t j = 0; j < m.Cols(); j++)
		n += m.At(i, j) * m.At(i, j);
	return sqrt(n);
}

#define DEG false
#define RAD true


/*!
\brief ���� ����� ���������.

\details ��������� ���� ����� ��������� �� �������.

\f$ Angle(\vec{a}, \vec{b}) = \frac{(\vec{a}, \vec{b})}{|\vec{a}||\vec{b}|}\f$

\param[in] v1 ������.
\param[in] v2 ������.
\param[in] rad ������ ������. 
���� �������� ��������� RAD, �� ����� ����� ��������� � ��������
���� �������� ��������� DEG, �� ����� ����� ��������� � ��������.


\return ���� ����� ��������� � �������� ��� � ��������.
*/
template<typename V1, typename V2, 
	typename = std::enable_if_t<(GetT(V1, isVectorPolicy) && GetT(V2, isVectorPolicy)), int>
>
auto Angle(const V1& v1,const V2& v2, bool rad = true) {
	return acos(ScalarProduct(v1, v2) / (EuclidNorm(v1) * EuclidNorm(v2))) * (rad ? 1 : 180 * M_1_PI );
}

/*!
\brief ����������.

\details ���������� ���������� ���������� �������.


\param[in] m �������.
\param[in] beg ������ ���� EntryPosition, ����������� ���������� ������ �������� ���� ����������.
\param[in] end ������ ���� EntryPosition, ����������� ���������� ������� ������� ���� ����������.

\return ���������� �� ��������� � ��������� �� beg �� end ������������.
*/
template<typename M, 
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
auto SubMatrix(const M& m, EntryPosition beg, EntryPosition end) {
	if (beg.col >= m.Cols() || beg.row >= m.Rows() || end.col >= m.Cols() || end.row >= m.Rows()) 
		throw std::exception(__FUNCSIG__" error: Index out of ragne.");
	if (beg.col > end.col || beg.row > end.row) 
		throw std::exception(__FUNCSIG__" error: The range was incorrect.");

	typename RowType<M>::T sub(end.row - beg.row + 1, end.col - beg.col + 1);

	for (size_t i = beg.row; i <= end.row; i++)
		for (size_t j = beg.col; j <= end.col; j++)
			sub.At(i - beg.row, j - beg.col) = m.At(i, j);
	return sub;
}

/*!
\brief ���������� ��� ������� � ����.

\details ���������� ��� ������� � ����, �������� �������� ������ ������� ������ �� ������.


\param[in] m1 �������.
\param[in] m2 �������.


\return �������, ���������� ������������ ����������.
*/
template<typename T1, typename T2, 
	typename = std::enable_if_t<(GetT(T1, isMatrix) && GetT(T2, isMatrix)), int>
>
auto Assign(const T1& m1,const T2& m2) {
	if (m1.Rows() != m2.Rows()) throw std::exception(__FUNCSIG__" error: Matrices with differ size of Cols cannot be assiged.");

	using C = std::common_type_t<GetT(T1, T), GetT(T2, T)>;
	Matrix<C> m(m1.Rows(), m1.Cols() + m2.Cols());
	for (size_t i = 0; i < m1.Rows(); i++) for (size_t j = 0; j < m1.Cols(); j++)
		m.At(i, j) = m1.At(i, j);
	for (size_t i = 0; i < m2.Rows(); i++) for (size_t j = 0; j < m2.Cols(); j++)
		m.At(i, j + m1.Cols()) = m2.At(i, j);
	return m;
}

/*!
\brief ��������� �������� �������
\details ��������� �������� ������� � ������� ������ ������.


\param[in] m �������.

\throw std::exception ���� ���������� ������� ���� �� ���������� ��� ���� �������� ������� � ������ �� ����������.

\return �������� �������.
*/
template<typename M, typename = std::enable_if_t<(GetT(M, isMatrix)), int>>
auto Inverse(const M& m) {
	if (m.Rows() != m.Cols()) std::exception(__FUNCSIG__" error: Inverse matrix cannot be found for not square matrix.");
	auto reduced = RowReduce(Assign(m, Identity<GetT(M, T)>(m.Rows())));

	auto r = Rg(SubMatrix(m, EntryPosition{ 0, 0 }, EntryPosition{ m.Rows() - 1, m.Cols() - 1 }));

	if (r != m.Rows()) throw std::exception(__FUNCSIG__" error: There is no inverse matrix 'cos of the rank is not equal to amount of Rows.");
	return SubMatrix(reduced, EntryPosition{ 0, m.Cols() }, EntryPosition{ m.Rows() - 1, 2 * m.Cols() - 1});
}




template<typename M,
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
std::fstream& operator>>(std::fstream& out, M& m) {
	m.ReadFromFStream(out);
	return out;
}

template<typename M,
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
std::fstream& operator<<(std::fstream& out,const M& m) {
	m.Save(out);
	return out;
}

template<typename M,
	typename = std::enable_if_t<(GetT(M, isMatrix)), int>
>
auto Generate(size_t rows, size_t cols, size_t left_b = 0, size_t right_b = 100) {
	typename RowType<M>::T m(rows, cols);
	for (size_t i = 0; i < rows; i++) for (size_t j = 0; j < cols; j++)
		m.At(i, j) = rand() % (right_b - left_b + 1) + left_b;
	return m;

}