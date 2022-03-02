/*
�� ������ ������������ ��������, ���� std::vector ��� ������� ������������ ���������� ����������� ���� ������ ��� ������������� ������ � �������� ������������ �����������.

������� ����� �������. (Yes)
������ ������������� ��� ������� ������� 1xN ��� Nx1.(Yes)

�������� ��������� ��� ������ ��������� ����� - ���������, ������������, ������� � ������ ����������� �������, ������������ �������. (Yes)

������������ ���������� ���������� ��� ���������� �������� ��������, ��������� � ��������� �������� � ������. (Yes)

��� ���������� ���������� ������� ������� ����������� ������� � ���������� � ��������� ����� (h � cpp).

������������ �������� ���������� ��� ��������� ��������� ��������.(Yes)

����������� �������� ��� ��������� � ���������: (Yes)

�������� ������, ��������� ������ � ��������� ������ �� ����� (Yes)
��������� ������ (Yes)
������������ ������� (Yes)

���������� �������� ������ � ����� (<<) (YEs)
*/

/*
��������� ���������� �������, ������������� � ������ ��1. ����������� �������� ��� ��������� � ���������:
���� ������� (Yes)

������������ ������� (������� ������) (Yes)

��������� ������������ �������� (Yes)

����� ������� (��������� �����, ������������ �����) (Yes)

����� ������� (����� ����������) (Yes)
*/

/*
��������� ���������� �������, ������������� � ������ ��1 � ��2. ����������� �������� ��� ��������� � ���������:

���� ����� ��������� (Yes)

���� ������� (� ������� ��������� ������) (Yes)

�������� ������� (���� ����������) (Yes)

���������������� ������� (Yes)
*/

/*
��������� ���������� �������, ������������� � ������ ��1-��3.
����������� ���������������� ��� �������� ������ � �������� �� ����� � ���������� ������ � �������� � ����. (Yes)
����������� 2 ������ ������: ������/������ � ��������� � �������� ����. (Yes)
��� ���������� ������ ������������ ���������� ���������� >> � <<. (Yes)
��� ��������� ������ ����������� � ������� �������������� ������. (Yes)
������������ �������� ���������� ��� ��������� ��������� �������� ��� ������ � �������. (Yes)
*/

#include <iostream>
#include "Debug.h"
#include "Matrix.h"
#include "MatrixUtility.h"

const char* const outline = "=====================================================================";

void subcounter(const char* str, bool reset = false) {
	static char c = 'a';
	if (reset) { c = 'a'; return; }
	std::cout << '\t' << c << ')' << str << std::endl;
	c++;
}

void counter(const char* str, bool reset = false) {
	static int c = 1;
	if (reset) { c = 1; return; }
	std::cout << c << ") " << str << std::endl;
	subcounter("", true);
	c++;
}

void begin_section(const char* title) {
	std::cout << outline << std::endl;
	std::cout << '=' << std::setw(strlen(outline) - 2)<<std::left << title << '=' << std::endl;
	std::cout << outline << std::endl << std::endl;
	counter("", true);
}

void end_section() {
	std::cout << outline << std::endl << std::endl;
}


 

int main() {
	setlocale(LC_ALL, "rus");
	//=====================================================
	// ���� ������������� �� ���������
	//=====================================================
	begin_section("������������ ��������� �������� �� ���������");
	{
		counter("����������� ������� �� ���������: Matrix<int> m");
		Matrix<int> m; //������� �� int ������� 1x1, ����������� ������
		std::cout << m;

		counter("����������� �������-������ �� ���������: RowVector<float> rv");
		RowVector<float> rv; //������-������ ������� 1, ����������� ������
		std::cout << rv;

		counter("����������� �������-������� �� ���������: RowVector<float> cv");
		RowVector<float> cv; //������-������ ������� 1, ����������� ������
		std::cout << cv;
	}
	end_section();
	//=====================================================

	//=====================================================
	// ��������������� ��������� �������� �� �� �������.
	//=====================================================
	begin_section("������������ ��������� �������� � ����� �����������");
	{
		counter("����������� ������� � ����� �����������: Matrix<int> m(3,3)");
		Matrix<int> m1(3, 3); //������� �� int ������� 3x3, ����������� ������
		std::cout << m1;

		counter("����������� ������� � ����� ����������: Matrix<int> m(3) = Matrix<int> m(3,1)");
		Matrix<int> m2(3); //������� �� int ������� 3x1, ����������� ������
		std::cout << m2;

		counter("����������� �������-������� � ����� ����������: ColVector<float> cv(3)");
		ColVector<float> cv(3); //�������-�������, ������� 3, ����������� ������
		std::cout << cv;

		counter("����������� �������-������ � ����� ����������: RowVector<float> rv(3)");
		RowVector<float> rv(3); //�������-������, ������� 3, ����������� ������
		std::cout << rv;
	}
	end_section();
	//=====================================================

	//=====================================================
	// ������������ �����������
	//=====================================================
	begin_section("������������ ����������� ��������� ��������");
	{
		counter("����������� ����������� �� ��������� ��� ������: Matrix<double> a = Matrix<double>{3,3}");
		Matrix<double> a = Matrix<int>{ 3,3 };
		std::cout << a;

		counter("����������� ����������� ������� � ������ ����� ���������: Matrix<int> b = Matrix<double>{ 2,3 }");
		Matrix<int> b = Matrix<float>{ 2,3 };
		std::cout << b;

		counter("����������� ����������� �� ��������� ��� ��������: RowVector<int> rv = RowVector<float>{ 4 }");
		RowVector<int> rv = RowVector<double>{ 4 };
		std::cout << rv;

		//RowVector<int> rv = ColVector<float>{ 4 }; //������!
		//std::cout << rv;

	}
	end_section();
	//=====================================================

	//=====================================================
	// ������������ c C-�������� � �������� ���������
	//=====================================================
	begin_section("������������ ��������� ��������� c C-�������� � �������� ���������");
	{
		int m_[] = { 1,2,3,4 };

		counter("�������� ������� �� ��������� �������: Matrix<float> m(m_, 2, 2), ��� int m_[] = { 1,2,3,4 };");
		Matrix<float> m(m_, 2, 2);
		std::cout << m;

		counter("�������� �������� �� ��������� �������: RowVector<int> v(m_, 3);");
		RowVector<int> v(m_, 3);
		std::cout << v;

	}
	end_section();
	//=====================================================

	//=====================================================
	// ������������ c �������� � �������� ���������
	//=====================================================
	begin_section("������������ ��������� ��������� c C-�������� � �������� ���������");
	{

		counter("�������� ������� �� ��������� ������� �� ��������: Matrix<float> m(m_), ��� std::vector<std::vector<int>> m_ = { {1,2},{2,3} };");
		std::vector<std::vector<int>> m_ = { {1,2, 0},{2,3, 1} };
		Matrix<float> m(m_);
		std::cout << m;

		subcounter("����������� ����������, ���� ������ �� ��������� �������: Matrix<float> m(incorrect_m_), ��� std::vector<std::vector<int>> incorrect_m_ = { {1,2},{2,3,1} }.");
		std::vector<std::vector<int>> incorrect_m_ = { {1,2},{2,3,1} };
		try {
			Matrix<float> m(incorrect_m_);
		}
		catch (std::exception& e) {
			std::cout << e.what() << std::endl;
		}

		counter("�������� ������� �� ����������� �������: RowVector<double> v(v_), ��� std::vector<float> v_ = { 1.0, 0.2, 3 }.");
		std::vector<float> v_ = { 1.0, 0.2, 3 };
		RowVector<double> v(v_);
		std::cout << v;
	}
	end_section();
	//=====================================================

	//=====================================================
	// ������������ c� ������� �������������.
	//=====================================================
	begin_section("������������ ��������� ��������� c� ������� �������������.");
	{
		counter("�������� ������� �� ���������� ������ �������������.");
		Matrix<float> m = {
			{1.2, 2.2, 3.2},
			{1.2, 0.3, 0.5},
			{9.8, 0.3, 0.9}
		};
		std::cout << m;

		counter("�������� ������� �� ����������� ������ �������������.");
		RowVector<int> v = { 1,2,3 };
		std::cout << v;

	}
	end_section();
	//=====================================================

	//=====================================================
	// ������� ���������
	//=====================================================
	begin_section("������� ���������.");
	{
		counter("������� ������: .Rows(), .Cols()");
		Matrix<float> m = {
			{1.2, 2.2, 3.2},
			{1.2, 0.3, 0.5},
			{9.8, 0.3, 0.9}
		};
		std::cout << m << std::endl;
		std::cout << "Rows: " << m.Rows() << ", Cols: " << m.Cols() << std::endl;

		counter("������� ��������: .Rows(), .Cols(), .Size()");
		RowVector<float> v = m[1];

		std::cout << v << std::endl;
		std::cout << "Rows: " << v.Rows() << ", Cols: " << v.Cols() << ", Size: " << v.Size() << std::endl;

		counter("������ � ��������� � ������ ������� At. ��������� ���������: m2.At(1, 1) = 0;");
		Matrix<float> m2 = {
			{1, 2, 3},
			{3, 5, 6},
			{7, 8, 9}
		};
		std::cout << m2 << std::endl;
		m2.At(1, 1) = 0;
		std::cout << m2 << std::endl;

		counter("������� At ��� ��������: v2.At(0, 3) = 0.");
		RowVector<int> v2 = { 1,2,3,4 };
		std::cout << v2 << std::endl;
		v2.At(0, 3) = 0;
		std::cout << v2 << std::endl;

		counter("operator[] ��� ������: m3[1]");
		Matrix<float> m3 = {
			{1, 2, 3},
			{3, 5, 6},
			{7, 8, 9}
		};
		std::cout << m3 << std::endl;
		std::cout << "������ ������: " << m3[1] << std::endl;


		counter("������� .GetRow(), .GetCol(): std::cout << m3.GetRow(1) << m3.GetCol(1);");
		std::cout << m3.GetRow(1) << m3.GetCol(1);

		counter("operator[] ��� ��������: v3[0] = 0.");
		RowVector<int> v3 = { 1,2,3,4 };
		std::cout << v3 << std::endl;
		v3[1] = 0;
		std::cout << v3 << std::endl;

	}
	end_section();
	//=====================================================

	//=====================================================
	// ������ ���������� � ����, �������� �� �����. ���������� � ���������� �� ���������� ������.
	//=====================================================
	begin_section("������ ���������� � ����, �������� �� �����. ���������� � ���������� ������ �� ������.");
	{
		counter("������� .Save � .ReadFromFile ��� ��������� � �������� ������");
		Matrix<float> m = {
			{1.2, 2.2, 3.2},
			{1.2, 0.3, 0.5},
			{9.8, 0.3, 0.9}
		};
		m.Save("materials/m_save.txt"); //��������� ������
		m.Save("materials/m_save.binm");//�������� ������

		Matrix<double> b;
		b.ReadFromFile("materials/m_save.txt");
		std::cout << "loaded from .txt: " << b << std::endl;
		b.ReadFromFile("materials/m_save.binm");
		std::cout << "loaded from .binm" << b << std::endl;

		counter("���������� � ���������� �� ������ �� ������.");
		std::fstream out("materials/stream.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
		if (!out.is_open()) return 1;

		Matrix<float> outM =
		{ {1,2,3},
			{3,3,4},
			{2,3,1}
		};

		Matrix<int> inM;

		out << outM;
		out.seekg(0);
		out >> inM;
		std::cout << inM;
		out.close();

		counter("����� ������� .Save() � ..ReadFromFile � �������������� ����������.");
		m.Save("materials/m_save.strage_ext1", TXT);
		m.Save("materials/m_save.strage_ext2", BINM);

		m.ReadFromFile("materials/m_save.strage_ext1", TXT);
		std::cout << m << std::endl;

		m.ReadFromFile("materials/m_save.strage_ext2", BINM);
		std::cout << m << std::endl;

	}
	end_section();
	//=====================================================

	//=====================================================
	// ����������� ���� ������.
	//=====================================================
	begin_section("����������� ���� ������");
	{
		counter("��������� ������� ������� n");
		std::cout << Identity<int>(10);

		counter("������-���������� �������.");
		std::vector<int> m1_ = 
		{1, 2, 3,
		    1, 2,
		       2};

		std::cout << UpperTriangular<int>(3, m1_);

		counter("�����-���������� �������.");
		std::vector<int> m2_ =
		{ 1,
		  1, 2,
		  2, 9 ,1 };

		std::cout << LowerTriangular<int>(3, m2_);

		counter("����������� �������");
		std::vector<int> m3_ =
		{
			    1,     2,   3,
			  /*2*/    2,   1,
			  /*3*/  /*1*/  0 
		};

		std::cout << Symmetric(3, m3_);
	}
	end_section();
	//=====================================================
	
	//=====================================================
	// ������ � ���������.
	//=====================================================
	begin_section("������ � ���������");
	{
		counter("��������� ������.");
		{
			Matrix<float> a = { {1.0,2.0,3.0}, {0.2, 1.3, 0.4}, {0.0 ,0.0, 1.0} };

			std::cout << "a = " << a;

			Matrix<double> b = { {4.0,2.0}, {2.2, 1.3}, {2.0 ,0.0} };

			std::cout << "b = " << b << std::endl;

			std::cout << "a * b = " << a * b << std::endl;

			ColVector<int> v = { 1,2,3 };

			std::cout << "c = " << v << std::endl;

			std::cout << "a * v = " << a * v << std::endl;
		}

		counter("��������, ��������� ������");
		{
			Matrix<int> A = { {1,2}, {3,4} };
			Identity<int> B(2);

			std::cout << "A = " << A <<std::endl;
			std::cout << "B = " << B << std::endl;

			std::cout << "A + B = " << A + B << std::endl;
			std::cout << "A - B = " << A - B << std::endl;

		}

		counter("��������, ��������� ��������");
		{
			ColVector<int> A = { 1,2,3 };

			std::cout << "A = " << A << std::endl;

			std::cout << "A + 0.5 * A = " << A + 0.5 * A << std::endl;
			std::cout << "A - 0.3 * A = " << A - 0.3 * A << std::endl;
		}

		counter("Create, ToMatrix, ScalarProduct, Transpose, HadamardProduct");
		{
			auto v = Create<ColVector<float>>(3, 1);
			auto m = Create<Matrix<float>>(3, 3);

			m = { {1,2,3}, {0, 3,4}, {2,3,4} };
			v = { 1,2,3 };
			

			std::cout << "m = " << m << std::endl;
			std::cout << "v = " << v << std::endl;

			std::cout << "ToMatrix(m) * ToMatrix(v) = " << ToMatrix(m) * ToMatrix(v) << std::endl;
			

			std::cout << "ScalarProduct(v, v) = "<< ScalarProduct(v, v) << std::endl;

			std::cout << "Transpose(m) = " << Transpose(m) << std::endl;

			std::cout << "HadamardProduct(m, m) = " << HadamardProduct(m, m) << std::endl;
		}

		counter("RowReduce, Rg, Det");
		{
			Matrix<double> m =
			{ {1., 0., 0.},
			  {4., 0., 0.},
			  {1., 0., 0.} 
			};
		

			double det;
			size_t rank;
			std::cout << RowReduce(m, &det, &rank) << std::endl;
			std::cout << "Det = " << det << std::endl;
			std::cout << "Rg = " << rank << std::endl;

			std::cout << "Det = " << Det(m) << std::endl;
			std::cout << "Rg = " << Rg(m) << std::endl;

		}

		counter("EuclidNorm, MaxNorm, FrobeniusNorm");
		{
			RowVector<int> v = { 1,2,3 };
			std::cout << "v = " << v << std::endl;

			std::cout << "EuclidNorm = " << EuclidNorm(v) << std::endl;

			std::cout << "MaxNorm = " << MaxNorm(v) << std::endl;

			std::cout << "v^T * v = " << Transpose(v) * v << std::endl;

			std::cout << "FrobeniusNorm(v^T * v) = " << FrobeniusNorm(Transpose(v) *  v) << std::endl;
		}
		
		counter("Angle");
		{
			RowVector<int> a = { 1,0,1 };
			ColVector<int> b = { 0,0,1 };

			std::cout << a << std::endl;
			std::cout << b << std::endl;
			std::cout << "deg = " << Angle(a, b, DEG) << std::endl;
			std::cout << "rad = " << Angle(a, b, RAD) << std::endl;
		}

		counter("SubMatrix, Assign, Inverse");
		{
			Matrix<float> m =
			{
				{1, 0},
				{0, 1},
				{4, 1}
			};

			ColVector<float> v = { 1, 2, 3 };
			auto a = Assign(m, v);
			std::cout << a;

			auto b = Assign(SubMatrix(a, EntryPosition{ 1, 0 }, EntryPosition{ 2, 1 }), ColVector<float>{ 1,2 });
			std::cout << b;

			Matrix<float> inv = Inverse(a); 
			std::cout <<"m^(-1) = " << inv << std::endl;


			std::cout <<"m * m^(-1) = " << a * inv << std::endl;
		}
	}
	end_section();


	begin_section("�������������� ����.");
	{
		Matrix<double> m =
		{ {1., 0., 0.},
		  {4., 5., 6.},
		};

		std::cout << Transpose(m);
		

		RowVector<int> rv = { 1,2,3 };
		std::cout << Transpose(rv);

	}
	end_section();


	//=====================================================
}
