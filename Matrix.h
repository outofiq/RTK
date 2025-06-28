#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

namespace Matrix_hx
{
	class Matrix
	{
	private:
		vector<vector<double>> data;
	public:
		//���캯��
		Matrix(size_t rows, size_t cols, double value);
		//��ȡ���������
		size_t rows() const;
		//��ȡ���������
		size_t cols() const;
		//��ӡ����
		void printMatrix() const;
		//����()����������ڷ��ʾ���Ԫ��
		double& operator()(size_t row, size_t col);
		//����()����������ڷ��ʾ���Ԫ��(const�汾)
		const double& operator()(size_t row, size_t col) const;
		//����+�������ʵ�־���ӷ�
		Matrix operator+(const Matrix& other) const;
		//����-�������ʵ�־������
		Matrix operator-(const Matrix& other) const;
		//����*�������ʵ�־���˷�
		Matrix operator*(const Matrix& other) const;
		//ʵ�־����ת��
		Matrix transpose() const;
		//��������(���ø�˹��Ԫ���������������С���һ�г���һ����������һ�м�����һ�еı���)
		//��������
		void swapRows(size_t row1, size_t row2);
		//��һ�г���һ������
		void multiplyRow(size_t row, double scalar);
		//��һ�м�����һ�еı���
		void addRowToRow(size_t fromRow, size_t toRow, double scalar);
		//��������
		Matrix inverse();
		//ɾ����
		void deletRow(int rowIndex);
		//ɾ����
		void deleteColumn(int colIndex);
		//����^�������ʵ����ά�������
		Matrix operator^(const Matrix& other) const;
	};
}
