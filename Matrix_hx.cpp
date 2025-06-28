#include "Matrix.h"

using namespace Matrix_hx;

/*
**********************************************************
�������������ʼ������
������rows    ��������
      cols	  ��������
      value   ����ÿ��Ԫ�صĳ�ʼֵ
�������ܣ����첢��ʼ������
**********************************************************
*/
Matrix::Matrix(size_t rows, size_t cols, double value = 0.0)
{
    data.resize(rows, vector<double>(cols, value));
}

/*
**********************************************************
����������������
����ֵ�����������
�������ܣ���þ��������
**********************************************************
*/
size_t Matrix::rows() const
{
    return data.size();
}

/*
**********************************************************
����������������
����ֵ�����������
�������ܣ���þ��������
**********************************************************
*/
size_t Matrix::cols() const
{
    return data.empty() ? 0 : data[0].size();
}

/*
**********************************************************
����������ӡ����
�������ܣ��������ӡ�������Ļ
**********************************************************
*/
void Matrix::printMatrix() const
{
    for (const auto& row : data)
    {
        for (const auto& elem : row)
        {
            cout << elem << " ";
        }
        cout << std::endl;
    }
}

/*
**********************************************************
������������()
������row    ����ĳһ�е���������1�е�������0��
      col	 ����ĳһ�е���������1�е�������0��
����ֵ�������row+1�е�col+1�е�Ԫ��ֵ
�������ܣ����������()��ͨ������������ȡ����ĳһԪ�ص�ֵ���ҿɸı�
**********************************************************
*/
double& Matrix::operator()(size_t row, size_t col)
{
    if (row >= rows() || col >= cols())
    {
        throw out_of_range("Index out of range");
    }
    return data[row][col];
}

/*
**********************************************************
������������()
������row    ����ĳһ�е���������1�е�������0��
      col	 ����ĳһ�е���������1�е�������0��
����ֵ�������row+1�е�col+1�е�Ԫ��ֵ
�������ܣ����������()��ͨ������������ȡ����ĳһԪ�ص�ֵ�������ܸı�
**********************************************************
*/
const double& Matrix::operator()(size_t row, size_t col) const
{
    if (row >= rows() || col >= cols())
    {
        throw out_of_range("Index out of range");
    }
    return data[row][col];
}

/*
**********************************************************
������������+
������other    ��һ������
����ֵ������������ӵõ��ľ���
�������ܣ����������+��ʵ�־�������
**********************************************************
*/
Matrix Matrix::operator+(const Matrix& other) const
{
    if (rows() != other.rows() || cols() != other.cols())
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(rows(), cols(), 0);
    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }
    return result;
}

/*
**********************************************************
������������-
������other    ��һ������(����)
����ֵ��������������õ��ľ���
�������ܣ����������-��ʵ�־�������
**********************************************************
*/
Matrix Matrix::operator-(const Matrix& other) const
{
    if (rows() != other.rows() || cols() != other.cols())
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(rows(), cols(), 0);
    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            result(i, j) = (*this)(i, j) - other(i, j);
        }
    }
    return result;
}

/*
**********************************************************
������������*
������other    �˺��ұߵľ���
����ֵ������������˵õ��ľ���
�������ܣ����������*��ʵ�־�������
**********************************************************
*/
Matrix Matrix::operator*(const Matrix& other) const
{
    if (cols() != other.rows())
    {
        throw invalid_argument("Matrix A's columns must be equal to Matrix B's rows");
    }
    Matrix result(rows(), other.cols(), 0);
    for (size_t i = 0; i < result.rows(); ++i)
    {
        for (size_t j = 0; j < result.cols(); ++j)
        {
            result(i, j) = 0.0;
            for (size_t k = 0; k < cols(); ++k)
            {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}

/*
**********************************************************
��������ת��
����ֵ��ת�Ⱥ�ľ���
�������ܣ�ʵ�־����ת��
**********************************************************
*/
Matrix Matrix::transpose() const
{
    Matrix result(cols(), rows(), 0);
    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            result(j, i) = (*this)(i, j);
        }
    }
    return result;
}

/*
**********************************************************
��������������
������row1    ����ĳһ�е�����
      row2    ������һ�е�����
�������ܣ����������ĳ����
**********************************************************
*/
void Matrix::swapRows(size_t row1, size_t row2)
{
    swap_ranges(data[row1].begin(), data[row1].end(), data[row2].begin());
}

/*
**********************************************************
�����������г���һ������
������row     ����ĳһ�е�����
      scalar  ���˵ı���
�������ܣ��������ĳһ�г���һ������
**********************************************************
*/
void Matrix::multiplyRow(size_t row, double scalar)
{
    for (auto& elem : data[row])
    {
        elem *= scalar;
    }
}

/*
**********************************************************
��������һ�м�����һ�еı���
������fromrow    ����
      torow      ���ı����
      scalar     ���� 
�������ܣ���ĳһ�г���һ�������ӵ���һ���ϣ��ı���һ�е�ֵ
**********************************************************
*/
void Matrix::addRowToRow(size_t fromRow, size_t toRow, double scalar)
{
    for (size_t col = 0; col < data[0].size(); ++col)
    {
        data[toRow][col] += data[fromRow][col] * scalar;
    }
}

/*
**********************************************************
����������������
����ֵ������������
�������ܣ�ʹ�ø�˹��Ԫ������þ���������
**********************************************************
*/
Matrix Matrix::inverse()
{
    if (rows() != cols())
    {
        throw invalid_argument("Matrix must be square for inversion");
    }
    Matrix A = *this;//��������ĸ���������ı�ԭ����
    size_t n = A.rows();
    if (n <= 0)
    {
        throw invalid_argument("Matrix's rows or cols must more than 0 ");
    }
    Matrix inverse(n, n, 0.0);

    // ������λ����
    for (size_t i = 0; i < n; ++i)
    {
        inverse(i, i) = 1.0;
    }

    // ��˹-Լ����Ԫ��
    for (size_t i = 0; i < n; ++i)
    {
        // Ѱ��������Ԫ��
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k)
        {
            if (abs(A(k, i)) > abs(A(maxRow, i)))
            {
                maxRow = k;
            }
        }
        if (A(maxRow, i) == 0)
        {
            throw runtime_error("Matrix is singular and cannot be inverted");
        }
        if (maxRow != i)
        {
            A.swapRows(i, maxRow);
            inverse.swapRows(i, maxRow);
        }
        //���Խ�Ԫ�ع淶��Ϊ1�������������
        double div = A(i, i);
        A.multiplyRow(i, 1.0 / div);
        inverse.multiplyRow(i, 1.0 / div);

        for (size_t j = 0; j < n; ++j)
        {
            if (i != j)
            {
                double scalar = -A(j, i);
                A.addRowToRow(i, j, scalar);
                inverse.addRowToRow(i, j, scalar);
            }
        }
    }

    return inverse;
}

/*
**********************************************************
��������ɾ����
������rowIndex    ����ĳһ�е�����
�������ܣ��������ĳһ��ɾ�����Ӷ����پ��������
**********************************************************
*/
void Matrix::deletRow(int rowIndex)
{
    if (rowIndex >= 0 && rowIndex < data.size()) 
    {
        data.erase(data.begin() + rowIndex);
    }
}

/*
**********************************************************
��������ɾ����
������colIndex    ����ĳһ�е�����
�������ܣ��������ĳһ��ɾ�����Ӷ����پ��������
**********************************************************
*/
void Matrix::deleteColumn(int colIndex) 
{
    if (colIndex >= 0 && colIndex < data[0].size()) 
    {
        for (auto& row : data) 
        {
            row.erase(row.begin() + colIndex);
        }
    }
}

/*
**********************************************************
������������~
������other    ��һ����ά����
����ֵ��������ά������˵õ�����ά����
�������ܣ����������~��ʵ����ά�����Ĳ��
**********************************************************
*/
Matrix Matrix::operator^(const Matrix& other) const
{
    if (rows() != 1 || cols() != 3 || other.rows() != 1 || other.cols() != 3)
    {
        throw invalid_argument("only tri-vector can calculate cross multiplication.");
    }
    Matrix result(1, 3, 0);
    result(0, 0) = (*this)(0, 1) * other(0, 2) - (*this)(0, 2) * other(0, 1);
    result(0, 1) = (*this)(0, 2) * other(0, 0) - (*this)(0, 0) * other(0, 2);
    result(0, 2) = (*this)(0, 0) * other(0, 1) - (*this)(0, 1) * other(0, 0);

    return result;
}
