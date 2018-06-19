#include<math.h>
#include <iostream>
using namespace std;
#ifndef _DATASTRUCT_LMATRIX_H_
#define _DATASTRUCT_LMATRIX_H_

#ifndef LTEMPLATE
#define LTEMPLATE template<typename Type>
#endif

#ifndef IN
#define IN
#endif

#ifndef INOUT
#define INOUT
#endif

#ifndef OUT
#define OUT
#endif


LTEMPLATE
class LMatrix
{
public:
	
	static bool MUL(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C);

	
	static bool DOTMUL(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C);

	
	static LMatrix<Type> DOTMUL(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B);

	static bool SCALARMUL(IN const LMatrix<Type>& A, IN const Type& B, OUT LMatrix<Type>& C);

	
	static bool DOTDIV(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C);

	
	static bool SUB(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C);

	
	static bool ADD(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C);

	
	static bool T(IN const LMatrix<Type>& A, OUT LMatrix<Type>& B);
	
	static LMatrix<Type> INV(IN const LMatrix<Type>& A);


public:
	
	LMatrix();

	
	~LMatrix();

	
	LMatrix(IN unsigned int row, IN unsigned int col);


	LMatrix(IN unsigned int row, IN unsigned int col, IN const Type& initValue);

	
	LMatrix(IN unsigned int row, IN unsigned int col, IN const Type* pDataList);

	
	LMatrix(IN const LMatrix<Type>& rhs);


	LMatrix<Type>& operator = (IN const LMatrix<Type>& rhs);

	
	LMatrix<Type> operator * (IN const LMatrix<Type>& B) const;

	
	LMatrix<Type> operator - (IN const LMatrix<Type>& B) const;

	
	LMatrix<Type> operator + (IN const LMatrix<Type>& B) const;

	
	Type*& operator[](IN unsigned int row);

	
	const Type* operator[](IN unsigned int row) const;

	LMatrix<Type> T() const;


	LMatrix<Type> ScalarMul(IN const Type& B) const;

	
	LMatrix<Type> GetRow(IN unsigned int row) const;


	void GetRow(IN unsigned int row, OUT LMatrix<Type>& rowVector) const;


	LMatrix<Type> GetColumn(IN unsigned int col) const;

	
	void GetColumn(IN unsigned int col, OUT LMatrix<Type>& colVector) const;


	void Reset(IN unsigned int row, IN unsigned int col);


	void Reset(IN unsigned int row, IN unsigned int col, IN const Type& initValue);

	LMatrix<Type> INV();

	void eye(int n);

	void show() const
	{
		for (unsigned int i = 0; i < RowLen; i++)
		{
			cout << endl;
			for (unsigned int j = 0; j < ColumnLen; j++)
				cout << Data[i][j] << " ";
		}
	}

public:
	unsigned int RowLen; ///< 行长度
	unsigned int ColumnLen; ///< 列长度
	Type** Data; ///< 数据列表
private:
	Type * m_data; ///< 实际存储的数据
};

LTEMPLATE
LMatrix<Type>::LMatrix()
	: RowLen(0), ColumnLen(0), Data(0), m_data(0)
{

}

LTEMPLATE
LMatrix<Type>::LMatrix(IN unsigned int row, IN unsigned int col)
	: RowLen(0), ColumnLen(0), Data(0), m_data(0)
{
	this->Reset(row, col);
}

LTEMPLATE
LMatrix<Type>::LMatrix(IN unsigned int row, IN unsigned int col, IN const Type& initValue)
	: RowLen(0), ColumnLen(0), Data(0), m_data(0)
{
	this->Reset(row, col);

	unsigned int size = row * col;
	for (unsigned int i = 0; i < size; i++)
	{
		this->m_data[i] = initValue;
	}
}

LTEMPLATE
LMatrix<Type>::LMatrix(IN unsigned int row, IN unsigned int col, IN const Type* pDataList)
	: RowLen(0), ColumnLen(0), Data(0), m_data(0)
{
	this->Reset(row, col);

	unsigned int size = row * col;
	for (unsigned int i = 0; i < size; i++)
	{
		this->m_data[i] = pDataList[i];
	}

}

LTEMPLATE
LMatrix<Type>::~LMatrix()
{
	if (this->Data)
	{
		delete[] this->Data;
		this->Data = 0;
	}

	if (this->m_data)
	{
		delete[] this->m_data;
		this->m_data = 0;
	}

	this->RowLen = 0;
	this->ColumnLen = 0;
}

LTEMPLATE
LMatrix<Type>::LMatrix(IN const LMatrix<Type>& rhs)
	: RowLen(0), ColumnLen(0), Data(0), m_data(0)
{
	this->Reset(rhs.RowLen, rhs.ColumnLen);

	unsigned int size = rhs.RowLen * rhs.ColumnLen;
	for (unsigned int i = 0; i < size; i++)
	{
		this->m_data[i] = rhs.m_data[i];
	}

}

LTEMPLATE
LMatrix<Type>& LMatrix<Type>::operator = (IN const LMatrix<Type>& rhs)
{
	this->Reset(rhs.RowLen, rhs.ColumnLen);

	unsigned int size = rhs.RowLen * rhs.ColumnLen;
	for (unsigned int i = 0; i < size; i++)
	{
		this->m_data[i] = rhs.m_data[i];
	}

	return *this;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::operator * (IN const LMatrix<Type>& B) const
{
	LMatrix<Type> C;
	MUL(*this, B, C);
	return C;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::operator - (IN const LMatrix<Type>& B) const
{
	LMatrix<Type> C;
	SUB(*this, B, C);
	return C;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::operator + (IN const LMatrix<Type>& B) const
{
	LMatrix<Type> C;
	ADD(*this, B, C);
	return C;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::ScalarMul(IN const Type& B) const
{
	LMatrix<Type> C;
	SCALARMUL(*this, B, C);
	return C;
}

LTEMPLATE
Type*& LMatrix<Type>::operator[](IN unsigned int row)
{
	return this->Data[row];
}

LTEMPLATE
const Type* LMatrix<Type>::operator[](IN unsigned int row) const
{
	return this->Data[row];
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::T() const
{
	LMatrix<Type> B;
	T(*this, B);
	return B;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::GetRow(IN unsigned int row) const
{
	LMatrix<Type> rowVector(1, this->ColumnLen);
	this->GetRow(row, rowVector);

	return rowVector;
}

LTEMPLATE
void LMatrix<Type>::GetRow(IN unsigned int row, OUT LMatrix<Type>& rowVector) const
{
	rowVector.Reset(1, this->ColumnLen);
	for (unsigned int i = 0; i < this->ColumnLen; i++)
	{
		rowVector.Data[0][i] = this->Data[row][i];
	}
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::GetColumn(IN unsigned int col) const
{
	LMatrix<Type> columnVector(this->RowLen, 1);
	this->GetColumn(col, columnVector);

	return columnVector;
}

LTEMPLATE
void LMatrix<Type>::GetColumn(IN unsigned int col, OUT LMatrix<Type>& colVector) const
{
	colVector.Reset(this->RowLen, 1);
	for (unsigned int i = 0; i < this->RowLen; i++)
	{
		colVector.Data[i][0] = this->Data[i][col];
	}
}

LTEMPLATE
void LMatrix<Type>::Reset(IN unsigned int row, IN unsigned int col)
{
	if ((this->RowLen != row) || this->ColumnLen != col)
	{
		if (this->Data)
		{
			delete[] this->Data;
			this->Data = 0;
		}

		if (this->m_data)
		{
			delete[] this->m_data;
			this->m_data = 0;
		}


		if (row * col > 0)
		{
			this->RowLen = row;
			this->ColumnLen = col;

			this->Data = new Type*[this->RowLen];
			this->m_data = new Type[this->RowLen * this->ColumnLen];
			for (unsigned int i = 0; i < this->RowLen; i++)
			{
				this->Data[i] = &this->m_data[this->ColumnLen * i];
			}
		}
		else
		{
			this->RowLen = 0;
			this->ColumnLen = 0;
		}
	}
}

LTEMPLATE
void LMatrix<Type>::Reset(IN unsigned int row, IN unsigned int col, IN const Type& initValue)
{
	this->Reset(row, col);

	unsigned int size = row * col;
	for (unsigned int i = 0; i < size; i++)
	{
		this->m_data[i] = initValue;
	}
}

LTEMPLATE
bool LMatrix<Type>::MUL(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C)
{
	if (A.ColumnLen != B.RowLen)
		return false;

	C.Reset(A.RowLen, B.ColumnLen);

	for (unsigned int i = 0; i < C.RowLen; i++)
	{
		for (unsigned int j = 0; j < C.ColumnLen; j++)
		{
			C.Data[i][j] = A.Data[i][0] * B.Data[0][j];
			for (unsigned int k = 1; k < A.ColumnLen; k++)
			{
				C.Data[i][j] += A.Data[i][k] * B.Data[k][j];
			}
		}
	}

	return true;
}

LTEMPLATE
bool LMatrix<Type>::DOTMUL(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C)
{
	if ((A.RowLen != B.RowLen) || (A.ColumnLen != B.ColumnLen))
		return false;

	C.Reset(A.RowLen, A.ColumnLen);

	for (unsigned int i = 0; i < C.RowLen; i++)
	{
		for (unsigned int j = 0; j < C.ColumnLen; j++)
		{
			C.Data[i][j] = A.Data[i][j] * B.Data[i][j];
		}
	}

	return true;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::DOTMUL(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B)
{
	LMatrix<Type> C;
	DOTMUL(A, B, C);
	return C;
}

LTEMPLATE
bool LMatrix<Type>::SCALARMUL(IN const LMatrix<Type>& A, IN const Type& B, OUT LMatrix<Type>& C)
{
	C.Reset(A.RowLen, A.ColumnLen);
	for (unsigned int row = 0; row < A.RowLen; row++)
	{
		for (unsigned int col = 0; col < A.ColumnLen; col++)
		{
			C.Data[row][col] = A.Data[row][col] * B;
		}
	}

	return true;
}

LTEMPLATE
bool LMatrix<Type>::SUB(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C)
{
	if ((A.RowLen != B.RowLen) || (A.ColumnLen != B.ColumnLen))
		return false;

	C.Reset(A.RowLen, A.ColumnLen);

	for (unsigned int i = 0; i < C.RowLen; i++)
	{
		for (unsigned int j = 0; j < C.ColumnLen; j++)
		{
			C.Data[i][j] = A.Data[i][j] - B.Data[i][j];
		}
	}

	return true;
}

LTEMPLATE
bool LMatrix<Type>::ADD(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C)
{
	if ((A.RowLen != B.RowLen) || (A.ColumnLen != B.ColumnLen))
		return false;

	C.Reset(A.RowLen, A.ColumnLen);

	for (unsigned int i = 0; i < C.RowLen; i++)
	{
		for (unsigned int j = 0; j < C.ColumnLen; j++)
		{
			C.Data[i][j] = A.Data[i][j] + B.Data[i][j];
		}
	}

	return true;
}

LTEMPLATE
bool LMatrix<Type>::DOTDIV(IN const LMatrix<Type>& A, IN const LMatrix<Type>& B, OUT LMatrix<Type>& C)
{
	if ((A.RowLen != B.RowLen) || (A.ColumnLen != B.ColumnLen))
		return false;

	C.Reset(A.RowLen, A.ColumnLen);

	for (unsigned int i = 0; i < C.RowLen; i++)
	{
		for (unsigned int j = 0; j < C.ColumnLen; j++)
		{
			C.Data[i][j] = A.Data[i][j] / B.Data[i][j];
		}
	}

	return true;
}

LTEMPLATE
bool LMatrix<Type>::T(IN const LMatrix<Type>& A, OUT LMatrix<Type>& B)
{
	B.Reset(A.ColumnLen, A.RowLen);
	for (unsigned int i = 0; i < A.RowLen; i++)
	{
		for (unsigned int j = 0; j < A.ColumnLen; j++)
			B.Data[j][i] = A.Data[i][j];
	}

	return true;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::INV(IN const LMatrix<Type>& A)
{
	LMatrix<float> TEMP(A.RowLen, A.ColumnLen * 2);
	unsigned int i, j, k;
	for (i = 0; i < A.ColumnLen; i++)
	{
		for (j = 0; j < 2 * A.ColumnLen; j++)
		{
			if (j < A.ColumnLen)
			{
				TEMP[i][j] = A[i][j];
			}
			else if (j == i + A.ColumnLen)
				TEMP[i][j] = 1.0;
			else
				TEMP[i][j] = 0.0;
		}
	}
	int maxI = 0;
	for (i = 1; i < A.ColumnLen; i++)
	{
		if (fabs(TEMP[maxI][0]) < fabs(TEMP[i][0]))
			maxI = i;
	}
	if (maxI != 0)
	{
		double temp;
		for (j = 0; j < 2 * A.ColumnLen; j++)
		{
			temp = TEMP[0][j];
			TEMP[0][j] = TEMP[maxI][j];
			TEMP[maxI][j] = temp;
		}
	}
	double temp2;
	for (i = 0; i < A.ColumnLen; i++)
	{
		if (TEMP[i][i] != 0)
			temp2 = 1.0 / TEMP[i][i];
		else
		{
			cout << "此矩阵无逆!" << endl;
			return A;
		}
		for (j = 0; j < 2 * A.ColumnLen; j++)
			TEMP[i][j] *= temp2;
		for (j = 0; j < A.ColumnLen; j++)
		{
			if (j != i)
			{
				double temp3 = TEMP[j][i];
				for (k = 0; k < 2 * A.ColumnLen; k++)
					TEMP[j][k] -= temp3 * TEMP[i][k];
			}
		}
	}
	LMatrix<float> output(A.RowLen, A.ColumnLen);
	for (i = 0; i < A.RowLen; i++)
	{
		for (j = 0; j < A.ColumnLen; j++)
		{
			output[i][j] = TEMP[i][j + A.ColumnLen];
		}
	}
	return output;
}

LTEMPLATE
LMatrix<Type> LMatrix<Type>::INV()
{
	LMatrix<float> A = *this;
	LMatrix<float> TEMP(A.RowLen, A.ColumnLen * 2);
	unsigned int i, j, k;
	for (i = 0; i < A.ColumnLen; i++)
	{
		for (j = 0; j < 2 * A.ColumnLen; j++)
		{
			if (j < A.ColumnLen)
			{
				TEMP[i][j] = A[i][j];
			}
			else if (j == i + A.ColumnLen)
				TEMP[i][j] = 1.0;
			else
				TEMP[i][j] = 0.0;
		}
	}
	int maxI = 0;
	for (i = 1; i < A.ColumnLen; i++)
	{
		if (fabs(TEMP[maxI][0]) < fabs(TEMP[i][0]))
			maxI = i;
	}
	if (maxI != 0)
	{
		double temp;
		for (j = 0; j < 2 * A.ColumnLen; j++)
		{
			temp = TEMP[0][j];
			TEMP[0][j] = TEMP[maxI][j];
			TEMP[maxI][j] = temp;
		}
	}
	double temp2;
	for (i = 0; i < A.ColumnLen; i++)
	{
		if (TEMP[i][i] != 0)
			temp2 = 1.0 / TEMP[i][i];
		else
		{
			cout << "此矩阵无逆!" << endl;
			return A;
		}
		for (j = 0; j < 2 * A.ColumnLen; j++)
			TEMP[i][j] *= temp2;
		for (j = 0; j < A.ColumnLen; j++)
		{
			if (j != i)
			{
				double temp3 = TEMP[j][i];
				for (k = 0; k < 2 * A.ColumnLen; k++)
					TEMP[j][k] -= temp3 * TEMP[i][k];
			}
		}
	}
	LMatrix<float> output(A.RowLen, A.ColumnLen);
	for (i = 0; i < A.RowLen; i++)
	{
		for (j = 0; j < A.ColumnLen; j++)
		{
			output[i][j] = TEMP[i][j + A.ColumnLen];
		}
	}
	return output;
}

LTEMPLATE
void LMatrix<Type>::eye(int n)
{
	this->Reset(n, n);
	for (unsigned int i = 0; i < n; i++)
	{
		for (unsigned int j = 0; j < n; j++)
		{
			if (i != j)
				this->Data[i][j] = 0;
			else
				this->Data[i][j] = 1;
		}
	}
}


#endif
