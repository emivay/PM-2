#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
	bool built;
	unsigned int rows;
	unsigned int cols;
	std::vector<std::vector<double>> matrix;
public:
	//Constructors
	Matrix();
	Matrix(unsigned int, unsigned int);
	Matrix(unsigned int, unsigned int, double);
	Matrix(unsigned int, unsigned int, std::vector<std::vector<double>>);
	Matrix(std::vector<std::vector<double>>);
	//Getters
	const bool isSingular(void) const;
	const bool isBuilt(void) const;
	const unsigned int getRows(void) const;
	const unsigned int getCols(void) const;
	const std::vector<std::vector<double>> getMatrix(void) const;
	const std::vector<double> getRow(unsigned int) const;
	const Matrix getMatrixRow(unsigned int) const;
	const std::vector<double> getCol(unsigned int) const;
	const Matrix getMatrixCol(unsigned int) const;
	const double at(unsigned int, unsigned int) const;
	//Modifiers
	bool eraseRow(unsigned int);
	bool eraseRows(unsigned int, unsigned int);
	bool eraseCol(unsigned int);
	bool eraseCols(unsigned int, unsigned int);
	bool insertRow(std::vector<double>, unsigned int);
	bool replaceRow(std::vector<double>, unsigned int);
	bool replaceMatrixRow(Matrix, unsigned int);
	bool insertCol(std::vector<double>, unsigned int);
	bool replaceCol(std::vector<double>, unsigned int);
	bool replaceMatrixCol(Matrix, unsigned int);
	bool insertMatrixInRow(Matrix m, unsigned int);
	bool insertMatrixInCol(Matrix m, unsigned int); 
	void roundMatrix(double);
	void set(unsigned int, unsigned int, double);
	bool resize(unsigned int, unsigned int);
	//Read and print
	void read(void);
	void print(void) const;
	void info(void) const;
	//Operations
	void scalarProduct(double);
	double determinant(void) const;
	static Matrix identity(unsigned int);
	static Matrix random(unsigned int, unsigned int, int, int);
	Matrix subMatrix(unsigned int, unsigned int, unsigned int, unsigned int) const;
	static Matrix scalarProduct(Matrix, double);
	static Matrix matrixSum(Matrix, Matrix);
	static Matrix matrixProduct(Matrix, Matrix);
	static Matrix transpose(Matrix);
	static Matrix adjoint(Matrix);
	static Matrix inverse(Matrix);
};

#endif //MATRIX_H
