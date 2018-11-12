#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include "Matrix.h"
using namespace std;

//CONSTRUCTORS

Matrix::Matrix() : built(false), rows(0), cols(0) {
	
	matrix = std::vector<std::vector<double>>(0, std::vector<double>(0, 0));
	
}

Matrix::Matrix(unsigned int matrixRows, unsigned int matrixCols) : built(true), rows(matrixRows), cols(matrixCols){
	
	matrix = std::vector<std::vector<double>>(matrixRows, std::vector<double>(matrixCols, 0));
		
}

Matrix::Matrix(unsigned int matrixRows, unsigned int matrixCols, double value) : built(true), rows(matrixRows), cols(matrixCols) {
	
	matrix = std::vector<std::vector<double>>(matrixRows, std::vector<double>(matrixCols, value));
	
}

Matrix::Matrix(unsigned int matrixRows, unsigned int matrixCols, std::vector<std::vector<double>> newMatrix) : built(true), rows(matrixRows), cols(matrixCols) {
	
	matrix = std::vector<std::vector<double>>(matrixRows, std::vector<double>(matrixCols, 0));
	for(unsigned int i = 0; i < matrixRows; i++)
		for(unsigned int j = 0; j < matrixCols; j++)
			matrix[i][j] = newMatrix[i][j];
	
}

Matrix::Matrix(std::vector<std::vector<double>> newMatrix) : built(true), rows(newMatrix.size()), cols(newMatrix[0].size()) {
	
	matrix = std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0));

	for(unsigned int i = 0; i < rows; i++)
		for(unsigned int j = 0; j < cols; j++)
			matrix[i][j] = newMatrix[i][j];
		
}


//GETTERS

const bool Matrix::isSingular(void) const {
	return (determinant() != 0 ? true : false);
}

const bool Matrix::isBuilt(void) const {
	return built;
}

const unsigned int Matrix::getRows(void) const {
	return rows;
}

const unsigned int Matrix::getCols(void) const {
	return cols;
}

const std::vector<std::vector<double>> Matrix::getMatrix(void) const {
	return this->matrix;
}

const std::vector<double> Matrix::getRow(unsigned int row) const {

	if(row < this->getRows())
		return matrix[row];
	return std::vector<double>(this->getCols(), 0);
}

const Matrix Matrix::getMatrixRow(unsigned int row) const {
	
	Matrix matrixRow;
	
	if(row < 0 or row >= this->getRows())
		return matrixRow;
		
	int size = this->getCols();
	matrixRow = Matrix(1,size);
	for(int j = 0; j < size; j++)
		matrixRow.set(0,j,this->matrix[row][j]);
	
	return matrixRow;
}

const std::vector<double> Matrix::getCol(unsigned int col) const {
	std::vector<double> column(this->getRows(), 0);

	if(col > 0 and col < this->getCols()) {
		for(unsigned int i = 0; i < this->getRows(); i++)
			column[i] = matrix[i][col];
	}
	
	return column;
}

const Matrix Matrix::getMatrixCol(unsigned int col) const {
		
	Matrix matrixCol;
	
	if(col < 0 or col >= this->getCols())
		return matrixCol;
	
	int size = this->getRows();
	matrixCol = Matrix(size,1);
	for(int i = 0; i < size; i++)
		matrixCol.set(i,0,matrix[i][col]);
	
	return matrixCol;
		
}

const double Matrix::at(unsigned int i, unsigned int j) const {
	
	if(i < 0 or i >= this->getRows())
		return 0;
	
	if(j < 0 or j >= this->getCols())
		return 0;
	
	return matrix[i][j];
}


//MODIFIERS

bool Matrix::eraseRow(unsigned int row) {
	
	if(row < this->getRows()) {
		matrix.erase(matrix.begin() + row);
		rows = this->getRows() - 1;
		return true;
	}
	return false;
}

bool Matrix::eraseRows(unsigned int ini_row, unsigned int fin_row) {
	
	if(ini_row <= fin_row and fin_row < this->getRows()) {
		matrix.erase(matrix.begin() + ini_row, matrix.begin() + fin_row + 1);
		rows = this->getRows() - (fin_row + 1 - ini_row);
		return true;
	}
	return false;
	
}

bool Matrix::eraseCol(unsigned int col) {
		
	if(col < this->getCols()) {
		for(unsigned int i = 0; i < this->getRows(); i++)
			matrix[i].erase(matrix[i].begin() + col);
		cols = this->getCols() - 1;
		return true;
	}
	return false;
}

bool Matrix::eraseCols(unsigned int ini_col, unsigned int fin_col) {
	
	if(ini_col <= fin_col and fin_col < this->getCols()) {
		
		for(unsigned int i = 0; i < this->getRows(); i++)
			matrix[i].erase(matrix[i].begin() + ini_col, matrix[i].begin() + fin_col + 1);
		
		cols = this->getCols() - (fin_col + 1 - ini_col);
		return true;
	}
	return false;
	
} 

bool Matrix::insertRow(std::vector<double> newRow, unsigned int row) {
	
	if(newRow.size() != this->getCols())
		return false;
	
	if(row < 0 or row > this->getRows())
		return false;
	
	matrix.insert(matrix.begin() + row, newRow);
	rows++;
	return true;
}

bool Matrix::replaceRow(std::vector<double> newRow, unsigned int row) {
	
	if(row < 0 or row >= this->getRows())
		return false;
		
	if(newRow.size() != this->getCols())
		return false;
		
	for(unsigned int j = 0; j < this->getCols(); j++)
		matrix[row][j] = newRow[j];
		
	return true;
}

bool Matrix::replaceMatrixRow(Matrix matrixRow, unsigned int row) {
	
	if(row < 0 or row >= this->getRows())
		return false;
		
	if(matrixRow.getRows() != 1)
		return false;
		
	if(matrixRow.getCols() != this->getCols())
		return false;
		
	for(unsigned int j = 0; j < this->getCols(); j++)
		matrix[row][j] = matrixRow.at(0,j);
		
	return true;
}

bool Matrix::insertCol(std::vector<double> newCol, unsigned int col) {
	
	if(newCol.size() != this->getRows())
		return false;
		
	if(col < 0 or col > this->getCols())
		return false;
		
	for(unsigned int i = 0; i < this->getRows(); i++)
		matrix[i].insert(matrix[i].begin() + col, newCol[i]);
	cols++;
	return true;
}

bool Matrix::replaceCol(std::vector<double> newCol, unsigned int col) {
	
	if(col < 0 or col >= this->getCols())
		return false;
	
	if(newCol.size() != this->getRows())
		return false;
	
	for(unsigned int i = 0; i < this->getRows(); i++)
		matrix[i][col] = newCol[i];
		
	return true;
}

bool Matrix::replaceMatrixCol(Matrix matrixCol, unsigned int col) {
	
	if(col < 0 or col >= this->getCols())
		return false;
		
	if(matrixCol.getCols() != 1)
		return false;
		
	if(matrixCol.getRows() != this->getRows())
		return false;
		
	for(unsigned int i = 0; i < this->getRows(); i++)
		matrix[i][col] = matrixCol.at(i,0);
		
	return true;
}

bool Matrix::insertMatrixInRow(Matrix m, unsigned int row) {
	
	if(row < 0 or row > this->getRows())
		return false;
		
	if(this->getCols() != m.getCols())
		return false;
		
	for(unsigned int i = 0; i < m.getRows(); i++) {
		std::vector<double> row_vector = m.getRow(i);
		this->insertRow(row_vector, i + row);
	}
	
	return true;
}

bool Matrix::insertMatrixInCol(Matrix m, unsigned int col) {
	
	if(col < 0 or col > this->getCols())
		return false;
		
	if(this->getRows() != m.getRows())
		return false;
		
	for(unsigned int j = 0; j < m.getCols(); j++) {
		std::vector<double> col_vector = m.getCol(j);
		this->insertCol(col_vector, j + col);
	}
		
	return true;
}

void Matrix::roundMatrix(double difference) {
	
	for(unsigned int i = 0; i < this->getRows(); i++)
		for(unsigned int j = 0; j < this->getCols(); j++) {
			double round_value = std::round(matrix[i][j]);
			if(abs(round_value - matrix[i][j]) <= difference)
				matrix[i][j] = round_value;
		}
}

void Matrix::set(unsigned int i, unsigned int j, double x) {
	
	if(i < this->getRows() and j < this->getCols())
		matrix[i][j] = x;
	
}

bool Matrix::resize(unsigned int newRows, unsigned int newCols) {
	
	if(newRows != 0 and newCols != 0) {
		rows = newRows;
		cols = newCols;
		matrix.resize(newRows);
		for(unsigned int i = 0; i < newRows; i++)
			matrix[i].resize(newCols);
		return true;
	}
	return false;
	
}


//READ AND PRINT

void Matrix::read(void) {

	for(unsigned int i = 0; i < this->getRows(); i++)
		for(unsigned int j = 0; j < this->getCols(); j++)
			std::cin >> matrix[i][j];

}

void Matrix::print(void) const {
	
	for(unsigned int i = 0; i < rows; i++) {
		for(unsigned int j = 0; j < cols; j++)
			std::cout << matrix[i][j] << (j < cols - 1 ? "\t" : "");
		std::cout << std::endl;
	}
	
}

void Matrix::info(void) const {
	
	std::cout << "Dimension: " << this->getRows() << " x " << this->getCols() << std::endl;
	print();
}

//OPERATIONS

void Matrix::scalarProduct(double scalar) {
	
	for(unsigned int i = 0; i < this->getRows(); i++)
		for(unsigned int j = 0; j < this->getCols(); j++)
			matrix[i][j] *= scalar;
	
}

double Matrix::determinant(void) const {
	
	unsigned int rows = this->getRows();
	unsigned int cols = this->getCols();
	
	if(rows != cols)
		return 0;
		
	if(rows == 1)
		return matrix[0][0];
		
	if(rows == 2)
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		
	Matrix copy(this->getMatrix());
	copy.eraseCol(0);

	double det = 0;

	for(unsigned int i = 0; i < rows; i++) {
		if(matrix[i][0] != 0) {
			Matrix aux = copy;
			aux.eraseRow(i);
			det += (i % 2 == 0 ? 1 : -1) * matrix[i][0] * aux.determinant();
		}
	}
	
	return det;
}

Matrix Matrix::identity(unsigned int size) {
	
	Matrix id;
	
	if(size > 0) {
		std::vector<std::vector<double>> id_matrix(size, std::vector<double>(size));
		for(unsigned int i = 0; i < size; i++)
				id_matrix[i][i] = 1;
		id = Matrix(id_matrix); 
	}
	
	return id;
	
}

Matrix Matrix::random(unsigned int rows, unsigned int cols, int lower, int upper) {
	
	Matrix m;
	
	if(lower >= upper)
		return m;
	
	if(rows == 0 or cols == 0)
		return m;
	
	srand(time(NULL));
	std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
	for(unsigned int i = 0; i < rows; i++)
		for(unsigned int j = 0; j < cols; j++)
			matrix[i][j] = rand() % (upper - lower + 1) + lower;

	return m = Matrix(matrix);
}

Matrix Matrix::subMatrix(unsigned int i, unsigned int j, unsigned int rows, unsigned int cols) const {
	
	Matrix sub;
	
	if((i + rows - 1) < this->getRows() and (j + cols - 1) < this->getCols()) {
		sub = Matrix(rows, cols);
		for(unsigned int a = 0; a < rows; a++)
			for(unsigned int b = 0; b < cols; b++)
				sub.set(a, b, matrix[i+a][j+b]);		
	}
	
	return sub;	
}

Matrix Matrix::scalarProduct(Matrix m1, double scalar) {
	
	Matrix m;

	if(not m1.isBuilt())
		return m;
	
	std::vector<std::vector<double>> matrix = m1.getMatrix();
	
	for(unsigned int i = 0; i < matrix.size(); i++)
		for(unsigned int j = 0; j < matrix[i].size(); j++)
			matrix[i][j] *= scalar;
			
	return Matrix(matrix);
}

Matrix Matrix::matrixSum(Matrix m1, Matrix m2) {
	
	unsigned int rows = m1.getRows();
	unsigned int cols = m1.getCols();
	
	Matrix m;
	
	if(rows == 0 or cols == 0)
		return m;
		
	if(rows != m2.getRows() or cols != m2.getCols())
		return m;
		
	std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
	for(unsigned int i = 0; i < rows; i++)
		for(unsigned int j = 0; j < cols; j++)
			matrix[i][j] = m1.at(i,j) + m2.at(i,j);
	
	return m = Matrix(matrix);
}

Matrix Matrix::matrixProduct(Matrix m1, Matrix m2) {
	
	Matrix m;
	
	if(not m1.isBuilt() or not m2.isBuilt())
		return m;
	
	unsigned int m1_rows = m1.getRows();
	unsigned int m1_cols = m1.getCols();
	unsigned int m2_cols = m2.getCols();
	
	if(m1_cols != m2.getRows())
		return m;

	std::vector<std::vector<double>> matrix(m1_rows, std::vector<double>(m2_cols));
	
	for(unsigned int i = 0; i < m1_rows; i++)
		for(unsigned int j = 0; j < m2_cols; j++)
			for(unsigned int k = 0; k < m1_cols; k++)
				matrix[i][j] += m1.at(i,k) * m2.at(k,j);
				
	m = Matrix(matrix);
	
	return m;
}

Matrix Matrix::transpose(Matrix m) {
	
	Matrix tran;
	
	unsigned int rows = m.getCols();
	unsigned int cols = m.getRows();
	
	if(rows == 0 or cols == 0)
		return tran;
		
	tran = Matrix(rows, cols);
	for(unsigned int i = 0; i < rows; i++)
		for(unsigned int j = 0; j < cols; j++)
			tran.set(i,j,m.at(j,i));
		
	return tran;
}

Matrix Matrix::adjoint(Matrix m) {
	
	unsigned int rows = m.getRows();
	unsigned int cols = m.getCols();
	
	if(rows != cols)
		return Matrix(rows, rows);
		
	std::vector<vector<double>> matrix(rows, std::vector<double>(cols));
	
	for(unsigned int i = 0; i < rows; i++) {
		Matrix sub_matrix = m;
		sub_matrix.eraseRow(i);
		
		for(unsigned int j = 0; j < rows; j++) {
			Matrix aux = sub_matrix;
			aux.eraseCol(j);
			double sign = ((i % 2 == 0) ? 1 : -1) * ((j % 2 == 0) ? 1 : -1);
			matrix[i][j] = sign * aux.determinant();
		}
	}
	
	return Matrix(matrix);
}

Matrix Matrix::inverse(Matrix m) {
	
	unsigned int rows = m.getRows();
	unsigned int cols = m.getCols();
	
	if(rows != cols)
		return Matrix(rows, rows);

	double det = m.determinant();
	
	if(det == 0)
		return Matrix(rows, rows);
	
	return Matrix::scalarProduct(Matrix::adjoint(Matrix::transpose(m)), 1/det);
}

