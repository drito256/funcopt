#ifndef MATRIX_H
#define MATRIX_H

#include "helper.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string_view>
#include <iomanip>

class Matrix{
    public:
        Matrix(size_t rows, size_t columns, std::unique_ptr<double[]> data);
        Matrix(size_t rows, size_t columns, std::vector<double> data);
        Matrix(const std::string& filename);
        Matrix(const Matrix& other);
        Matrix(std::vector<double> v);
        Matrix() = delete;
       ~Matrix() = default;

       size_t getRows() const;
       size_t getColumns() const;
    
       void print() const;
       void save(const std::string& filename);
       void read_file(const std::string& filename);

       // ostale metode za supstituciju , dekompoziciju, determinantu
       Matrix transpose() const;
       Matrix inverse();
       Matrix subs_backward(const Matrix& b) const;
       Matrix subs_forward(const Matrix& y) const;
       Matrix LU_decomp() const;
       std::pair<Matrix, Matrix> LUP_decomp();
       Matrix extract_L() const;
       Matrix extract_U() const;
       Matrix solve_w_LU(const Matrix& vec) const;
       Matrix solve_w_LUP(const Matrix& vec);  // not const cause of swap_count
       double det();

       // assignment operator
       Matrix& operator=(const Matrix& matrix);

       // arithmetic operators
       Matrix& operator+=(const Matrix& m1);
       Matrix& operator-=(const Matrix& m1);
       Matrix& operator*=(const double scalar);
       Matrix& operator/=(const double scalar);
       
       double operator()(int row, int col) const;
       double& operator()(int row, int col);
       
       // comparison operator
       friend bool operator==(const Matrix& m1, const Matrix& m2);
       
       //arithemtic operators
       friend Matrix operator+(const Matrix& m1, const Matrix& m2);
       friend Matrix operator-(const Matrix& m1, const Matrix& m2);
       friend Matrix operator*(const Matrix& m1, const Matrix& m2);
       friend Matrix operator*(const Matrix& m, const double scalar);

    private:
        size_t m_rows{0}, m_columns{0};
        size_t m_swap_count{0};
        std::unique_ptr<double[]> m_data; // mozda shared

        void swap_rows(int row1, int row2);
};

#endif
