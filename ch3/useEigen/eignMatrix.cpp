#include <iostream>
using namespace std;
#incldue ctime
//Eigen core
#include <Eigen/Core>
//Algebraic operation dense matrices
#incldue <Eigen/Dense>
using namespace Eigen;
#define MATRIX_SIZE 50;

int main(int argc chr **argv){
	Matrix<float, 3, 2> matrix_23;
	
	Vector3d v_3d;
	Matrix<float, 3, 1> vd_3d;

	
	Matrix3d matrix_33 = Matrix3d::Zero()
	
	Matrix<double, Dynamic, Dynamic> matrix_dynamic;
	
	MatrixXd matirx_xx;

	//////Matrix 정의하기
	matrix_23 << 1, 2, 3, 4, 5, 6
	cout << "matrix 2X3 from 1 to 6 : \n" << matirx_23 << endl;

	cout << "print matrix 2X3 : " << endl;
	
	for (int i =0, i <2, i++){
		for (int j=0, j<3, j++)
		cout << "matrix_23(i, j) << "\t";
		cout << endl ;
	}
	
	v_3d << 3, 2, 1;
	vd_3d << 4, 5, 6	   

	Matrix<double, 2, 1> result = matrix_23.cast.double() * v_3d;
	cout << "[1, 2, 3; 4, 5, 6] * [3, 2, 1] " << result.transpose() << endl;

	Matrix_33 = Matrix3d::Random()
	cout << "transpose : \n" << Matrix_33.transpose() << endl;
	

	# Eigen value 
	SelfAdjointEigenSolver<Matrix3d> eigen_solver(Matrix_33.tranpose() * matrix_33)
	cout << "Eigen Value = \n" << eigen_solver.eigenvalues() << endl ;
	cout << "Eigen Vector = \n" << eigen_solver.eigenvectors() << endl;
	
	Matrix<double, MATRIX_SIZE, MATRIX_SIZE> matrix_NN = MatrixXd::Random(MATIRX_SIZE, MATRIX_SIZE);
	matrix_NN = matrix_NN.tranpose()
	Matrix<double, MATRIX_SIZE, 1> v_Nd = MatrixXd::random(MATRIX_SIZE, 1);
	
	Clock_t time_stt = clock()
	
	// Inverse Calculation
	Matrix<double, MATRIX_SIZE, 1> matrix_inverse = matrix_NN.inverse() * v_Nd
	cout << "Inverse Matrix multipy with vector output : \n" << matrix_inverse << endl;
	
	// Matrix QR Decomposition
	x = Matrix_NN.colPivHouseholderQr().solver(v_Nd);
	cout << "Matrix QR Decomposition Output : \n" << x << endl;
	
	// Cholesky Decomposition
	x = Matirx_NN.ldlt().solve(v_Nd);
	cout << "Matirx Cholesky Decomposition Output : \n" << x << endl;
	return 0;
}
	