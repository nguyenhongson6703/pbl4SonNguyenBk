#include<iostream>
#include<math.h>
#include<cmath>
#include<Eigen/Dense>
using namespace std;
#define max 5
using namespace Eigen;
// phan tich eigencomposition la phai ma tran vuong moi thuc hien duoc 
// nhap ma tran vuong cap n
void nhapMaTran(float x[max][max], int& n) {
    do {
        cout << "\n Nhap vao cap ma tran vuong:";
        cin >> n;
    } while (n <= 0);
    cout << "\n Nhap vao ma tran: ";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "a[" << i  + 1 << "][" << j +1  << "] = ";
            cin >> x[i][j];
        }


    }

}
// xuat ma tran vuong cap n
void xuatMaTran(float x[max][max], int n) {
   
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "  a[" << i + 1 << "][" << j + 1 << "] = " << x[i][j] << " ";


        }
        cout << "\n";
    }
}
// nhan hai ma tran vuong
void  MultiMatrix(float A[max][max], float B[max][max], float C[max][max], int n){

		for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                C[i][j] = 0;
                for (int k = 0; k < n; k++)
                {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
		
	}

	

// ham swap cac cot
void swapColum(float matrix[max][max], int size, int m, int n) {
    for (int i = 0; i < size; i++) {
        float temp = matrix[i][m];
        matrix[i][m] = matrix[i][n];
        matrix[i][n] = temp;

    }
}
// ham swap cac dong
void swapRow(float matrix[max][max], int size, int m, int n) {
    for (int i = 0; i < size; i++) {
        float temp = matrix[m][i];
        matrix[m][i] = matrix[n][i];
        matrix[n][i] = temp;
    }
}
// ham tim vi tri 0 tren hang
int FindNonZeroCell(float matrix[max][max], int size, int i) {
    for (int j = 0; j < size; j++) {
        if (matrix[i][j] != 0) {
            return j;
        }
    }
    return -1;
}
// ham chuyen thanh ma tran tam giac tren
int toUpperTriangle(float matrix[max][max], int size) {
    int j = 0, rank = size;
    for (int i = 0; i < rank; i++) {
        while (rank > 0) {
            j = FindNonZeroCell(matrix, size,i);
            if (j == -1) {
                swapRow(matrix, size, i, --rank);
            }
            else {
                break;
            }
        }
        if (j != i) {
            swapColum(matrix, size, i, j);
        }
        float b = matrix[i][i];
        for (j = i + 1; j < rank; j++) {
            float a = matrix[j][i];
            matrix[j][i] = 0;
            int k;
            for (k = i + 1; k < rank; k++) {
                matrix[j][k] -= (a * matrix[i][k]) / b;
            }

        }
        while (rank > 0 && FindNonZeroCell(matrix, size, rank - 1) == -1) {
            rank--;
        }

    }
    return rank;

}
// ham tinh dinh thuc
float matrixDet(float matrix[max][max], int size) {
    //int rank = toUpperTriangle(matrix, size);
    float temp[max][max];
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            temp[i][j] = matrix[i][j];
        }
    }
    int rank = toUpperTriangle(temp, size);
    if (rank < size) {
        return 0;
    }
    float det = 1;
    for (int i = 0; i < size; i++) {
        det *= temp[i][i];

    }
    return det;

}

// ham tim ma tran nghich dao

MatrixXd RoudingMatrixValues(MatrixXd& matrix){
  	double min = 0.000001;
    for (int i = 0; i < matrix.rows(); i++) {
    	for (int j = 0; j < matrix.cols(); j++) {
    	    if (abs(matrix(i, j)) < min)
    	        matrix(i, j) = 0;
    	}
	}
  	return matrix;
}
void LamTronMatran(float a[max][max], int n){
	float min = 0.000001;
	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			if(abs(a[i][j])< min){
				a[i][j] = 0;
			}
		}
	}
}
MatrixXd ToInverse(const MatrixXd &matrix) {
	MatrixXd m = matrix.inverse();
	return RoudingMatrixValues(m);
}

void XacDinhPVaP1(float a[max][max], float P[max][max],float P1[max][max], int n, float D[max][max]){
	MatrixXd A;
	A.resize(n,n);
	for(int i = 0;i<n;i++){
		for(int j = 0;j<n;j++){
			A(i,j) = a[i][j];
		}
	}
	EigenSolver<MatrixXd> solver(A);
    
  
    
    MatrixXd eigenvectors = solver.eigenvectors().real();
    MatrixXd PInvert(n,n);
    PInvert = ToInverse(eigenvectors);
    

    int row = eigenvectors.rows();
    int col = eigenvectors.cols();
   
    for(int i = 0;i<row;i++){
    	for(int j = 0;j<col ;j++){
    		P[i][j] = eigenvectors(i,j);
    		P1[i][j] = PInvert(i,j);
    	
		}
	}
	
	
	float temp1[max][max];
	
	MultiMatrix(P1,a,temp1,n);
	MultiMatrix(temp1,P,D,n);
	
	

	
}

int main(){
	int n;
	float A[max][max];
	nhapMaTran(A,n);
	cout << "\n Ma tran A \n";
	xuatMaTran(A,n);
	float D[max][max];
	float P[max][max];
	float P1[max][max];
	XacDinhPVaP1(A,P,P1,n,D);
	cout << "\n Kiem tra gia tri";
	cout << "\n Ma tran P\n";
	LamTronMatran(P,n);
	xuatMaTran(P,n);
	cout << "\n Ma tran P^-1:\n";
	LamTronMatran(P1,n);
	xuatMaTran(P1,n);
	cout << "\n Ma tran D:\n";
	LamTronMatran(D,n);
	xuatMaTran(D,n);
	cout << "\n Kiem tra lai A*P\n";
	float temp1[max][max];
	MultiMatrix(A, P, temp1, n);
	LamTronMatran(temp1,n);
	xuatMaTran(temp1,n);
	float temp2[max][max];
	MultiMatrix(P,D,temp2,n);
	cout << "\n Kiem tra P*D\n";
	LamTronMatran(temp2,n);
	xuatMaTran(temp2,n);
	
	
}










