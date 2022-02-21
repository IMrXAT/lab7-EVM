#include <iostream>
#include <cblas.h>
#include <ctime>
#include <cstdlib>
#include <cstring>
#define SIZE 2048
#define M 10
using namespace std;


float* MatrixB(float* matrix);
float* MatrixR(const float* matrixI, float* matrixB, float* matrix);
float*  NeededMatrix(float*matrix, float* matrixR, float* matrixB);


float* Mul(float* matrix1, float*matrix2);
float*Sub( float*matrix1, const float* matrix2);
float*Add(float* matrix1, const float* matrix2);
float* Transposition(const float* matrix);


void PrintMatrix(float* matrix);
//float A_2_count(const float *matrix);
float *MatrixI();

// 1)searching B = matrix^t/(A_1*A_2)
// 1.1)searching A_1 and A_2
// 2)searching R = I-BA
// 3)searching matrix^-1 = (I+R+R^2+R^3+...)B

float A_1_count(float *matrix);

float A_2_count(float *matrix);

float* MatrixReversal(float* matrix){ // N - matrix SIZE, M - num iterations
    auto *mem = (float*) malloc(4*SIZE*SIZE*sizeof(float ));
    float* matrixA = mem;
    float* matrixB = mem+1*SIZE*SIZE;
    float* matrixR = mem+2*SIZE*SIZE;// its R
    float* matrixM = mem+3*SIZE*SIZE;
    float* matrixI;
    float* matrixT = mem+4*SIZE*SIZE;

    matrixI = MatrixI();
    //1
    //PrintMatrix(matrixI);
    float A_1, A_2;
    A_1 = A_1_count(matrix);
    A_2 = A_2_count(matrix);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, SIZE, SIZE, SIZE,
                1.0f/A_1/A_2, matrix, SIZE, matrix, SIZE, 0, matrixR, SIZE);

    cblas_saxpy(SIZE * SIZE, -1.0f, matrixI, 1, matrixR, 1);
    memcpy(matrixM, matrixR, SIZE * SIZE * sizeof(float));
    for (int i = 0; i < M; ++i) {
        cblas_saxpy(SIZE * SIZE, 1.0f, matrixM, 1, matrixI, 1);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, SIZE, SIZE, SIZE,
                    1.0f, matrixM, SIZE, matrixR, SIZE, 0, matrixM, SIZE);
        //matrixM = matrixT;
    }
    //4
//    cout << "matrixR\n";
//    PrintMatrix(matrixR);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, SIZE, SIZE, SIZE,
                1.0f, matrixA, SIZE, matrixI, SIZE, 0, matrixT, SIZE);
    return matrixT;
}

float *MatrixI() {
    float *matrix;
    matrix =(float*) calloc(SIZE*SIZE, sizeof(float));
    for (int i = 0; i < SIZE; i++){
        matrix[i*SIZE+i] = 1;
    }
    return matrix;
}

void PrintMatrix(float* matrix) {
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            cout << matrix[i*SIZE+j] << " ";
        }
        cout <<endl;
    }
    cout << "\n";
}

//////////////strange matrix and steps to win
float* MatrixB(float* matrix){

    return matrix;

}

float A_2_count(float *matrix) {
    float sums[SIZE];
    for (int i = 0; i < SIZE; ++i) {
        sums[i] = cblas_sasum(SIZE, matrix + i, SIZE);
    }
    return sums[cblas_isamax(SIZE, sums, 1)];
}

float A_1_count(float *matrix) {
    float sums[SIZE];
    for (int i = 0; i < SIZE; ++i) {
        sums[i] = cblas_sasum(SIZE, matrix + SIZE * i, 1);
    }
    return sums[cblas_isamax(SIZE, sums, 1)];
}

////////////////////Matrix operations
float *initMatrix(){
    float* matrix;
    matrix = (float*) malloc(sizeof(float )*SIZE*SIZE);
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            matrix[i*SIZE+j] = (float)(rand() % 10);
        }
    }
    return  matrix;
}


int main(){
    float* matrix;
    matrix = initMatrix();
    //PrintMatrix(matrix);
    clock_t begin = clock();
    matrix = MatrixReversal(matrix);
    clock_t end = clock();
    //PrintMatrix(matrix);
    cout << (float)(end - begin)/CLOCKS_PER_SEC;
    return 0;
}
