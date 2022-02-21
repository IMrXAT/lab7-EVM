#include <iostream>
#include <utility>
#include <vector>
#include <ctime>
#define SIZE 2048
#define M 10
using namespace std;




vector<float> MatrixB(vector<float> matrixT);
vector<float>MatrixR(vector<float> matrixI, vector<float> matrixB, vector<float> matrix);
vector<float> NeededMatrix(vector<float> matrix, const vector<float>& matrixR, vector<float> matrixB);


vector<float> Mul(vector<float> matrix1, vector<float>matrix2);
vector<float>Sub(vector<float>matrix1, vector<float> matrix2);
vector<float>Add(vector<float> matrix1, vector<float> matrix2);
vector<float> Transposition(vector<float> matrix);


void PrintMatrix(vector<float> matrix);

// 1)searching A^t
// 2)searching B = A^t/(A_1*A_2)
// 2.1)searching A_1 and A_2
// 3)searching R = I-BA
// 4)searching A^-1 = (I+R+R^2+R^3+...)B
vector<float> MatrixReversal(const vector<float>& matrix){ // N - matrix SIZE, M - num iterations
    vector<float> matrixB(SIZE*SIZE); // its B
    vector<float>matrixI(SIZE*SIZE);
    vector<float> matrixR(SIZE*SIZE); // its R
    vector<float> matrixT(SIZE*SIZE);

    for (int i = 0; i < SIZE; i++){
        matrixI[i * SIZE + i] = 1;
    }
    //1
    matrixT = Transposition(matrix);
    //2
    matrixB = MatrixB(matrixT);
    //3
    matrixR = MatrixR(matrixI, matrixB, matrix);
    //4
    return NeededMatrix(matrixI, matrixR, matrixB);
}

void PrintMatrix(vector<float> matrix) {
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            cout << matrix[i*SIZE+j] << " ";
        }
        cout <<endl;
    }
    cout << "\n\n";
}

//////////////strange matrix and steps to win
vector<float> NeededMatrix(vector<float> matrix, const vector<float>& matrixR, vector<float>matrixB) {
    vector<float> matrixR2 = matrixR;
    for (int i = 1; i <= M; i++){
        cout << i;
        matrix = Add(matrix, matrixR2);
        matrixR2 = Mul(matrixR2, matrixR);
    }
    matrix = Mul(matrix, std::move(matrixB));
    return matrix;
}

vector<float>MatrixR(vector<float> matrixI, vector<float> matrixB, vector<float> matrix) {
    vector<float> matrix2;
    matrix2 = std::move(matrixI);
    return Sub(matrix2, Mul(std::move(matrix), std::move(matrixB)));
}



vector<float> MatrixB(vector<float> matrixT){
    float A_1 = -1, A_2 = -1, tmp1, tmp2;
    //2.1
    for (int i = 0; i < SIZE; i++){
        tmp1 = 0;
        tmp2 = 0;
        for (int j = 0; j < SIZE; j++){
            tmp1 += matrixT[i * SIZE + j];
            tmp2 += matrixT[j * SIZE + i];
        }
        if (tmp1 > A_1){
            A_1 = tmp1;
        }
        if (tmp2 > A_2){
            A_2 = tmp2;
        }
    }
    cout << A_1 << " "<< A_2 << endl;
    vector<float> matrix(SIZE*SIZE);
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++) {
            matrix[i*SIZE+j] = matrixT[i * SIZE + j] / (A_1 * A_2);
        }
    }
    return matrix;
}

////////////////////Matrix operations
vector<float> initMatrix(){
    vector<float> matrix(SIZE*SIZE);

    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            matrix[i*SIZE+j] = (float)(rand() % 10);
            //cin >> matrix[i*SIZE+j];
        }

    }
    return matrix;
}

vector<float> Transposition(vector<float> matrix) {
    for (int i = 0; i < SIZE; i++){
        for (int j = i; j < SIZE; j++){
            swap(matrix[i*SIZE+j], matrix[j*SIZE+i]);
        }
    }
    return matrix;
}

vector<float> Mul(vector<float> matrix1, vector<float> matrix2) {
    vector<float> matrix(SIZE*SIZE);
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            for (int k = 0; k < SIZE; k++){
                matrix[i*SIZE+j] += (matrix1[i*SIZE+k] * matrix2[k*SIZE+j]);
            }
        }
    }
    return matrix;
}

vector<float> Sub(vector<float> matrix1, vector<float> matrix2) {
    vector<float> matrix(SIZE*SIZE);
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            matrix[i*SIZE+j] = matrix1[i*SIZE+j] - matrix2[i*SIZE+j];
        }
    }
    return matrix;
}

vector<float> Add(vector<float> matrix1, vector<float> matrix2) {
    vector<float> matrix(SIZE*SIZE);
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            matrix[i*SIZE+j] = matrix1[i*SIZE+j] + matrix2[i*SIZE+j];
        }
    }
    return matrix;
}


int main(){
    vector<float> matrix(SIZE*SIZE);
    matrix = initMatrix();
    PrintMatrix(matrix);
    clock_t begin = clock();
    matrix = MatrixReversal(matrix);
    clock_t end = clock();
    cout << (end-begin)/CLOCKS_PER_SEC;
    PrintMatrix(matrix);
    return 0;
}
