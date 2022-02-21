#include <iostream>
#include <xmmintrin.h>
#include <ctime>
#define SIZE 2048
#define M 10
using namespace std;


float* MatrixB(float*  matrixT, float* matrix);
float* MatrixR(const float* matrixI, float* matrixB, float* matrix);
float*  NeededMatrix(float*matrix, float* matrixR, float* matrixB);


float* Mul(float* matrix1, float*matrix2);
float*Sub( float*matrix1, const float* matrix2);
float*Add(float* matrix1, const float* matrix2);
float* Transposition(const float* matrix);


void PrintMatrix(float* matrix);
float A_count(const float* matrix);
float value_from_vector(__m128 summ);
//float A_2_count(const float *matrix);
float *MatrixI();

// 1)searching matrix^t
// 2)searching B = matrix^t/(A_1*A_2)
// 2.1)searching A_1 and A_2
// 3)searching R = I-BA
// 4)searching matrix^-1 = (I+R+R^2+R^3+...)B

float* MatrixReversal( float* matrix){ // N - matrix SIZE, M - num iterations
    float* matrixB;
    float* matrixR;// its R
    float* matrixT;
    float* matrixI;

    matrixI = MatrixI();
    //1
    //PrintMatrix(matrixI);
    matrixT = Transposition(matrix);
    //2
    matrixB = MatrixB(matrixT, matrix);
//    cout << "matrixB\n";
//    PrintMatrix(matrixB);
    //3
    matrixR = MatrixR(matrixI, matrixB, matrix);
    //4
//    cout << "matrixR\n";
//    PrintMatrix(matrixR);
    return NeededMatrix(matrixI, matrixR, matrixB);
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
float* NeededMatrix(float* matrix, float* matrixR, float*matrixB) {
    float* matrixR2;
    matrixR2 = (float *)malloc(sizeof(float)*SIZE*SIZE);
    for (int i = 0; i < SIZE*SIZE; i++){
        matrixR2[i] = matrixR[i];
    }

    for (int i = 1; i <= M; i++){
        //cout << "matrixR^"<< i <<endl;
        //PrintMatrix(matrixR2);
        matrix = Add(matrix, matrixR2);
        matrixR2 = Mul(matrixR2, matrixR);
        cout << i<< endl;

        //cout << "not A^-1\n";
        //PrintMatrix(matrix);
    }
    matrix = Mul(matrix, matrixB);
    return matrix;
}

float* MatrixR(const float* matrixI, float* matrixB, float* matrix) {
    float* matrix2;
    matrix2 = (float*)malloc(sizeof(float)*SIZE*SIZE);
    for (int i = 0; i < SIZE* SIZE; i++){
        matrix2[i] = matrixI[i];
    }
    return Sub(matrix2, Mul(matrix, matrixB)); //осуждаю себя, но было лень писать нормально
}


float* MatrixB(float* matrixT, float* matrix){
    float A_1, A_2;
    //2.1

    A_1 = A_count(matrixT);
    A_2 = A_count(matrix); // передавать начальную матрицу
    float tmp[4];
    for (int i = 0; i < 4; i++){
        tmp[i] = A_1*A_2;
    }
    __m128 AA =((__m128*)tmp)[0];
//    float* matrixB;
//    matrixB = (float *) malloc(sizeof(float) * SIZE * SIZE);
    auto* matrixB = (__m128*)matrixT;

    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE/4; j++) {
            matrixB[i*SIZE/4+j] = _mm_div_ps( ((__m128*)matrixT)[i*SIZE/4+j], AA );
            //matrix[i*SIZE+j] = matrixT[i * SIZE + j] / (A_1 * A_2);
        }
    }

    return matrixT;
}

float A_count(const float *matrix) {
    float max = -1, new_max;
    __m128 summ;
    auto *AA = (__m128 *) matrix;
    for (int i = 0; i < SIZE; i++) {
        summ = _mm_setzero_ps(); //обнуляет все 4 значения вектора
        for (int j = 0; j < SIZE / 4; j++) {
            summ = _mm_add_ps(summ, AA[i * SIZE/4 + j]); //сложение
        }
        new_max = value_from_vector(summ);
        if (new_max >= max) {
            max = new_max;
        }
    }
    return max;
}

float value_from_vector(__m128 summ) {
    __m128 p;
    p = _mm_movehl_ps(p, summ); //копирование старшей половины в младшую
    summ = _mm_add_ps(summ, p); //сложение всех "координат" вектора
    p = _mm_shuffle_ps(summ, summ, 1); // перестановка упакованных значений
    summ = _mm_add_ss(summ, p); //сложение
    float result;
    _mm_store_ss(&result, summ); //записать младшее значение
    return result;
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

float* Transposition(const float* matrix) {
    float* transposedMatrix;
    transposedMatrix = (float *) malloc(sizeof(float) * SIZE * SIZE);

    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            transposedMatrix[i+j*SIZE] = matrix[j+i*SIZE];
        }
    }
    return transposedMatrix;
}

float* Mul(float* matrix1, float* matrix2) {
    float* matrix;
    matrix = (float *) malloc(sizeof(float) * SIZE * SIZE);
    auto *BB = (__m128 *) matrix2;
    __m128 mul;
    __m128 sum;
    for (int i = 0; i < SIZE; ++i) {
        float *c = matrix + i * SIZE;
        for (int k = 0; k < SIZE; ++k) {
            int b = k * SIZE / 4;
            __m128 a = _mm_set1_ps(matrix1[i * SIZE + k]); // устанавливает четыре позиции в одно значение
            for (int j = 0; j < SIZE / 4; j++) {
                mul = _mm_mul_ps(a, BB[b + j]); // умножение всех компонент вектора
                sum = _mm_loadu_ps(c + j * 4); //Загрузить четыре значения по невыровненному адресу
                sum = _mm_add_ps(sum, mul);// сложение всех компонент вектора
                _mm_storeu_ps(c + j * 4, sum); //записать четыре значения по невыровненному адресу
            }
        }
    }
    return matrix;
}

float* Sub( float* matrix1, const float* matrix2) {

    __m128 *intrMatrix1, *intrMatrix2;
    intrMatrix1 = (__m128*)matrix1;
    intrMatrix2 = (__m128*)matrix2;
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j <SIZE/4; j++){
            intrMatrix1[i*SIZE/4+j] = _mm_sub_ps( intrMatrix1[i*SIZE/4+j], intrMatrix2[i*SIZE/4+j]);
        }
    }
    return matrix1;
}

float*Add(float*matrix1, const float* matrix2) {
    __m128 *intrMatrix1, *intrMatrix2;
    intrMatrix1 = (__m128*)matrix1;
    intrMatrix2 = (__m128*)matrix2;
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j <SIZE/4; j++){
            intrMatrix1[i*SIZE/4+j] = _mm_add_ps( intrMatrix1[i*SIZE/4+j], intrMatrix2[i*SIZE/4+j]);
        }
    }
    return matrix1;
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
