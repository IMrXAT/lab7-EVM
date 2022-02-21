#include <stdio.h>
#include <xmmintrin.h>
#define SIZE 4

float value_from_vector(__m128 s) {
    __m128 p;
    p = _mm_movehl_ps(p, s); //копирование старшей половины в младшую
    s = _mm_add_ps(s, p); //сложение всех "координат" вектора
    p = _mm_shuffle_ps(s, s, 1); // перестановка упакованных значений
    s = _mm_add_ss(s, p); //сложение
    float result;
    _mm_store_ss(&result, s); //записать младшее значение
    return result;
}

float *mul_matrix(const float *matrix1, const float *matrix2, int size) {
    float *matrix;
    matrix = (float *) calloc(size * size, sizeof(float));
    auto *BB = (__m128 *) matrix2;
    __m128 mul;
    __m128 sum;
    for (int i = 0; i < size; ++i) {
        float *c = matrix + i * size;
        for (int k = 0; k < size; ++k) {
            int b = k * size / 4;
            __m128 a = _mm_set1_ps(matrix1[i * size + k]); // устанавливает четыре позиции в одно значение
            for (int j = 0; j < size / 4; j++) {
                mul = _mm_mul_ps(a, BB[b + j]); // умножение всех компонент вектора
                sum = _mm_loadu_ps(c + j * 4); //Загрузить четыре значения по невыровненному адресу
                sum = _mm_add_ps(sum, mul);// сложение всех компонент вектора
                _mm_storeu_ps(c + j * 4, sum); //записать четыре значения по невыровненному адресу
            }
        }
    }
    return matrix;
}


float *sub_Matrix(const float *A, const float *B, int N) {
    auto *sub = (float *) calloc(N * N, sizeof(float));
    for (int i = 0; i < N * N; i++) {
        sub[i] = A[i] - B[i];
    }
    return sub;
}
float *add_Matrix(const float *A, const float *B, int N) {
    auto *sub = (float *) calloc(N * N, sizeof(float));
    for (int i = 0; i < N * N; i++) {
        sub[i] = A[i] + B[i];
    }
    return sub;
}

float *trans_Matrix(const float *A, int N) {
    auto *A_T = (float *) malloc(sizeof(float) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A_T[i + j * N] = A[j + i * N];
        }
    }
    return A_T;
}
float A_1_counting(const float *A, int N) {
    float max = -1, new_max;
    __m128 summ;
    auto *AA = (__m128 *) A;
    for (int i = 0; i < N; i++) {
        summ = _mm_setzero_ps(); //обнуляет все 4 значения вектора
        for (int j = 0; j < N / 4; j++) {
            summ = _mm_add_ps(summ, AA[i * N / 4 + j]); //сложение
        }
        new_max = value_from_vector(summ);
        if (new_max >= max) {
            max = new_max;
        }
    }
    return max;
}

float A_2_counting(const float *matrix, int size) {
    float max = -1, new_max = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            new_max += matrix[i + j * size];
        }
        if (new_max >= max) {
            max = new_max;
        }
        new_max = 0;
    }
    return max;

}
float *B_creating(const float *A_T, float A_1, float A_inf, int N) {
    auto *B = (float *) malloc(sizeof(float) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[j + i * N] = A_T[j + i * N] / (A_1 * A_inf);
        }
    }
    return B;
}
float *R_creating(float *A, float *B, int N) {
    float *I = (float *) calloc(N * N, sizeof(float));
    for (int i = 0; i < N; i++) {
        I[i + i * N] = 1;
    }
    float *BA = mul_matrix(B, A, N);
    float *R = sub_Matrix(I, BA, N);
    return R;
}


float *series_creating(float *R, int N, int M) {
    auto *matrixI = (float *) calloc(N * N, sizeof(float));
    for (int i = 0; i < N; i++) {
        matrixI[i + i * N] = 1;
    }
    auto *series = (float *) calloc(N * N, sizeof(float));
    series = add_Matrix(series, matrixI, N);
    auto *const_R = (float *) calloc(N * N, sizeof(float));
    const_R = add_Matrix(const_R, R, N);
    for (int i = 0; i < M - 1; i++) {
        R = mul_matrix(R, const_R, N);
        series = add_Matrix(series, R, N);
    }
    return series;
}

int main() {
//    FILE *fin;
//    fin = fopen("in.txt", "r");
    int M;
    //fscanf(fin, "%d %d", &N, &M);
    float *A;
    A = (float *) malloc(sizeof(float) * SIZE * SIZE);
    //get_Matrix(A, fin, SIZE);
    float *A_T;
    A_T = trans_Matrix(A, SIZE);
    float A_1 = A_1_counting(A, SIZE);
    float A_inf = A_2_counting(A, SIZE);
    float *B;
    B = B_creating(A_T, A_1, A_inf, SIZE);
    float *R;
    R = R_creating(A, B, SIZE);
    float *series = series_creating(R, SIZE, M);
    series = mul_matrix(series, B, SIZE);
    float *ans = mul_matrix(A, series, SIZE);
    free(series);
    free(R);
    free(B);
    free(A_T);
    free(A);
    return 0;

}