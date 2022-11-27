#ifndef STDLIB_H
#define STDLIB_H
#include<stdlib.h>
#endif
#ifndef STRING_H
#define STRING_H
#include<string.h>
#endif
#ifndef X86INTRIN_H
#define X86INTRIN_H
#include<x86intrin.h>
#endif
#ifndef OMP_H
#define OMP_H
#include<omp.h>
#endif
#ifndef MATRIX
#define MATRIX
#include"matrix.h"
#endif

inline size_t get_position(const struct Matrix *now_used, const size_t row, const size_t col) {
    return (now_used -> col)*row+col;
}

inline float get_num(const struct Matrix *now_used, const size_t row, const size_t col) {
    if( row >= ( now_used -> row ) || col >= ( now_used -> col ) || row < 0 || col < 0 ) return 0.0;
    return now_used -> data[ get_position (now_used, row, col) ];
}
inline void set_num(struct Matrix *now_used, const float num, const size_t row, const size_t col) {
    if( row >= (now_used -> row) || col >= (now_used -> col) || row<0 || col<0 ) return ;
    now_used -> data[ get_position (now_used, row, col) ]=num;
}

//构造函数
struct Matrix *create_file(struct Matrix *dest, FILE *source){
    long long int r,c;
    fscanf(source,"%lld%lld",&r,&c);
    if(r <= 0 || c <= 0) {
        printf("illegal, r<=0 or c<=0\n");
        return NULL;
    }
    dest -> col = (size_t)c;
    dest -> row = (size_t)r;
    dest -> data = (float *)malloc(r*c*sizeof(float));
    for(size_t i=0; i<r; ++i) {
        for(size_t j=0; j<c; ++j) {
            fscanf(source, "%f", &(dest -> data[get_position(dest,i,j)]));
        }
    }
    fclose(source);//关闭文件
    return dest;

}
void create_mat(const size_t row, const size_t col, struct Matrix *used) {
    if( row == 0 || col == 0) return;
    if( used == NULL ) return;
    used -> col = col;
    used -> row = row;
    used -> data = (float *)malloc(row*col*sizeof(float));
}
void delete_mat(struct Matrix *used) {
    if(used == NULL) return;
    free(used -> data);
    free(used);
}
struct Matrix *copy_mat(const struct Matrix *source, struct Matrix *dest) {
    if(source == NULL) return NULL;
    delete_mat(dest);
    size_t row=source->row;
    size_t col=source->col;
    create_mat(row, col, dest);
    memcpy(dest->data, source->data, row*col*sizeof(float));
    return dest;
}

//计算
inline void swap(float *a,float *b) {
    float tmp=*a;
    *a=*b;
    *b=tmp;
}

void inverse(struct Matrix *b) {
    struct Matrix *tmp = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    copy_mat(b,tmp);
    size_t row = b -> col;
    size_t col = b -> row;
    b -> row = row;
    b -> col = col;
    for(register size_t i=0; i<row; ++i) {
        for(register size_t j=0; j<col; ++j) {
            set_num(b, get_num(tmp, i, j), j, i);
        }
    }
    delete_mat(tmp);
}
struct Matrix *add_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if( a -> row != b -> row || a -> col != b -> col ) return NULL;
    delete_mat(ans);
    create_mat(a -> row, a -> col, ans);
    for(size_t r=0; r<ans -> row; r++) {
        for(size_t c=0; c<ans -> col; c++) {
            set_num(ans, get_num(a, r, c) + get_num(b, r, c), r, c);
        }
    }
    return ans;
}
struct Matrix *sub_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if( a -> row != b -> row || a -> col != b -> col ) return NULL;
    delete_mat(ans);
    create_mat(a -> row, a -> col, ans);
    for(size_t r=0; r<ans -> row; r++) {
        for(size_t c=0; c<ans -> col; c++) {
            set_num(ans, get_num(a, r, c) - get_num(b, r, c), r, c);
        }
    }
    return ans;
}

struct Matrix *mul_plain(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
//我知道这可能使您血压升高，但我们管这个叫‘baseline’
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if(a -> col != b -> row) return NULL;
    delete_mat(ans);
    create_mat(a -> row, b -> col, ans);
    for(size_t r = 0; r < ans -> row; r++) {
        for(size_t c = 0; c < ans -> col; c++) {
            float ans_tmp=0;
            for(size_t k = 0; k < a->col; k++) {
                ans_tmp += get_num(a, r, k) * get_num(b, k, c);
            }
            set_num(ans, ans_tmp, r, c);
        }
    }
    return ans;
}
struct Matrix *mul_ikj(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if(a -> col != b -> row) return NULL;
    delete_mat(ans);
    create_mat(a -> row, b -> col, ans);
    float ark = 0;
    register size_t r,c,k;
    for(r = 0; r < ans -> row; ++r) {
        for(k = 0; k < a -> col; ++k) {
            ark = get_num(a, r, k);
            for(c = 0; c < ans -> col; ++c) {
                set_num(ans, 
                ark*get_num(b, k, c)+ get_num(ans, r, c), 
                r, c);
            }
        }
    }
    return ans;
}
struct Matrix *mul_avx(const struct Matrix *a, struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if(a -> col != b -> row) return NULL;
    delete_mat(ans);
    create_mat(a -> row, b -> col, ans);
    inverse(b);
    register size_t i, j, k;
    register float tmp[8],ans_tmp;
    __m256 X,Y;
    __m256 acc;
    //经过转置，c[i][j]=a.row[i]*b.row[j];
    for(i = 0; i < a -> row; ++i) {
        for(j = 0; j < b -> row; ++j) {
            acc = _mm256_setzero_ps();
            for(k = 0; k + 8 < a -> col; k += 8) {
                X = _mm256_loadu_ps(a -> data + get_position(a,i,k));
                Y = _mm256_loadu_ps(b -> data + get_position(b,j,k));
                acc = _mm256_add_ps(acc, _mm256_mul_ps(X, Y));
            }
            _mm256_storeu_ps(tmp, acc);
            ans->data[i*ans->col+j]=tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7];
            for(; k < a -> col; ++k) {
                ans->data[i*ans->col+j] += get_num(a, i, k) * get_num(b, j, k);
            }
        }
    }
    inverse(b);
    return ans;
}
#define ROLL 4
struct Matrix *mul_unroll(const struct Matrix *a, struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if(a -> col != b -> row) return NULL;
    delete_mat(ans);
    create_mat(a -> row, b -> col, ans);
    inverse(b);
    register size_t i, j, k, l;
    register float tmp[ROLL][8],ans_tmp;
    __m256 X,Y;
    __m256 acc[ROLL];
    //经过转置，c[i][j]=a.row[i]*b.row[j];
    for(i = 0; i < a -> row; i+=ROLL) {
        for(j = 0; j < b -> row; ++j) {
            for(l = 0; l < ROLL; ++l) {
                acc[l] = _mm256_setzero_ps();
            }
            for(k = 0; k + 8 < a -> col; k += 8) {
                for(l = 0; l < ROLL; ++l) {
                    X = _mm256_loadu_ps(a -> data + get_position(a,i+l,k));
                    Y = _mm256_loadu_ps(b -> data + get_position(b,j,k));
                    acc[l] = _mm256_add_ps(acc[l], _mm256_mul_ps(X, Y));
                }
            }
            for(l = 0; l < ROLL; ++l) {
                _mm256_storeu_ps(tmp[l], acc[l]);
                ans->data[i*ans->col+j]=tmp[l][0]+tmp[l][1]+tmp[l][2]+tmp[l][3]+tmp[l][4]+tmp[l][5]+tmp[l][6]+tmp[l][7];
            }
            for(; k < a -> col; ++k) {
                ans->data[i*ans->col+j] += get_num(a, i, k) * get_num(b, j, k);
            }
        }
    }
    inverse(b);
    return ans;
}
struct Matrix *mul_omp(const struct Matrix *a, struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if(a -> col != b -> row) return NULL;
    delete_mat(ans);
    create_mat(a -> row, b -> col, ans);
    //inverse
    struct Matrix *ttmp = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    copy_mat(b,ttmp);
    size_t row = b -> col;
    size_t col = b -> row;
    b -> row = row;
    b -> col = col;
    #pragma omp parallel for
    for(register size_t i=0; i<row; ++i) {
        for(register size_t j=0; j<col; ++j) {
            set_num(b, get_num(ttmp, i, j), j, i);
        }
    }
    register size_t i, j, k;
    register float tmp[8],ans_tmp;
    __m256 X,Y;
    __m256 acc;
    //经过转置，c[i][j]=a.row[i]*b.row[j];
    for(i = 0; i < a -> row; ++i) {
        for(j = 0; j < b -> row; ++j) {
            acc = _mm256_setzero_ps();
            for(k = 0; k + 8 < a -> col; k += 8) {
                X = _mm256_loadu_ps(a -> data + get_position(a,i,k));
                Y = _mm256_loadu_ps(b -> data + get_position(b,j,k));
                acc = _mm256_add_ps(acc, _mm256_mul_ps(X, Y));
            }
            _mm256_storeu_ps(tmp, acc);
            ans->data[i*ans->col+j]=tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7];
            for(; k < a -> col; ++k) {
                ans->data[i*ans->col+j] += get_num(a, i, k) * get_num(b, j, k);
            }
        }
    }
    //inverse
    copy_mat(b,ttmp);
    size_t row = b -> col;
    size_t col = b -> row;
    b -> row = row;
    b -> col = col;
    #pragma omp parallel for
    for(register size_t i=0; i<row; ++i) {
        for(register size_t j=0; j<col; ++j) {
            set_num(b, get_num(ttmp, i, j), j, i);
        }
    }
    delete_mat(ttmp);
    return ans;
}
struct Matrix *mul_strassen(const struct Matrix *a, struct Matrix *b, struct Matrix *ans, size_t size) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if(a -> col != b -> row) return NULL;
    delete_mat(ans);
    create_mat(a -> row, b -> col, ans);
    size_t row=a->row, col=b->col;
    size_t nxtsize=size>>1;

    if(row<=128) {//分治结束，灵感来源：stl sort
        return mul_ikj(a,b,ans);
    }

    struct Matrix *A11 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *A12 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *A21 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *A22 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *B11 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *B12 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *B21 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *B22 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *C11 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *C12 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *C21 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *C22 = (struct Matrix*) malloc(sizeof(struct Matrix));

    create_mat(nxtsize, nxtsize, A11);
    create_mat(nxtsize, nxtsize, A12);
    create_mat(nxtsize, nxtsize, A21);
    create_mat(nxtsize, nxtsize, A22);
    create_mat(nxtsize, nxtsize, B11);
    create_mat(nxtsize, nxtsize, B12);
    create_mat(nxtsize, nxtsize, B21);
    create_mat(nxtsize, nxtsize, B22);
    create_mat(nxtsize, nxtsize, C11);
    create_mat(nxtsize, nxtsize, C12);
    create_mat(nxtsize, nxtsize, C21);
    create_mat(nxtsize, nxtsize, C22);

    struct Matrix *resA = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *resB = (struct Matrix*) malloc(sizeof(struct Matrix));

    struct Matrix *S1 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S2 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S3 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S4 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S5 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S6 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S7 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S8 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S9 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *S10 = (struct Matrix*) malloc(sizeof(struct Matrix));

    struct Matrix *P1 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *P2 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *P3 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *P4 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *P5 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *P6 = (struct Matrix*) malloc(sizeof(struct Matrix));
    struct Matrix *P7 = (struct Matrix*) malloc(sizeof(struct Matrix));

    create_mat(nxtsize, nxtsize, resA);
    create_mat(nxtsize, nxtsize, resB);
    create_mat(nxtsize, nxtsize, S1);
    create_mat(nxtsize, nxtsize, S2);
    create_mat(nxtsize, nxtsize, S3);
    create_mat(nxtsize, nxtsize, S4);
    create_mat(nxtsize, nxtsize, S5);
    create_mat(nxtsize, nxtsize, S6);
    create_mat(nxtsize, nxtsize, S7);
    create_mat(nxtsize, nxtsize, S8);
    create_mat(nxtsize, nxtsize, S9);
    create_mat(nxtsize, nxtsize, S10);
    create_mat(nxtsize, nxtsize, P1);
    create_mat(nxtsize, nxtsize, P2);
    create_mat(nxtsize, nxtsize, P3);
    create_mat(nxtsize, nxtsize, P4);
    create_mat(nxtsize, nxtsize, P5);
    create_mat(nxtsize, nxtsize, P6);
    create_mat(nxtsize, nxtsize, P7);

    #pragma omp parallel for
    for(size_t r = 0; r < nxtsize; ++r) {
        for(size_t c = 0; c < nxtsize; ++c) {
            A11 -> data[r * nxtsize + c] = a -> data[r * col + c];
            A12 -> data[r * nxtsize + c] = a -> data[r * col + c + ( row >> 1 )];
            A21 -> data[r * nxtsize + c] = a -> data[(r + nxtsize) * col + c];
            A22 -> data[r * nxtsize + c] = a -> data[(r + nxtsize) * col + c + nxtsize];
            B11 -> data[r * nxtsize + c] = b -> data[r * col + c];
            B12 -> data[r * nxtsize + c] = b -> data[r * col + c + (row >> 1)];
            B21 -> data[r * nxtsize + c] = b -> data[(r + nxtsize) * col + c];
            B22 -> data[r * nxtsize + c] = b -> data[(r + nxtsize) * col + c + (row >> 1)];
        }
    }

    S1 = mat_sub(B11, B22, S1);
    S2 = mat_sub(A11, A12, S2);
    S3 = mat_sub(A21, A22, S3);
    S4 = mat_sub(B21, B11, S4);
    S5 = mat_sub(A11, A22, S5);
    S6 = mat_sub(B11, B22, S6);
    S7 = mat_sub(A12, A22, S7);
    S8 = mat_sub(B21, B22, S8);
    S9 = mat_sub(A11, A21, S9);
    S10 = mat_sub(B11, B12, S10);

    P1 = mul_strassen(A11, S1, P1, nxtsize);
    P2 = mul_strassen(S2, B22, P2, nxtsize);
    P3 = mul_strassen(S3, B11, P3, nxtsize);
    P4 = mul_strassen(A22, S3, P4, nxtsize);
    P5 = mul_strassen(S5, S6, P5, nxtsize);
    P6 = mul_strassen(S7, S8, P6, nxtsize);
    P7 = mul_strassen(S9, S10, P7, nxtsize);

    resA = add_mat(P4, P5, resA);
    resB = sub_mat(resA, P2, resB);
    C11 = add_mat(resB, P6, C11);
    C12 = add_mat(P1, P2, C12);
    C21 = add_mat(P3, P4, C21);
    resA = add_mat(P1, P5, resA);
    resB = sub_mat(resA, P3, resB);
    C22 = sub_mat(resB, P7, C22);

    for(size_t r = 0; r < nxtsize; ++r) {
        for(size_t c = 0; c < nxtsize;  ++c) {
            ans -> data[r * col + c] = C11 -> data[r * nxtsize + c];
            ans -> data[r * col + c + ( row >> 1 )] = C12 -> data[r * nxtsize + c];
            ans -> data[(r + ( row >> 1 )) * col + c] = C21 -> data[r * nxtsize + c];
            ans -> data[(r + ( row >> 1 )) * col + c + ( row >> 1 )] = C22 -> data[r * nxtsize + c];
        }
    }

    delete_mat(A11);
    delete_mat(A12);
    delete_mat(A21);
    delete_mat(A22);
    delete_mat(B11);
    delete_mat(B12);
    delete_mat(B21);
    delete_mat(B22);
    delete_mat(C11);
    delete_mat(C12);
    delete_mat(C21);
    delete_mat(C22);
    delete_mat(resA);
    delete_mat(resB);
    delete_mat(S1);
    delete_mat(S2);
    delete_mat(S3);
    delete_mat(S4);
    delete_mat(S5);
    delete_mat(S6);
    delete_mat(S7);
    delete_mat(S8);
    delete_mat(S9);
    delete_mat(S10);
    delete_mat(P1);
    delete_mat(P2);
    delete_mat(P3);
    delete_mat(P4);
    delete_mat(P5);
    delete_mat(P6);
    delete_mat(P7);

    return ans;
}

struct Matrix *addScaler(const float num, const struct Matrix *used, struct Matrix *ans) {
    if(used == NULL) return NULL;
    ans = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    copyMatrix(used, ans);
    for(int r=1; r<=ans -> row; r++) {
        for(int c=1; c<=ans -> col; c++) {
            set_num(ans, get_num(ans, r, c) + num, r, c);
        }
    }
    return ans;
}
struct Matrix *subtractScaler(const float num, const struct Matrix *used, struct Matrix *ans) {
    if(used == NULL) return NULL;
    ans = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    copyMatrix(used, ans);
    for(int r=1; r<=ans -> row; r++) {
        for(int c=1; c<=ans -> col; c++) {
            set_num(ans, get_num(ans, r, c) - num, r, c);
        }
    }
    return ans; 
}
struct Matrix *multiplyScaler(const float num, const struct Matrix *used, struct Matrix *ans) {
    if(used == NULL) return NULL;
    ans = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    copyMatrix(used, ans);
    for(int r=1; r<=ans -> row; r++) {
        for(int c=1; c<=ans -> col; c++) {
            set_num(ans, get_num(ans, r, c) * num, r, c);
        }
    }
    return ans;
}

//最大最小
float max(const struct Matrix *used) {
    if(used == NULL) return 0.0;
    float ans = get_num(used, 1, 1);
    for(int r=1; r<=used -> row; r++) {
        for(int c=1; c<=used -> col; c++) {
            if(get_num(used, r, c) > ans) {
                ans =get_num(used, r, c);
            }
        }
    }
    return ans;
}
float min(const struct Matrix *used) {
    if(used == NULL) return 0.0;
    float ans = get_num(used, 1, 1);
    for(int r=1; r<=used -> row; r++) {
        for(int c=1; c<=used -> col; c++) {
            if(get_num(used, r, c) < ans) {
                ans =get_num(used, r, c);
            }
        }
    }
    return ans;
}

//输出
void print(const struct Matrix *used) {
    if(used == NULL) return;
    for(int r=1; r<=used -> row; r++) {
        for(int c=1; c<=used -> col; c++) {
            printf("%f ", get_num(used, r, c));
        }
        printf("\n");
    }
    return ;
}