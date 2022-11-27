#ifndef STDLIB_H
#define STDLIB_H
#include<stdlib.h>
#endif
#ifndef STRING_H
#define STRING_H
#include<string.h>
#endif
#ifndef OMP_H
#define OMP_H
#include<omp.h>
#endif
#ifndef MATRIX
#define MATRIX
#include"matrix.h"
#endif

inline size_t getPosition(const struct Matrix *now_used, const size_t row, const size_t col) {
    return (now_used -> col)*row+col;
}

inline float get_num(const struct Matrix *now_used, const size_t row, const size_t col) {
    if( row >= ( now_used -> row ) || col >= ( now_used -> col ) || row < 0 || col < 0 ) return 0.0;
    return now_used -> data[ getPosition (now_used, row, col) ];
}
inline void set_num(struct Matrix *now_used, const float num, const size_t row, const size_t col) {
    if( row >= (now_used -> row) || col >= (now_used -> col) || row<0 || col<0 ) return ;
    now_used -> data[ getPosition (now_used, row, col) ]=num;
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
            fscanf(source, "%f", &(dest -> data[getPosition(dest,i,j)]));
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

void inverse(struct Matrix *a) {
    struct Matrix tmp;
    copy_mat(a,&tmp);
    size_t row = a -> col;
    size_t col = a -> row;
    a -> row = col;
    a -> col = row;
    for(size_t i=0; i<row; ++i) {
        for(size_t j=0; j<col; ++j) {
            set_num(a,get_num(&tmp,j,i),i,j);
        }
    }
}
struct Matrix *add_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if( a -> row != b -> row || a -> col != b -> col ) return NULL;
    deleteMatrix(ans);
    createMatrix(a -> row, a -> col, ans);
    for(size_t r=0; r<ans -> row; r++) {
        for(size_t c=0; c<ans -> col; c++) {
            setNum(ans, getNum(a, r, c) + getNum(b, r, c), r, c);
        }
    }
    return ans;
}
struct Matrix *sub_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    if(a == NULL || b == NULL || a == ans || b == ans) return NULL;
    if( a -> row != b -> row || a -> col != b -> col ) return NULL;
    deleteMatrix(ans);
    createMatrix(a -> row, a -> col, ans);
    for(size_t r=0; r<ans -> row; r++) {
        for(size_t c=0; c<ans -> col; c++) {
            setNum(ans, getNum(a, r, c) - getNum(b, r, c), r, c);
        }
    }
    return ans;
}
struct Matrix *mul_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    deleteMatrix(ans);
    ans =  (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    if(a == NULL || b == NULL) return NULL;
    if(a -> col != b -> row) return NULL;
    createMatrix(a -> row, b -> col, ans);
    for(int r=1; r<=ans -> row; r++) {
        for(int c=1; c<=ans -> col; c++) {
            float ans_tmp=0;
            for(int k=1; k<=a->col;k++) {
                ans_tmp+=getNum(a, r, k) * getNum(b, k, c);
            }
            setNum(ans, ans_tmp, r, c);
        }
    }
    return ans;
}
struct Matrix *addScaler(const float num, const struct Matrix *used, struct Matrix *ans) {
    if(used == NULL) return NULL;
    ans = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    copyMatrix(used, ans);
    for(int r=1; r<=ans -> row; r++) {
        for(int c=1; c<=ans -> col; c++) {
            setNum(ans, getNum(ans, r, c) + num, r, c);
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
            setNum(ans, getNum(ans, r, c) - num, r, c);
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
            setNum(ans, getNum(ans, r, c) * num, r, c);
        }
    }
    return ans;
}

//最大最小
float max(const struct Matrix *used) {
    if(used == NULL) return 0.0;
    float ans = getNum(used, 1, 1);
    for(int r=1; r<=used -> row; r++) {
        for(int c=1; c<=used -> col; c++) {
            if(getNum(used, r, c) > ans) {
                ans =getNum(used, r, c);
            }
        }
    }
    return ans;
}
float min(const struct Matrix *used) {
    if(used == NULL) return 0.0;
    float ans = getNum(used, 1, 1);
    for(int r=1; r<=used -> row; r++) {
        for(int c=1; c<=used -> col; c++) {
            if(getNum(used, r, c) < ans) {
                ans =getNum(used, r, c);
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
            printf("%f ", getNum(used, r, c));
        }
        printf("\n");
    }
    return ;
}