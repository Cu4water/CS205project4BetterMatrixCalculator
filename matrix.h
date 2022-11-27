#ifndef STDIO_H
#define STDIO_H
#include<stdio.h>
#endif
#ifndef STDDEF_H
#define STDDEF_H
#include<stddef.h>
#endif

struct Matrix{
    size_t row,col;// 我们使用了[0,0]~[row-1,col-1]的内存,请使用get和set
    float *data;
};

//get, set
inline float get_num(const struct Matrix *now_used, const size_t row, const size_t col);
inline void set_num(struct Matrix *now_used, const float num, const size_t row, const size_t col);

//构造函数&析构函数
struct Matrix *create_file(struct Matrix *dest, FILE *source);
void create_mat(const size_t row, const size_t col, struct Matrix *used);
void delete_mat(struct Matrix *used);
struct Matrix *copy_mat(const struct Matrix *source, struct Matrix *dest);

//计算,请确保ans不等于a且ans不等于b
void inverse(struct Matrix *a);
struct Matrix *add_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);         //ans=a+b
struct Matrix *sub_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);         //ans=a-b
struct Matrix *mul_mat(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);         //ans=a*b,plain_version
struct Matrix *mul_ikj(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);         //寻址优化
struct Matrix *mul_avx(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);         //avx指令集
struct Matrix *mul_unroll(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);      //unroll
struct Matrix *mul_block(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);       //分治
struct Matrix *mul_omp(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);         //omp+block
struct Matrix *mul_strassen(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);    //strassen,会丢精度，不是omp的问题

struct Matrix *mul_improved(const struct Matrix *a, struct Matrix *b,struct Matrix *ans);           //转置优化无损乘

//project3的任务，但是在这里是石山
struct Matrix *addScaler(const float num, const struct Matrix *used, struct Matrix *ans);
struct Matrix *subtractScaler(const float num, const struct Matrix *used, struct Matrix *ans);
struct Matrix *multiplyScaler(const float num, const struct Matrix *used, struct Matrix *ans);
//最大最小
float max(const struct Matrix *used);
float min(const struct Matrix *used);

//输出
void print(const struct Matrix *used);
