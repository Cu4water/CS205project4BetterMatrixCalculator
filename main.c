#ifndef STDIO_H
#define STDIO_H
#include<stdio.h>
#endif
#ifndef STDLIB_H
#define STDLIB_H
#include<stdlib.h>
#endif
#ifndef MATRIX
#define MATRIX
#include"matrix.h"
#endif

int main() {
    struct Matrix *a = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    createMatrix(4, 4, a);
    int i=0;
    for(int r=1; r<=4; r++) {
        for(int c=1; c<=4; c++) {
            setNum(a, i, r, c);
            i++;
        }
    }
    struct Matrix *b = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    createMatrix(4, 4, b);
    i=0;
    for(int r=1; r<=4; r++) {
        for(int c=1; c<=4; c++) {
            setNum(b, i, r, c);
            i++;
        }
    }

    struct Matrix *ans = NULL;

    print(a);
    putchar('\n');

    ans=addMatrix(a, b, ans);
    print(ans);
    putchar('\n');

    ans=subtractMatrix(a, b, ans);
    print(ans);//ans是NULL也没关系哦，放心用就好啦
    putchar('\n');

    ans=multiplyMatrix(a, b, ans);
    print(ans);
    putchar('\n');

    ans=multiplyScaler(2, a, ans);
    print(ans);
    putchar('\n');

    printf("%f\n",max( multiplyMatrix(a, b, ans ) ) );
}