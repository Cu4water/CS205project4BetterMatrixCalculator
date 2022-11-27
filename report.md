# CS205 C/C++ Programming Report Of Project 3


**Name:** 杨博乔  
**SID:** 12112805  
## Part 1. analysis  
> 题目需要仅用C设计一个矩阵类，重点要易用而非效率，因此，需要强大的鲁棒性

显然，用户可能出错的大抵是不可计算，或指针错误这两种情况，指针错误是无法鉴别的（仅限于把int指针当一个matrix指针传到函数里），但是指针是NULL这种问题可以由开发者来帮助用户解决。  
此外，由于人对于指针的掌握程度是具有顺序的，因此对于参数都能传成NULL的用户，二级指针是可怕的，我使用特殊的回传参数避免了这个问题。  

## Part 2. code

```cpp
struct Matrix *addMatrix(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans);       //ans=a+b
```
这段代码会将ans作为返回值，对于错误的输入，则会返回NULL,其实现是  
```cpp
struct Matrix *addMatrix(const struct Matrix *a, const struct Matrix *b, struct Matrix *ans) {
    deleteMatrix(ans);
    ans = (struct Matrix *) malloc ( sizeof( struct Matrix ) );
    if(a == NULL || b == NULL) return NULL;
    if( a -> row != b -> row || a -> col != b -> col ) return NULL;
    createMatrix(a -> row, a -> col, ans);
    for(int r=1; r<=ans -> row; r++) {
        for(int c=1; c<=ans -> col; c++) {
            setNum(ans, getNum(a, r, c) + getNum(b, r, c), r, c);
        }
    }
    return ans;
}
```
同理，对于标量计算，我们有（这里仍以加法为例）  
```cpp
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
```  
  
其余计算部分代码类似，Matrix结构的实现为
```cpp
struct Matrix{
    int row,col;// 我们使用了[0,0]~[row-1,col-1]的内存，您可以调用[1,1]~[row,col]的数据,请使用get和set
    float *data;
};
```
矩阵的生成与释放
```cpp
void createMatrix(const int row, const int col, struct Matrix *used) {
    if( row <= 0 || col <= 0) return;
    if(used == NULL) return;
    used -> col = col;
    used -> row = row;
    used -> data = (float *) malloc( sizeof(float) * row * col );
}
void deleteMatrix(struct Matrix *used) {
    if(used == NULL) return;
    free(used -> data);//由于需要用到used域中的东西，所以需要预处理NULL
    free(used);
}
struct Matrix *copyMatrix(const struct Matrix *source, struct Matrix *dest) {
    if(source == NULL || dest == NULL) return NULL;
    deleteMatrix(dest);
    createMatrix(source -> row, source -> col, dest);
    for(int r=1; r<=source -> row; r++) {
        for(int c=1; c<=source -> col; c++){
            setNum(dest, getNum(source, r, c), r, c);
        }
    }
    return dest;
}
```
get和set
```cpp
int getPosition(const struct Matrix *now_used, const int row, const int col) {
    return (now_used -> col) * (row - 1) + col - 1;
}

float getNum(const struct Matrix *now_used, const int row, const int col) {
    if( row > ( now_used -> row ) || col > ( now_used -> col ) || row <= 0 || col <= 0 ) return 0.0;
    return now_used -> data[ getPosition (now_used, row, col) ];
}
void setNum(struct Matrix *now_used, const float num, const int row, const int col) {
    if( row > (now_used -> row) || col > (now_used -> col) || row<=0 || col<=0 ) return ;
    now_used -> data[ getPosition (now_used, row, col) ]=num;
}
```
请注意，头文件中并未提供getPosition的定义，因为getPosition需要先**确认位置合法**  
此外，请注意矩阵求最大最小值，我为节省空间，采用了类似strlen的处理，其复杂度为r*c,请尽可能减少这个函数的使用  

## Part 3. Result
我留了一份示例代码在main.c中，您可以通过gcc编译器或我提供的makefile运行main.c
> a=b=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]的4x4矩阵  
> 输出为a, a+b, a-b, a*b, a*2, max(a*b)

## Part 4. Difficulties & Solutions
处理内存泄漏与用户可能提供的奇形怪状指针：本次project采取了特殊的函数返回值，即函数返回了一个传入的指针参数，这同时避免了对指针不熟练的其他开发人员使用二极指针这样 “艰难困苦” 的东西，也能有效处理用户乱传的如NULL指针（毕竟返回NULL就是有问题嘛），此外，对于所以函数中可能改变指针的一切操作，都会提前free/delete,杜绝除了用户作死（例如传了个常量NULL进去，这没救的）外的一切内存泄漏。  
本次project的鲁棒性也高到了，即使传入一个NULL作为结果指针，仍然不会段错误并且输出正确结果，原理是将用户可能遗忘或用混的malloc大量应用在函数中，尽可能的提高了鲁棒性（例如，main.c中的第二组输出）  
虽然，本次project也有不足之处，本来想着显卡加速矩阵乘，但是太忙，就没做。