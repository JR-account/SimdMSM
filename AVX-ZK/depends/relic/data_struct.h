#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H
#include "types.h"

#define STACK_SIZE  DBL_P_SIZE /*栈大小*/
#define TOP_OF_STACK -1 /*栈顶位置*/
typedef int ElementType; /*栈元素类型*/

#define QUEUE_SIZE  32
#define SUCCESS 1
#define FAILURE 0

/*define stack*/
typedef struct StackInfo
{
    int topOfStack; /*记录栈顶位置*/
    ElementType stack[STACK_SIZE]; /*栈数组，也可以使用动态数组实现*/
}StackInfo_st;

/*define circle_queue*/
typedef struct 
{
    point_buf pbuf[QUEUE_SIZE];
    int front;
    int rear;
} circle_queue;

typedef struct 
{
    pair_buf pbuf[QUEUE_SIZE];
    int front;
    int rear;
} pair_circle_queue;

/*stack*/
int stack_push(StackInfo_st *s,ElementType value);
int stack_pop(StackInfo_st *s,ElementType *value);
int stack_top(StackInfo_st *s,ElementType *value);
int stack_is_full(StackInfo_st *s);
int stack_is_empty(StackInfo_st *s);

/*circle_queue*/
void point52_copy(point52_t r, const point52_t p);
int queue_in(circle_queue *q, point52_t op, int id);
int queue_out(circle_queue *q, int num);
int queue_size(circle_queue *q);
int queue_is_empty(circle_queue *q);
int queue_is_full(circle_queue *q);

/*pair_circle_queue*/
void pair52_copy(pair52 r, const pair52 p);
int pair_queue_in(pair_circle_queue *q, pair op, int id);
int pair_queue_out(pair_circle_queue *q, int num);
int pair_queue_size(pair_circle_queue *q);
int pair_queue_is_empty(pair_circle_queue *q);
int pair_queue_is_full(pair_circle_queue *q);
#endif