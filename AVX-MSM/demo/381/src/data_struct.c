#include "data_struct.h"

/*入栈，0表示成，非0表示出错*/
int stack_push(StackInfo_st *s,ElementType value)
{
    if(stack_is_full(s))
        return 0;
    /*先增加topOfStack，再赋值*/
    s->topOfStack++;
    s->stack[s->topOfStack] = value;
    return 1;
}

/*出栈*/
int stack_pop(StackInfo_st *s,ElementType *value)
{
    /*首先判断栈是否为空*/
    if(stack_is_empty(s))
    {
        printf("stack_pop error!");
        return 0;
    }
    *value = s->stack[s->topOfStack];
    s->topOfStack--;
    return 1;
}

/*访问栈顶元素*/
int stack_top(StackInfo_st *s,ElementType *value)
{
    /*首先判断栈是否为空*/
    if(stack_is_empty(s))
        return 0;
    *value = s->stack[s->topOfStack];
    return 1;
}

/*判断栈是否已满，满返回1，未满返回0*/
int stack_is_full(StackInfo_st *s)
{
    return s->topOfStack == STACK_SIZE - 1;
}

/*判断栈是否为空，空返回1，非空返回0*/
int stack_is_empty(StackInfo_st *s)
{
    return s->topOfStack ==  - 1;
}


void point52_copy(point52_t r, const point52_t p) {
	for(int i = 0; i < HT_NWORDS; i++)
    {
        r->x[i] = p->x[i];
        r->y[i] = p->y[i];
        r->z[i] = p->z[i];
    }
}

/*circle_queue*/
int queue_in(circle_queue *q, point52_t op, int id) 
{
    if (queue_is_full(q))
    {
        printf("queue is full\n");
        return 0;
    }
    point52_copy(q->pbuf[q->rear].point, op);
    q->pbuf[q->rear].id = id;
    q->rear = (q->rear + 1) % QUEUE_SIZE;
    return 1;
}

int queue_out(circle_queue *q, int num) 
{
    if (queue_size(q) < num)
    {
        printf("queue is empty\n");
        return 0;
    }
    // int id = q->queue[q->front].id;
    q->front = (q->front + num) % QUEUE_SIZE;
    return 1;
}

int queue_size(circle_queue *q) {
    return (q->rear - q->front + QUEUE_SIZE) % QUEUE_SIZE;
}

int queue_is_empty(circle_queue *q) {
    return q->front == q->rear;
}

int queue_is_full(circle_queue *q) {
    return (q->rear + 1) % QUEUE_SIZE == q->front;
}

void pair52_copy(pair52 r, const pair52 p) {
	for(int i = 0; i < HT_NWORDS; i++)
    {
        r->g1->x[i] = p->g1->x[i];
        r->g1->y[i] = p->g1->y[i];
        r->g1->z[i] = p->g1->z[i];
        r->g2->x[0][i] = p->g2->x[0][i];
        r->g2->x[1][i] = p->g2->x[1][i];
        r->g2->y[0][i] = p->g2->y[0][i];
        r->g2->y[1][i] = p->g2->y[1][i];
        r->g2->z[0][i] = p->g2->z[0][i];
        r->g2->z[1][i] = p->g2->z[1][i];
    }
}

int pair_queue_in(pair_circle_queue *q, pair op, int id)
{
    if (pair_queue_is_full(q))
    {
        printf("pair_queue is full\n");
        return 0;
    }
    // void pair52_copy(pair52 r, const pair52 p)
    pair52_copy(q->pbuf[q->rear].point, op);
    q->pbuf[q->rear].id = id;
    q->rear = (q->rear + 1) % QUEUE_SIZE;
    return 1;
}
int pair_queue_out(pair_circle_queue *q, int num)
{
    if (pair_queue_size(q) < num)
    {
        printf("pair_queue is empty\n");
        return 0;
    }
    q->front = (q->front + num) % QUEUE_SIZE;
    return 1;

}
int pair_queue_size(pair_circle_queue *q)
{
    return (q->rear - q->front + QUEUE_SIZE) % QUEUE_SIZE;

}

int pair_queue_is_empty(pair_circle_queue *q)
{
    return q->front == q->rear;
}

int pair_queue_is_full(pair_circle_queue *q)
{
    return (q->rear + 1) % QUEUE_SIZE == q->front;
}