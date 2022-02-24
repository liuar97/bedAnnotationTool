#include <pthread.h>
#include <vector>
struct job_info{
    const void **input;
    void **output;
    int size;
    void *args;
};
void **job_distributer(const void **list, int list_size, int segment_size, int thread_num, void *(*processor)(void*), void *args);