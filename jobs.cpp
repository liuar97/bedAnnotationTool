#include <pthread.h>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "jobs.h"
void **job_distributer(const void **list, int list_size, int segment_size, int thread_num, void *(*processor)(void*), void *args){
    if (thread_num==0){
        auto job_info=new struct job_info;
        auto ret=new void*[list_size];
        job_info->args=args;
        job_info->input=list;
        job_info->output=ret;
        job_info->size=list_size;
        processor((void *) job_info);
        delete job_info;
        return (void**)ret;
    }
    auto threads=new pthread_t[thread_num];
    auto job_info=new struct job_info[thread_num];
    auto available=new bool[thread_num];
    for (int i=0; i<thread_num; ++i){
        available[i]=true;
        job_info[i].args=args;
    }
    auto ret=new void*[list_size];
    int n_submit=0;
    int n_finish=0;
    while (n_finish<list_size) {
        for (int i=0; i<thread_num; ++i) {
            if (available[i] && n_submit < list_size) {
                job_info[i].input=list+n_submit;
                job_info[i].output=ret+n_submit;
                job_info[i].size=std::min(list_size-n_submit, segment_size);
                pthread_create(&threads[i], nullptr, processor, (void *) &job_info[i]);
                available[i]=false;
                n_submit+=job_info[i].size;
            }
        }
        for (int i=0; i<thread_num; ++i) {
            if (!available[i] && pthread_tryjoin_np(threads[i], nullptr)==0) {
                available[i]=true;
                n_finish+=job_info[i].size;
            }
        }
    }
delete []threads;
delete []job_info;
delete []available;
return (void**)ret;
}
