#include "btools.h"
#include "jobs.h"
#include <getopt.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
struct intersect_parameters{
    char* bed1;
    char* bed2;
    char* outfile;
    uint32_t bin_size;
    uint32_t shift_val;
    bool strandness;
    uint32_t thread_num;
    uint32_t byte_size;
    uint32_t segment_size;
};
static struct intersect_parameters parameters;
struct run_args{
    bed_index* bi=nullptr;
    bed_dict* bd=nullptr;
};
void *get_overlap(void * j_info){
    auto job_info=(struct job_info *) j_info;
    auto input=(char**) job_info->input;
    auto output=(std::pair<bed6_t*, std::vector<bed6_t*>>**) job_info->output;
    auto args=(struct run_args *) job_info->args;
    auto it=new struct bed_index_it;
    for (int i=0; i<job_info->size; ++i) {
        bed6_t* bed=new_bed(input[i]);
        bed6_t* overlap_bed;
        output[i]=new std::pair<bed6_t*, std::vector<bed6_t*>>;
        output[i]->first=bed;
        init_search(it, bed, args->bi, args->bd, parameters.shift_val);
        while((overlap_bed=next(it, parameters.strandness))!=nullptr) output[i]->second.emplace_back(overlap_bed);
    }
    delete it;
    return nullptr;
}
void intersect_run(){
    auto bed_1=read_bed(parameters.bed2);
    auto bed_idx=create_bed_index(bed_1, parameters.shift_val);
    int byte_size=parameters.byte_size;
    int segment_size=parameters.segment_size;
    char *memory=new char[byte_size];
    char *available=memory;
    char **lines=new char*[1];
    uint32_t line_size=1;
    struct run_args args;
    args.bd=bed_1;
    args.bi=bed_idx;
    std::ifstream input;
    input.open(parameters.bed1);
    input.getline(available, 2000);
    while (!input.eof()){
        int n=0;
        while(!input.eof()){
            if (n>=line_size){
                auto old_lines=lines;
                lines=new char*[line_size<<1u];
                memcpy(lines, old_lines, sizeof(char *)*line_size);
                line_size=line_size<<1u;
                delete []old_lines;
            }
            lines[n++]=available;
            available+=strlen(available)+1;
            if (byte_size-(available-memory)<=2000) break;
            input.getline(available, 2000);
        }
        auto ret=(std::pair<bed6_t*, std::vector<bed6_t*>> **) job_distributer((const void **)lines, n, segment_size?segment_size:(int)floor(n/(2*parameters.thread_num)), parameters.thread_num-1,  &get_overlap, (void *) &args);
        for (int i=0; i<n; ++i) {
            for (int j=0; j<ret[i]->second.size(); ++j){
                print_bed(ret[i]->first);
                std::cout<<'\t';
                print_bed(ret[i]->second[j]);
                std::cout<<"\t"<<std::min(ret[i]->first->chromEnd, ret[i]->second[j]->chromEnd)-std::max(ret[i]->first->chromStart, ret[i]->second[j]->chromStart)<<std::endl;
            }
        }
        for (int i=0; i<n; ++i) delete ret[i]->first;
        for (int i=0; i<n; ++i) delete ret[i];
        available=memory;
        input.getline(available, 2000);
    }
    input.close();
    delete []memory;
    delete []lines;
    delete_bed_index(bed_idx);
    delete_bed_dict(bed_1);
}
void intersect_usage(const char* message)
{
    fprintf(stderr, "%s", message);
}
int intersect_main(int argc, char *argv[]){
    if (argc == 1) {intersect_usage("Invalid arguments!\n"); exit(0);}
    const char *shortOptions = "hvsa:b:o:p:B:M:S:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "strand" , no_argument , NULL, 's' },
                    { "bed1" , required_argument , NULL, 'a' },
                    { "bed2" , required_argument , NULL, 'b' },
                    { "out" , required_argument , NULL, 'o' },
                    { "thread", required_argument, NULL, 'p'},
                    { "bin_size", required_argument, NULL, 'B'},
                    { "byte_size", required_argument, NULL, 'M'},
                    { "segment_size", required_argument, NULL, 'S'},
                    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
            };
    char c;
    parameters.strandness=false;
    parameters.bin_size=128;
    parameters.shift_val=7;
    parameters.thread_num=1;
    parameters.byte_size=10000000;
    parameters.segment_size=0;
    while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
        switch(c){
            case 'h':
                intersect_usage("");
                break;
            case 's':
                parameters.strandness=true;
                break;
            case 'a':
                parameters.bed1=optarg;
                break;
            case 'b':
                parameters.bed2=optarg;
                break;
            case 'o':
                parameters.outfile=optarg;
                break;
            case 'p':
                parameters.thread_num=strtol(optarg, NULL, 10);
                break;
            case 'B':
                parameters.bin_size=strtol(optarg, NULL, 10);
                parameters.shift_val=floor(log(parameters.bin_size)/log(2));
                break;
            case 'M':
                parameters.byte_size=std::max(2000, (int)strtol(optarg, NULL, 10));
                break;
            case 'S':
                parameters.segment_size=strtol(optarg, NULL, 10);
                break;
            default:
                intersect_usage("Invalid arguments!\n");
                exit(0);
        }
    if (argc != optind) {intersect_usage("Invalid arguments!\n"); exit(0);}
    if (parameters.bed1==NULL || parameters.bed2==NULL) {intersect_usage("Input file not provided!\n");exit(0);}
    intersect_run();
    return 0;
}