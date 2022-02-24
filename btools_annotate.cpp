#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include "btools.h"
#include "jobs.h"
#include "tools.h"
struct annotate_parameters{
    char* bed1;
    char* bed2;
    char* priority;
    char* outfile;
    uint32_t bin_size;
    uint32_t shift_val;
    bool strandness;
    uint32_t thread_num;
    uint32_t max_line;
    uint32_t segment_size;
};
static struct annotate_parameters parameters;
struct annotation_t{
    bed_dict* raw_annotation;
    bed_dict* transcripts;
    bed_index* transcripts_idx;
    bed_dict* exons;
    bed_dict* features;
    std::unordered_map <bed6_t*, char*>* transcript2tid;
    std::unordered_map <char*, bed6_t*, cstrhash, cstrcmp>* tid2transcript;
    std::unordered_map <char*, int32_t, cstrhash, cstrcmp>* priority;
};
#define TYPE_NONE 0
#define TYPE_PARTIAL 1
#define TYPE_INCLUDE 2
struct annotation_info_t{
    int32_t annotation_type=TYPE_NONE;
    char* attributes=nullptr;
    char* transcript_type=nullptr;
    int32_t priority=0;
    int32_t exon_length=0;
    int32_t intron_length=0;
    int32_t cds_length=0;
    int32_t utr5_length=0;
    int32_t utr3_length=0;
};
static int type_score[3]={0, 1, 2};
bool is_prior(annotation_info_t* first, annotation_info_t* second){
    if (type_score[first->annotation_type]>type_score[second->annotation_type]) return true;
    if (type_score[second->annotation_type]>type_score[first->annotation_type]) return false;
    if (first->exon_length>0 && second->exon_length==0) return true;
    if (first->exon_length==0 && second->exon_length>0) return false;
    if (first->priority>second->priority) return true;
    if (second->priority>first->priority) return false;
    if (first->exon_length>second->exon_length) return true;
    if (second->exon_length>first->exon_length) return false;
    if (first->intron_length>second->intron_length) return true;
    if (second->intron_length>first->intron_length) return false;
    if (first->cds_length>second->cds_length) return true;
    if (second->cds_length>first->cds_length) return false;
    if (first->utr5_length>second->utr5_length) return true;
    if (second->utr5_length>first->utr5_length) return false;
    if (first->utr3_length>second->utr3_length) return true;
    if (second->utr3_length>first->utr3_length) return false;
    return false;
}
#define min(a, b) (((a)>(b))?(b):(a))
#define max(a, b) (((a)>(b))?(a):(b))
#define intersect_length(a, b) (min((a)->chromEnd, (b)->chromEnd)-max((a)->chromStart, (b)->chromStart))

void get_annotation_info(bed6_t* bed, bed6_t* transcript, annotation_t* annotation, annotation_info_t* annotation_info){
    annotation_info->annotation_type=TYPE_NONE;
    annotation_info->attributes=nullptr;
    annotation_info->transcript_type=nullptr;
    annotation_info->priority=0;
    annotation_info->exon_length=0;
    annotation_info->intron_length=0;
    annotation_info->cds_length=0;
    annotation_info->utr5_length=0;
    annotation_info->utr3_length=0;
    if (transcript->chromStart<=bed->chromStart && transcript->chromEnd>=bed->chromEnd) annotation_info->annotation_type=TYPE_INCLUDE;
    else annotation_info->annotation_type=TYPE_PARTIAL;
    annotation_info->attributes=transcript->name;
    char transcript_type[200];
    strextract(transcript->name, 5, '|', transcript_type, 200);
    auto p=annotation->priority->find(transcript_type);
    annotation_info->transcript_type=p->first;
    annotation_info->priority=p->second;

    char *transcript_id = annotation->transcript2tid->at(transcript);
    for (auto it: *(annotation->exons->at(transcript_id))) {
        if (strcmp(it->name, "exon")==0) annotation_info->exon_length+=max(0, intersect_length(bed, it));
        if (strcmp(it->name, "intron")==0) annotation_info->intron_length+=max(0, intersect_length(bed, it));
    }
    if (annotation->features->find(transcript_id)!=annotation->features->end()) {
        for (auto it: *(annotation->features->at(transcript_id))) {
            if (strcmp(it->name, "cds")==0) annotation_info->cds_length+=max(0, intersect_length(bed, it));
            if (strcmp(it->name, "utr5")==0) annotation_info->utr5_length+=max(0, intersect_length(bed, it));
            if (strcmp(it->name, "utr3")==0) annotation_info->utr3_length+=max(0, intersect_length(bed, it));
        }
    }
}
static const char *get_s[4]={".", "intron", "exon", "exon:intron"};
static const char *get_s1[8]={".", "utr3", "utr5", "utr5:utr3", "cds", "cds:utr3", "cds:utr5", "cds:utr5:utr3"};
static const char *get_t[4]={"none", "partial", "include", "error"};
#define get_index(i)(((i->exon_length>0)<<1u)+(i->intron_length>0))
#define get_index1(i)(((i->cds_length>0)<<2u)+((i->utr5_length>0)<<1u)+(i->utr3_length>0))

void print_annotation_info(struct annotation_info_t* info){
    if (info->attributes==nullptr) printf("none\tintergenic\t.\t.");
    else printf("%s\t%s\t%s\t%s", get_t[info->annotation_type], info->attributes,
            get_s[get_index(info)], get_s1[get_index1(info)] );
}
struct annotation_t* read_annotation(const char* file1, const char* file2){
    auto ret=new struct annotation_t;
    auto raw_annotation=read_bed(file1);
    auto transcripts=new bed_dict;
    auto exons=new bed_dict;
    auto features=new bed_dict;
    auto transcript2tid=new std::unordered_map <bed6_t*, char*>;
    auto tid2transcript=new std::unordered_map <char*, bed6_t*, cstrhash, cstrcmp>;
    char type[2000];
    char tid[2000];
    for (auto &it: *raw_annotation){
        for (int i=0; i<it.second->size(); ++i) {
            bed6_t *bed=it.second->at(i);
            if (strextract(bed->name, 6, '|', type, 2000)!=nullptr && strextract(bed->name, 3, '|', tid, 2000)!=nullptr){
                if (strcmp(type, "transcript")==0) {
                    if (transcripts->find(bed->chrom)==transcripts->end()) (*transcripts)[strnew(bed->chrom)]=new bed_vector;
                    transcripts->at(bed->chrom)->emplace_back(bed);
                    char* new_tid=strnew(tid);
                    (*transcript2tid)[bed]=new_tid;
                    (*tid2transcript)[new_tid]=bed;
                }
                if (strcmp(type, "exon")==0) {
                    delete bed->name;
                    bed->name=strnew("exon");
                    if (exons->find(tid)==exons->end()) (*exons)[strnew(tid)]=new bed_vector;
                    exons->at(tid)->emplace_back(bed);
                }
                if (strcmp(type, "intron")==0) {
                    delete bed->name;
                    bed->name=strnew("intron");
                    if (exons->find(tid)==exons->end()) (*exons)[strnew(tid)]=new bed_vector;
                    exons->at(tid)->emplace_back(bed);
                }
                if (strcmp(type, "utr5")==0) {
                    delete bed->name;
                    bed->name=strnew("utr5");
                    if (features->find(tid)==features->end()) (*features)[strnew(tid)]=new bed_vector;
                    features->at(tid)->emplace_back(bed);
                }
                if (strcmp(type, "utr3")==0) {
                    delete bed->name;
                    bed->name=strnew("utr3");
                    if (features->find(tid)==features->end()) (*features)[strnew(tid)]=new bed_vector;
                    features->at(tid)->emplace_back(bed);
                }
                if (strcmp(type, "cds")==0) {
                    delete bed->name;
                    bed->name=strnew("cds");
                    if (features->find(tid)==features->end()) (*features)[strnew(tid)]=new bed_vector;
                    features->at(tid)->emplace_back(bed);
                }
            }

        }
    }
    auto priority=new std::unordered_map <char*, int32_t, cstrhash, cstrcmp>;
    char line[2000];
    char *items[2];
    std::ifstream input(file2);
    input.getline(line, 2000);
    while (!input.eof()){
        strsplit(line, items, 2, ' ');
        if (priority->find(items[0])==priority->end()) (*priority)[strnew(items[0])]=strtol(items[1], nullptr, 10);
        input.getline(line, 2000);
    }
    input.close();
    ret->raw_annotation=raw_annotation;
    ret->transcripts=transcripts;
    ret->transcripts_idx=create_bed_index(transcripts, parameters.shift_val);
    ret->exons=exons;
    ret->features=features;
    ret->transcript2tid=transcript2tid;
    ret->tid2transcript=tid2transcript;
    ret->priority=priority;
    return ret;
}
void delete_annotation(struct annotation_t* annotation){
    delete_bed_dict(annotation->raw_annotation);
    delete_bed_index(annotation->transcripts_idx);
    for (auto it: *annotation->transcripts) delete it.first;
    for (auto it: *annotation->exons) delete it.first;
    for (auto it: *annotation->features) delete it.first;
    for (auto it: *annotation->tid2transcript) delete it.first;
    for (auto it: *annotation->priority) delete it.first;
    delete annotation->transcripts;
    delete annotation->exons;
    delete annotation->features;
    delete annotation->tid2transcript;
    delete annotation->transcript2tid;
    delete annotation->priority;
    delete annotation;
}
struct annotation_info_t* get_best_annotation_info(bed6_t* bed, struct annotation_t* annotation, struct annotation_info_t* ret){
    struct bed_index_it it;
    struct annotation_info_t info;
    bed6_t* transcript;
    init_search(&it, bed, annotation->transcripts_idx, annotation->transcripts, parameters.shift_val);
    while((transcript=next(&it, parameters.strandness))!=nullptr) {
        get_annotation_info(bed, transcript, annotation, &info);
        if (info.priority==0) continue; //skip those transcript with zero priority
        if (is_prior(&info, ret)) *ret=info;
    }
    return ret;
}
void *get_annotations(void * j_info){
    auto job_info=(struct job_info *) j_info;
    auto input=(bed6_t**) job_info->input;
    auto output=(struct annotation_info_t**) job_info->output;
    auto annotation=(struct annotation_t*) job_info->args;
    for (int i=0; i<job_info->size; ++i) {
        bed6_t* bed=input[i];
        output[i]=new annotation_info_t;
        get_best_annotation_info(bed, annotation, output[i]);
    }
    return nullptr;
}
void annotate_run(){
    int segment_size=parameters.segment_size;
    auto annotation=read_annotation(parameters.bed2, parameters.priority);
    std::ifstream input;
    char line[2000];
    auto beds=new bed6_t*[parameters.max_line];
    input.open(parameters.bed1);
    input.getline(line, 2000);
    while (!input.eof()){
        int n=0;
        while(!input.eof()){
            beds[n++]=new_bed(line);
            if (n==parameters.max_line) break;
            input.getline(line, 2000);
        }
        auto ret=(struct annotation_info_t **) job_distributer((const void **)beds, n, segment_size?segment_size:(int)floor(n/(2*parameters.thread_num)), parameters.thread_num-1,  &get_annotations, (void *) annotation);
        for (int i=0; i<n; ++i){
            auto bed=beds[i];
            auto info=ret[i];
            print_bed(bed);
            printf("\t");
            print_annotation_info(info);
            printf("\n");
        }
        for (int i=0; i<n; ++i) delete ret[i];
        delete []ret;
        for (int i=0; i<n; ++i) delete beds[i];
        input.getline(line, 2000);
    }
    input.close();

    delete_annotation(annotation);
}
void annotate_usage(const char* info){
    printf("%s", info);
}
int annotate_main(int argc, char *argv[]){
    if (argc == 1) {annotate_usage("Invalid arguments!\n"); exit(0);}
    const char *shortOptions = "hvim:r:a:b:o:p:B:M:S:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "ignore_strand" , no_argument , NULL, 'i' },
                    { "mode" , required_argument , NULL, 'm' },
                    { "bed1" , required_argument , NULL, 'a' },
                    { "bed2" , required_argument , NULL, 'b' },
                    { "priority" , required_argument , NULL, 'r' },
                    { "out" , required_argument , NULL, 'o' },
                    { "thread", required_argument, NULL, 'p'},
                    { "bin_size", required_argument, NULL, 'B'},
                    { "max_line", required_argument, NULL, 'M'},
                    { "segment_size", required_argument, NULL, 'S'},
                    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
            };
    char c;

    parameters.strandness=true;
    parameters.bin_size=128;
    parameters.shift_val=7;
    parameters.max_line=10000;
    parameters.segment_size=0;
    parameters.thread_num=1;

    while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
        switch(c){
            case 'h':
                annotate_usage("");
                break;
            case 'i':
                parameters.strandness=false;
                break;
            case 'm':
                if (strcmp(optarg, "include_only")==0){
                    type_score[0]=0;
                    type_score[1]=-1;
                    type_score[2]=2;
                }
                if (strcmp(optarg, "equal")==0){
                    type_score[0]=0;
                    type_score[1]=1;
                    type_score[2]=1;
                }
                if (strcmp(optarg, "include_first")==0){
                    type_score[0]=0;
                    type_score[1]=1;
                    type_score[2]=2;
                }
                break;
            case 'a':
                parameters.bed1=optarg;
                break;
            case 'b':
                parameters.bed2=optarg;
                break;
            case 'r':
                parameters.priority=optarg;
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
                parameters.max_line=(int)strtol(optarg, NULL, 10);
                break;
            case 'S':
                parameters.segment_size=strtol(optarg, NULL, 10);
                break;
            default:
                annotate_usage("Invalid arguments!\n");
                exit(0);
        }
    if (argc != optind) {annotate_usage("Invalid arguments!\n"); exit(0);}
    if (parameters.bed1==NULL || parameters.bed2==NULL) {annotate_usage("Input file not provided!\n");exit(0);}
    annotate_run();
    return 0;

}