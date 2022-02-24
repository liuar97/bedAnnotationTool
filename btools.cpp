#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "jobs.h"
#include "btools.h"
struct bed6_t* new_bed(char * line){
    auto bed=new struct bed6_t;
    char *results[6];
    strsplit(line, results, 6, '\t');
    bed->chrom=strnew(results[0]);
    bed->chromStart=strtol(results[1], nullptr, 10);
    bed->chromEnd=strtol(results[2], nullptr, 10);
    bed->name=strnew(results[3]);
    bed->score=strnew(results[4]);
    bed->strand=results[5][0];
    return bed;
}
void delete_bed(bed6_t* ptr){
    delete ptr->chrom;
    delete ptr->name;
    delete ptr->score;
    delete ptr;
}
bool compare_bed(const bed6_t* bed1, const bed6_t* bed2){
    return bed1->chromStart<bed2->chromStart;
}
void print_bed(bed6_t* bed){
    std::cout<<bed->chrom<<'\t'<<bed->chromStart<<'\t'<<bed->chromEnd<<'\t'<<bed->name<<'\t'<<bed->score<<'\t'<<bed->strand;
}
bed_dict* read_bed(const char* filename){
    auto bd = new bed_dict;
    char line[2000];
    std::ifstream input(filename);
    input.getline(line, 2000);
    while (!input.eof()){
        auto bed=new_bed(line);
        if (bd->find(bed->chrom)==bd->end()) (*bd)[strnew(bed->chrom)]=new bed_vector;
        bd->at(bed->chrom)->emplace_back(bed);
        input.getline(line, 2000);
    }
    input.close();
    for (auto &it: *bd) sort(it.second->begin(), it.second->end(), compare_bed);
    return bd;
}
void delete_bed_dict(bed_dict * bd){
    for (auto const &it: *bd) {
        delete it.first;
        for (auto bed: *it.second) {
            delete_bed(bed);
        }
    }
    delete bd;
}
bed_index *create_bed_index(bed_dict *beds, int shift_val){
    auto bi = new bed_index;
    for (auto const &it: *beds){
        const auto bv=it.second;
        int max_coord=0;
        for (int i=0; i<bv->size(); ++i) max_coord=std::max(max_coord, (*bv)[i]->chromEnd-1);
        int bin_num=((uint32_t) max_coord >> shift_val)+1;
        auto idx=new std::vector <int32_t>;
        idx->resize(bin_num);
        int pos=0;
        for (int i=0; i<bv->size(); ++i){
            int bin=(uint32_t) ((*bv)[i]->chromEnd-1) >> shift_val;
            for ( ; pos<=bin; ++pos) (*idx)[pos]=i;
        }
        (*bi)[strnew(it.first)]=idx;
    }
    return bi;
}
void delete_bed_index(bed_index* bi){
    for (auto const & it: *bi){
        delete it.first;
        delete it.second;
    }
    delete bi;
}
void init_search(struct bed_index_it* it, bed6_t* bed, bed_index* bi, bed_dict* bd, uint32_t shift_val){
    it->not_found=false;
    if (bi->find(bed->chrom)==bi->end()) {it->not_found=true; return;}
    int bin_num=(uint32_t) bed->chromStart >> shift_val;
    auto biv=bi->at(bed->chrom);
    if (bin_num >= biv->size()) {it->not_found=true; return;}
    int32_t idx=biv->at(bin_num);
    auto bdv=bd->at(bed->chrom);
    while (idx<bdv->size() && !isOverlap(bdv->at(idx), bed)) idx++;
    if (idx >= bdv->size()) {it->not_found=true; return;}
    it->bed=bed;
    it->bv=bdv;
    it->idx=idx;
}
bed6_t *next(bed_index_it* it, bool strandness){
    if (it->not_found) return nullptr;
    auto bv=it->bv;
    auto bed=it->bed;
    while (it->idx<bv->size() && (bv->at(it->idx)->chromStart<bed->chromEnd) &&
           !(isOverlap(bv->at(it->idx), bed) && (isSameStrand(bv->at(it->idx), bed) || !strandness))) it->idx++;
    if (it->idx >= bv->size() || (bv->at(it->idx)->chromStart>=bed->chromEnd)) return nullptr;
    return bv->at(it->idx++);
}
