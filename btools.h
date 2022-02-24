#include <cstdint>

#include <iostream>
#include <vector>
#include <unordered_map>

#include <string.h>
#include "cstr.h"
#include "tools.h"
struct bed6_t{
    char* chrom;
    int32_t chromStart;
    int32_t chromEnd;
    char* name;
    char* score;
    char strand;
};

typedef std::vector<struct bed6_t*> bed_vector;
typedef std::unordered_map<const char *, bed_vector*, cstrhash, cstrcmp> bed_dict;
typedef std::unordered_map<const char *, std::vector<int32_t> *, cstrhash, cstrcmp> bed_index;

struct bed_index_it{
    bed_vector* bv=nullptr;
    uint32_t idx=0;
    bed6_t* bed=nullptr;
    bool not_found=false;
};

#define isOverlap(a, b) ((std::min((a)->chromEnd, (b)->chromEnd))>(std::max((a)->chromStart,(b)->chromStart)))
#define isSameStrand(a, b) (((a)->strand)==((b)->strand))
struct bed6_t* new_bed(char *);
void delete_bed(bed6_t* );
bool compare_bed(const bed6_t*, const bed6_t*);
void print_bed(bed6_t* );
bed_dict* read_bed(const char* );
void delete_bed_dict(bed_dict *);
bed_index *create_bed_index(bed_dict *, int);
void delete_bed_index(bed_index* );
void init_search(struct bed_index_it* , bed6_t* , bed_index* , bed_dict*, uint32_t);
bed6_t *next(bed_index_it* , bool);

int annotate_main(int , char*[]);


