#include "cstr.h"
#include <unordered_map>
#include <vector>
struct exon_t{
    char* seq_id;
    char* source;
    char* feature;
    int32_t start;
    int32_t end;
    char* score;
    char strand;
    int32_t phase;
    std::unordered_map<char*, std::vector<char *>, cstrhash, cstrcmp> attributes;
};
struct transcript_t{
    char* seq_id;
    char* source;
    char* feature;
    int32_t start;
    int32_t end;
    char* score;
    char strand;
    int32_t phase;
    std::unordered_map<char*, std::vector<char *>, cstrhash, cstrcmp> attributes;
    std::vector<struct exon_t*> exons;
    std::vector<struct exon_t*> regions;
};
struct gene_t{
    char* seq_id;
    char* source;
    char* feature;
    int32_t start;
    int32_t end;
    char* score;
    char strand;
    int32_t phase;
    std::unordered_map<char*, std::vector<char *>, cstrhash, cstrcmp> attributes;
    std::vector<struct transcript_t*> transcripts;
};