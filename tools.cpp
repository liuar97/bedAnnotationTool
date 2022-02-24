#include <string.h>
char ** strsplit(char * line, char ** results, int length,char c='\t'){
    char *start=line;
    char *end=nullptr;
    int i=0;
    while ((end=strchr(start, c))!=nullptr && i<length){
        end[0]='\0';
        results[i]=start;
        start=end+1;
        i=i+1;
    }
    if (i<length && start[0]!='\0') {
        results[i]=start;
        i=i+1;
    }
    for (;i<length;++i) results[i]=nullptr;
    return results;
}
char *strextract(char* line, int index, char c,char* ret, int length){
    char* start=line;
    char* end;
    for (int i=0; i<index && start!=nullptr; i++) start=strchr(start+1, c);
    if (start!=nullptr) {
        end=strchr(start+1, c);
        if (end==nullptr) end=strchr(start, '\0');
        if (end!=nullptr && end-start<=length) {
            strncpy(ret, start+1, end-start-1);
            ret[end-start-1]='\0';
            return ret;
        } else return nullptr;
    }
    return nullptr;
}
char * strnew(const char *str){
    char *cstr=new char[strlen(str)+1];
    strcpy(cstr, str);
    return cstr;
}
