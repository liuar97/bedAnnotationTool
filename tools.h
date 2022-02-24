#ifndef TOOLS_H
#define TOOLS_H
char ** strsplit(char * line, char ** results, int length,char c='\t');
char * strnew(const char *str);
char *strextract(char*, int, char, char*, int);
#endif
