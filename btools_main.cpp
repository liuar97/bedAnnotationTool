#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include "btools.h"
#include "btools_annotate.h"
#include <cstdint>
#include <cmath>
#include <iostream>

void intersect_main(int, char*[]);
void usage(){
    std::cout<<"Error!"<<std::endl;
    exit(0);
}
int main(int argc, char *argv[])
{
    if (argc==1) usage();
    if (strcmp(argv[1], "intersect")==0) intersect_main(argc-1, argv+1);
    if (strcmp(argv[1], "annotate")==0) annotate_main(argc-1, argv+1);
    return(0);
}


