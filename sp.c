/* give me multiplkes of 3 for */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *nn(char *fname, )
{
    size_t fnsz=1+strlen(argv[1]);
    char *tp;
    char insertstr[10]="_filtered";
    char *newn=calloc(fnsz+10, sizeof(char));
    tp=strrchr(argv[1], '.');
    sprintf(newn, "%.*s%s%s", (int)(tp-argv[1]), argv[1], insertstr, tp);

    printf("%s\n", newn); 
    free(newn);
   return 0;
}
int main(int argc, char *argv[])
{
   /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
	if(argc!=2) {
		printf("Error. Pls supply argument (name of file).\n");
		exit(EXIT_FAILURE);
	}
    size_t fnsz=1+strlen(argv[1]);
    char *tp;
    char insertstr[10]="_filtered";
    char *newn=calloc(fnsz+10, sizeof(char));
    tp=strrchr(argv[1], '.');
    sprintf(newn, "%.*s%s%s", (int)(tp-argv[1]), argv[1], insertstr, tp);

    printf("%s\n", newn); 
    free(newn);
   return 0;
}
