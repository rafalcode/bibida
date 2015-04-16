/* give me multiplkes of 3 for */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
   /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
	if(argc!=2) {
		printf("Error. Pls supply argument (name of file).\n");
		exit(EXIT_FAILURE);
	}
    int i, h=atoi(argv[1]);
    for(i=1;i<h;++i) 
        printf("%d ", i*3);
    printf("\n"); 
   return 0;
}
