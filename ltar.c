#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/param.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <libtar.h>
#include <dirent.h>

#define TAR_MAXPATHLEN 256
int use_gnu = 1;

static int tarupthelist(char *tarfile, char *rootdir, libtar_list_t *l)
{
    TAR *t;
    char *pathname;
    char buf[TAR_MAXPATHLEN];

    if (tar_open(&t, tarfile, NULL, O_WRONLY | O_CREAT, 0644, (use_gnu ? TAR_GNU : 0)) == -1) {
        fprintf(stderr, "tar_open(): %s\n", strerror(errno));
        return -1;
    }

    libtar_listptr_t lp;
    libtar_listptr_reset(&lp);
    while (libtar_list_next(l, &lp) != 0) {
        pathname = (char *)libtar_listptr_data(&lp);
        if (pathname[0] != '/' && rootdir != NULL)
            snprintf(buf, sizeof(buf), "%s/%s", rootdir, pathname);
        else
            strncpy(buf, pathname, sizeof(buf));
        if (tar_append_tree(t, buf, pathname) != 0) {
            fprintf(stderr, "tar_append_tree(\"%s\", \"%s\"): %s\n", buf, pathname, strerror(errno));
            tar_close(t);
            return -1;
        }
    }

    if (tar_append_eof(t) != 0) {
        fprintf(stderr, "tar_append_eof(): %s\n", strerror(errno));
        tar_close(t);
        return -1;
    }

    if (tar_close(t) != 0) {
        fprintf(stderr, "tar_close(): %s\n", strerror(errno));
        return -1;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    DIR *dirp;
    struct dirent *direntp;

    if(argc != 2) {
        fprintf( stderr, "Usage: %s tarfilename\n", argv[0]);
        exit(1);
    }
    if ((dirp = opendir(".")) == NULL) { /* yes, we're assuming the CWD */
        perror("open CWD");
        exit(2);
    }
    libtar_list_t *l;

    l = libtar_list_new(LIST_QUEUE, NULL); /* First arg an int, a typey-index of some sort, second arg is "libtar_cmpfunc_t" a fnptr which we set to NULL here */
    while ((direntp = readdir(dirp)) != NULL)
        if(direntp->d_type == DT_REG)
            libtar_list_add(l, direntp->d_name);

    tarupthelist(argv[1], ".", l);

    closedir(dirp);
    return 0;
}
