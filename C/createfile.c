#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

int main() {
    int N = 2028;
    char dirName[200];
    sprintf(dirName,"./output/time/%d/",N);
    if (mkdir(dirName, 0700) == -1) {
        perror("Erro ao criar diretório");
        return 1;
    }

    return 0;
}
