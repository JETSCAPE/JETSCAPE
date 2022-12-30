#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char* argv[]){
    // Creating FIFO FILE
    char *myfifo="/tmp/myfifo";

    if (mkfifo(myfifo,0777)==-1){
        if (errno != EEXIST){
            printf("File already exist, could not create fifo file\n");
            return 1;
        }
    }
    //Opening FIFO FILE With Write Mode
    printf("Opening...\n");
    // int fd=open(myfifo,O_RDWR);
    int fd=open(myfifo,O_WRONLY);
    if (fd==-1){
        return 3;
    }
    // int fd=open(myfifo,O_WRONLY);
    printf("Opened\n");
    //Writing info FIFO FILE in the write Mode
    int x='Z';
    if (write(fd,&x,sizeof(x))==-1){
        return 2;
    }
    printf("Written\n");
    close(fd);
    printf("Closed\n");
    return 0;
}