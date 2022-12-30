
// C program to implement one side of FIFO
// This side writes first, then reads
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
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
    int fd=open(myfifo,O_RDWR);
    // int fd=open(myfifo,O_WRONLY);
    if (fd==-1){
        return 3;
    }

    printf("Opened\n");
    ifstream hepmc_file;
    hepmc_file.open("../../../build/test_out.hepmc");
    string current_line;
    char* current_line_char_array;
    //Writing info FIFO FILE in the write Mode
    if (hepmc_file.is_open())
    {
        cout<<"HepMc File is opened";
        cout<< "Reading hepMc file content"<<endl;
        while ( hepmc_file ) {
            // Take the next line from the hepmc file
            getline (hepmc_file, current_line);
            cout<<current_line<<endl;
            current_line_char_array = &current_line[0];

            // Write the line on FIFO
            if(write(fd, current_line_char_array, strlen(current_line_char_array)+1)==-1){
                return 2;
            }
        }
        hepmc_file.close();
    }
    printf("Written\n");
    // and close it
    close(fd);
    printf("Closed\n");
    return 0;
}