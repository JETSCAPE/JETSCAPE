## Using JETSCAPE via Docker

Docker is a software tool that allows one to deploy an application in a portable environment. 
A docker "image" can be created for the application, allowing any user to run a docker "container" from this image.
We have prepared a docker image for the JETSCAPE environment, which allows you to use JETSCAPE on macOS or linux without
installing a long list of pre-reqs or worrying about interference with software you already have installed.

### Step 1: Install Docker

#### macOS

1. Install Docker Desktop for Mac: https://docs.docker.com/docker-for-mac/install/
2. Open Docker, go to Preferences --> Advanced and 
    - (i) Set CPUs to max that your computer has (`sysctl -n hw.ncpu`),
    - (ii) Set memory to what you are willing to give Docker.

#### linux

1. Install Docker: https://docs.docker.com/install/
2. Allow your user to run docker (requires admin privileges): 
    ```
    sudo groupadd docker
    sudo usermod -aG docker $USER
    ```
    Log out and log back in.

### Step 2: Run JETSCAPE

The docker container will contain only the pre-requisite environment to build JETSCAPE, but will not actually contain JETSCAPE itself. Rather, we will create a directory on our own machine with the JETSCAPE code, and share this directory with the docker container. This will allow us to build and run JETSCAPE inside the docker container, but to easily edit macros and access the output files on our own machine. 

1. Make a directory on your machine (which will be shared with the docker container), and clone JETSCAPE into it. 
    ```
    mkdir ~/jetscape-user
    cd ~/jetscape-user
    git clone https://github.com/JETSCAPE/JETSCAPE.git
    ```

2. Create and start a docker container that contains all of the JETSCAPE pre-reqs: 

    **macOS:**
    ```
    docker run -it -v ~/jetscape-user:/home/jetscape-user --name myJetscape jdmulligan/jetscape-base:v1
    ```
    
    **linux:**
    ```
    docker run -it -v ~/jetscape-user:/home/jetscape-user --name myJetscape --user $(id -u):$(id -g) jdmulligan/jetscape-base:v1
    ```

    This is what the `docker run` command does:
    - `docker run` creates and starts a new docker container from a pre-defined image jdmulligan/jetscape-base:v1, which will be downloaded if necessary.
    - `-it` runs the container with an interactive shell.
    - `-v` mounts a shared folder between your machine (at ~/jetscape-user) and the container (at /home/jetscape-user), through which you can transfer files to and from the container. You can edit the location of the folder on your machine as you like.
    - `--name` (optional) sets a name for your container, for convenience. Edit it as you like.
    - `--user $(id -u):$(id -g)` (only needed on linux) runs the docker container with the same user permissions as the current user on your machine (since docker uses the same kernel as your host machine, the UIDs are shared). Note that the prompt will display "I have no name!", which is normal.

3. Build JETSCAPE as usual:
    ```
    cd JETSCAPE
    mkdir build
    cd build
    cmake ..
    make -j4
    ```

*That's it!* You are now inside the docker container, with JETSCAPE and all of its prequisites installed. You can run JETSCAPE executables or edit and re-compile code. Moreover, since we set up the jetscape-user folder to be shared between your host and the docker container, you can do text-editing etc. on your host machine, and then immediately build JETSCAPE in the docker container. Output files are also immediately accessible on your host machine for analysis.

Some useful commands:
- To see the containers you have running, and get their ID: `docker container ls` (`-a` to see also stopped containers)
- To stop the container: `docker stop <container>` or `exit`
- To re-start the container: `docker start -ai <container>`
- To put a running container into detatched mode: `Ctrl-p Ctrl-q`, and to re-attach: `docker attach <container>` 
- To delete a container: `docker container rm <container>`

========================================================================================

If you have any issues or suggestions, please send a mail to james.mulligan@berkeley.edu.
