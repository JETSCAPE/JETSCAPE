## Using JETSCAPE via Docker

Docker is a software tool that allows one to deploy an application in a portable environment. 
A docker "image" can be created for the application, allowing any user to run a docker "container" from this image.
We have prepared a docker image for the JETSCAPE environment, which allows you to use JETSCAPE on macOS or linux without
installing a long list of pre-reqs or worrying about interference with software you already have installed. Step-by-step instructions are provided below. 

For those unfamiliar with Docker: To illustrate what this will look like, consider the following standard workflow. 
In a terminal on your machine (call it Terminal 1), you will clone JETSCAPE &mdash; this terminal is on your "host" machine &mdash; 
just a standard, typical terminal. In another terminal (call it Terminal 2), you will invoke a command that runs a pre-defined docker container. 
Terminal 2 then lives entirely inside this docker container, completely separated from your host machine. It can *only* access the files that 
are inside that pre-defined docker container &mdash; and not any of the files on your host machine &mdash; unless we explicitly share a 
folder between them. The standard workflow that we envision is the following: You will share the folder containing JETSCAPE between the 
host machine and the docker container. Then, anytime you want to **build** or **run** JETSCAPE, you *must* do it inside the docker container. 
But anytime you want to edit text files (e.g. to construct your own configuration file), or analyze your output files, you can do this from your 
host machine (which we recommend). Simple as that: Depending which action you want to do, perform it either on the host machine, 
or in the docker container, as appropriate &mdash; otherwise it will not work.

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
    mkdir ~/jetscape-docker
    cd ~/jetscape-docker
    git clone https://github.com/JETSCAPE/JETSCAPE.git
    ```

2. Create and start a docker container that contains all of the JETSCAPE pre-reqs: 

    **macOS:**
    ```
    docker run -it -v ~/jetscape-docker:/home/jetscape-user --name myJetscape jetscape/base:v1.3
   ```

    **linux:**
    ```
    docker run -it -v ~/jetscape-docker:/home/jetscape-user --name myJetscape --user $(id -u):$(id -g) jetscape/base:v1.3
    ```

    This is what the `docker run` command does:
    - `docker run` creates and starts a new docker container from a pre-defined image jetscape/base:v1.3, which will be downloaded if necessary.
    - `-it` runs the container with an interactive shell.
    - `-v` mounts a shared folder between your machine (at ~/jetscape-docker) and the container (at /home/jetscape-user), through which you can transfer files to and from the container. You can edit the location of the folder on your machine as you like.
    - `--name` (optional) sets a name for your container, for convenience. Edit it as you like.
    - `--user $(id -u):$(id -g)` (only needed on linux) runs the docker container with the same user permissions as the current user on your machine (since docker uses the same kernel as your host machine, the UIDs are shared). Note that the prompt will display "I have no name!", which is normal.

3. Build JETSCAPE:
    ```
    cd JETSCAPE
    mkdir build
    cd build
    cmake ..
    make -j4     # Builds using 4 cores; adapt as appropriate
    ```

*That's it!* You are now inside the docker container, with JETSCAPE and all of its prequisites installed. 
You can run JETSCAPE executables or re-compile code. Moreover, since we set up the jetscape-docker folder to be shared between your 
host and the docker container, you can do text-editing etc. on your host machine, and then immediately build JETSCAPE in the docker container. 
Output files are also immediately accessible on your host machine for analysis.

Some useful commands:
- To see the containers you have running, and get their ID: `docker container ls` (`-a` to see also stopped containers)
- To stop the container: `docker stop <container>` or `exit`
- To re-start the container: `docker start -ai <container>`
- To put a running container into detatched mode: `Ctrl-p Ctrl-q`, and to re-attach: `docker attach <container>` 
- To delete a container: `docker container rm <container>`
