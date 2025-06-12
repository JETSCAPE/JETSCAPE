# JETSCAPE Tunes

Currently, the following tunes are available:

- PP [arXiv:1910.05481](https://arxiv.org/abs/1910.05481)

- AA (soft) [arXiv:2011.01430](https://arxiv.org/pdf/2011.01430)

- AA (hard) [arXiv:2204.01163](https://arxiv.org/pdf/2204.01163)

The different tunes can be run by:
```
./runJetscape ../config/publications_config/arXiv_number/jetscape_user_system_arXiv.xml
```


## Instructions for adding XML files

When a new JETSCAPE paper is published, please add the XML file(s) to `publications_config/arXiv_#` and name the corresponding XML files similar to this one from the PP tune: `jetscape_user_PP_1910.05481.xml`.
Specify the system (PP, AA, PA, ...) and put the arXiv number in the name.

If there is more information needed to reproduce the results an additional `README.md` file can be placed in the directory for the publication.

When XML files are added to the repository, an entry should be added in this file.

## Instructions for running tunes in Docker with specific versions of JETSCAPE

The JETSCAPE Collaboration maintains Docker images with fully installed versions of JETSCAPE and X-SCAPE for several releases. The images are available [here](https://hub.docker.com/r/jetscape/jetscape_full) for JETSCAPE and [here](https://hub.docker.com/r/jetscape/xscape_full) for X-SCAPE.

Use the DockerHub tag corresponding to the version you want to run.  For example, to run the PP tune with JETSCAPE 3.7.1, follow these steps:

1) At the Linux command prompt, navigate to the folder where `jetscape_user_PP_1910.05481.xml` is located.
```bash
cd config/publications_config/arXiv_1910.05481
```
2) Update the `<outputFilename>test_out</outputFilename>` line in the XML file to include the path to the host system.
```xml
<outputFilename>/home/jetscape-user/JETSCAPE/host/test_out</outputFilename>
```
3) Run the following command to execute the JETSCAPE simulation.
```bash
docker run --rm -w /home/jetscape-user/JETSCAPE/build --user $(id -u):$(id -g) --entrypoint /home/jetscape-user/JETSCAPE/build/runJetscape -v $(pwd):/home/jetscape-user/JETSCAPE/host jetscape/jetscape_full:beta_v0.11 ../host/jetscape_user_PP_1910.05481.xml
```

### Explanation of the Docker command
* The above command pulls the tagged Docker image if it isn't already downloaded.
* The `--rm` flag removes the container after it finishes running.
* The `-w` flag sets the working directory inside the container to the JETSCAPE build directory.
* The `--user $(id -u):$(id -g)` flag applies the same user and group IDs from the current host user to the container.
* The `-v $(pwd):/home/jetscape-user/JETSCAPE/host` flag mounts the current directory to the specified path inside the container. This allows the container to access the XML input file and write output files to the host.
* The `--entrypoint` flag specifies the entry point for the container, which is set to call the `runJetscape` executable.
* The `../host/jetscape_user_PP_1910.05481.xml` argument specifies the path to the XML input file inside the container, relative to the mounted directory.

Note that it is required to have Docker installed on your system but it is not required to have JETSCAPE installed, as the Docker image contains the full JETSCAPE installation.
