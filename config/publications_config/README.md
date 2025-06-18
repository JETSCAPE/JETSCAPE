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

## Running tunes in Docker with specific versions of JETSCAPE or X-SCAPE

The JETSCAPE Collaboration maintains Docker images with fully installed versions of JETSCAPE and X-SCAPE for the past several releases. The images are available [here](https://hub.docker.com/r/jetscape/jetscape_full) for JETSCAPE and [here](https://hub.docker.com/r/jetscape/xscape_full) for X-SCAPE.

Use the DockerHub repository and tag corresponding to the version you want to run. For example, to run the PP tune with JETSCAPE 3.7.1, follow these steps:


1) Update the `<outputFilename>test_out</outputFilename>` line in the XML file to include the path to the host system.

```xml
<!-- Use this path if running JETSCAPE -->
<outputFilename>/home/jetscape-user/JETSCAPE/host/test_out</outputFilename>
<!-- Use this path if running X-SCAPE -->
<outputFilename>/home/jetscape-user/X-SCAPE/host/test_out</outputFilename>
```

2) From a Linux bash shell, run the `runContainer.sh` script to execute the simulation. Pass the path to the XML tune and the image `<repository>:<tag>` for the version of JETSCAPE or X-SCAPE you want to run. The image will be downloaded if it isn't available locally.
```bash
./runContainer.sh arXiv_1910.05481/jetscape_user_PP_1910.05481.xml jetscape_full:v3.7.1
```

Note that it is required to have either Docker or Apptainer/Singularity installed on your system but it is not required to have JETSCAPE or X-SCAPE installed, as these images contain the full installation.

## Running tunes in Apptainer with specific versions of JETSCAPE or X-SCAPE

Apptainer (formerly Singularity) is especially useful on HPC clusters where Docker is likely unavailable.  If Docker is unavailable, the `runContainer.sh` script will then check for and use Apptainer.

Request an interactive job on the cluster and ensure necessary modules are loaded (for example, `module load apptainer`) before running the `runContainer.sh` script.
