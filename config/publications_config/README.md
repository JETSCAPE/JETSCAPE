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

The JETSCAPE Collaboration maintains Docker images with fully installed versions of JETSCAPE for the past several releases. The images are available [here](https://hub.docker.com/r/jetscape/jetscape_full).

Use the DockerHub tag corresponding to the version you want to run.  For example, to run the PP tune with JETSCAPE 3.7.1, follow these steps:


1) Update the `<outputFilename>test_out</outputFilename>` line in the XML file to include the path to the host system.

```xml
<outputFilename>/home/jetscape-user/JETSCAPE/host/test_out</outputFilename>
```

2) From a Linux or WSL bash shell, run the `runDocker.sh` script to execute the JETSCAPE simulation.  Pass the path the to XML tune and the image tag for the version of JETSCAPE you want to run. The image will be downloaded if it isn't available locally.
```bash
./runDocker.sh arXiv_1910.05481/jetscape_user_PP_1910.05481.xml beta_v0.11
```

Note that it is required to have Docker installed on your system but it is not required to have JETSCAPE installed, as the Docker image contains the full JETSCAPE installation.

## Instructions for running tunes in Apptainer with specific versions of JETSCAPE

Apptainer (formerly Singularity) is especially useful on HPC clusters where Docker is likely unavailable.
