# JETSCAPE Tunes

Currently, the following tunes are available:

- PP [arXiv:1910.05481](https://arxiv.org/abs/1910.05481)

- AA (soft) [arXiv:2011.01430](https://arxiv.org/pdf/2011.01430)

- AA (hard) [arXiv:2102.11337](https://arxiv.org/pdf/2102.11337)

The different tunes can be run by:
```
./runJetscape ../config/publications_config/arXiv_number/jetscape_user_system_arXiv.xml
```


## Instructions for adding XML files

When a new JETSCAPE paper is published, please add the XML file(s) to `publications_config/arXiv_#` and name the corresponding XML files similar to this one from the PP tune: `jetscape_user_PP_1910.05481.xml`.
Specify the system (PP, AA, PA, ...) and put the arXiv number in the name.

If there is more information needed to reproduce the results an additional `README.md` file can be placed in the directory for the publication.

When XML files are added to the repository, an entry should be added in this file.