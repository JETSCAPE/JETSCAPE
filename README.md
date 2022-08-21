# JETSCAPE 3.5.1

The [JETSCAPE](http://jetscape.org) simulation framework is an overarching computational envelope for developing complete event generators for heavy-ion collisions.
It allows for modular incorporation of a wide variety of existing and future software that simulates different aspects of a heavy-ion collision.
For a full introduction to JETSCAPE, please see [The JETSCAPE framework](https://arxiv.org/abs/1903.07706).

Please cite [The JETSCAPE framework](https://arxiv.org/abs/1903.07706) if you use this package for scientific work.

## Installation

Please see the [Installation Instructions](https://github.com/JETSCAPE/JETSCAPE/wiki/Doc.Installation).

## Running JETSCAPE

The main executable to generate JETSCAPE events is `runJetscape`, located in the `build/` directory.
To generate JETSCAPE events, you should pass an XML file specifying the settings with which you would like to run:

```
./runJetscape ../config/jetscape_user.xml
```

### The XML Configuration

All of the JETSCAPE settings are specified by two XML files:
- Main XML file: *you don't modify this*
  - Contains default values for every possible module and parameter
- User XML file: *you provide this*
  - Contains a list of which modules to run, and which default parameter values to override

An example User XML file is provided at `config/jetscape_user.xml`. 
You should adapt this as you like:
 - Set number of events to run
 - Set output format (`Writer` type) and filename 
 - Set which modules to include (in order of execution)
 - Set any default parameter values (from Main XML file) to override
 
The Main XML file is located at `config/jetscape_main.xml`, and contains a list of 
the default parameter settings which will be used for all activated modules (as specified by the User XML file),
if they are not overridden in the User XML file.

You can pass the path to your user XML file as a command-line argument to the `runJetscape` executable:
```
./runJetscape /path/to/my/user_config.xml
```

## JETSCAPE Output

JETSCAPE output can be generated in Ascii, gzipped Ascii, or HepMC format,
and contains a full list of particles and the parton shower history.
You must specify which format you would like to activate in your User XML file.

### Analysis of JETSCAPE Output

Analysis of JETSCAPE output is generally beyond the scope of the JETSCAPE framework, and is the responsibility of the user.
The JETSCAPE docker container includes ROOT, python, fastjet, and several other tools that may be useful.

An example reading an ascii output file is provided:

```bash
./build/readerTest
```

which reads in the generated showers does some DFS search and shows the output. You can generate an output graph format which can be easily visualized using graphViz or other tools like Gephi (GUI for free for Mac) or more adanvanced, graph-tools (Python) and many more. Furthermore, as a "closure" test, the FastJet core package (compiled in our JetScape library) is used to perform a simple jetfinding (on the "final" partons, in graph language, incoming partons in a vertex with no outgoing partons/childs), and since the "shower" is perfectly collinear the jet pT is identical to the hard process parton pT (modulo some random new partons/roots in the final state, see above).  

## JETSCAPE Tunes

Currently, there exists a pp tune [PP19](https://arxiv.org/abs/1910.05481), which can be run by:
```
./runJetscape ../config/jetscape_user_PP19.xml
```

Tuning of Pb-Pb is ongoing.
Several example hydro profiles can be downloaded using `examples/get_hydroSample*`.

## Developing modules

To develop a new JETSCAPE module, you should inherit from the relevant base class (InitialState, JetEnergyLoss, etc.) 
and implement the relevant initialization and execution functions, described in detail in [The JETSCAPE framework](https://arxiv.org/abs/1903.07706)
Section 3.3.

Additionally, you must register your module with the framework with the following steps:
- Add the following to your module .h:
  ```
  private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<MyClass> reg;
  ```
- Add the following to your module .cc: 
  ```
  // Register the module with the base class
  RegisterJetScapeModule<MyClass> MyClass::reg("CustomModuleBlahBlah");
  ```
where `MyClass` is the name of your class, and "CustomModuleBlahBlah" is the name that should be added to the XML configuration.
You can see any of the established modules, e.g.  `Matter`, as an example.

Important Note: In the case of custom modules, you *must* start your module name with "CustomModule..." 
in order for it to be recognized by the framework (for custom writers, you must start the name with "CustomWriter"). 

New modules should not use multiple inheritance, if avoidable.

Once these steps are done, one can just add the module name to the XML, and it will be automatically available to run in JETSCAPE.

## Troubleshooting

If you encounter a problem, please report the issue [here](https://github.com/JETSCAPE/JETSCAPE/issues).
Please be sure to include enough information so that we can reproduce your issue: your platform, JETSCAPE version,
configuration file, and anything else that may be relevant.

## Contributing to JETSCAPE

If you would like to contribute code to JETSCAPE (new module, feature, bug fix, etc.) please open 
a [Pull Request](https://github.com/JETSCAPE/JETSCAPE/pulls) with your changes, or an [Issue](https://github.com/JETSCAPE/JETSCAPE/issues) describing what you intend to do. For further details, see [Tips for git management](https://github.com/JETSCAPE/JETSCAPE/wiki/Tips-for-git-management).
