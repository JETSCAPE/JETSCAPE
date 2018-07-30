<html>
<head>
<title>A Python Interface</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='PythonInterface.php'>
 
<h2>A Python Interface</h2> 
 
<a href="https://www.python.org/" target="_top">Python</a> is a 
powerful, high-level 
interpreted language that is widely used within the particle physics 
community. It is particularly useful as an interactive language that 
provides fast proto-typing of code. An experimental Python interface 
to PYTHIA is available. This interface has been generated 
using <a href="http://www.swig.org/" target="_top">SWIG</a> 
(Simplified Wrapper and Interface Generator). Please note that this 
interface has not been extensively tested and consequently issues may 
arise. If bugs are found or additional features are required, as 
usual, please report them to the PYTHIA authors. 
 
<h3>Features</h3> 
 
An attempt has been made to translate all PYTHIA classes and functions 
to the Python interface as directly as possible. The following 
features are included in the Python interface: 
<ul> 
<li>The interface is compatible with Python version <code>2.1</code> 
and beyond, including <code>3.5</code>, the current version of Python.</li> 
<li>All PYTHIA classes and functions are 
available. See <code>main01.py</code> for a direct Python translation 
of the C++ <code>main01.cc</code> example.</li> 
<li>Most of the plugin classes are also available in the 
interface. See <code>main34.py</code> for a direct Python translation 
of the C++ <code>main34.cc</code> example which uses 
the <code>LHAupMadgraph</code> class 
from <code>include/Pythia8Plugins/LHAMadgraph.h</code>.</li> 
<li>When available, documentation through the 
built-in <code>help</code> function in Python is provided. Please note 
that this documentation is automatically generated, similar to 
the <a href="http://home.thep.lu.se/~torbjorn/doxygen/index.html">Doxygen</a> 
documentation. Consequently, the inline Python documentation is not a 
replacement for this manual.</li> 
<li>All operators defined in C++, e.g. <code>Vec4*double</code>, as 
well as reverse operators, e.g. <code>double*Vec4</code>, are 
available.</li> 
<li>Classes with defined <code>[]</code> operators are iterable, using 
standard Python iteration, e.g. <code>for prt in 
pythia.event</code>.</li> 
<li>Classes with a <code>&lt&lt</code> operator or a <code>list</code> 
function can be printed via the built-in <code>print</code> function 
in Python. Note this means that a string representation 
via <code>str</code> is also available for these classes in 
Python.</li> 
<li>Specific versions of templates needed by PYTHIA classes are 
available where the naming scheme is the template class name followed 
by its arguments (stripped of namespace specifiers); pointers to 
classes are prepended with <code>Ptr</code>. For 
example, <code>vector&ltint&gt</code> is available via the interface 
as <code>VectorInt</code>, <code>map&ltstring, Mode&gt</code> 
as <code>MapStringMode</code>, 
and <code>vector&ltProcessContainer*&gt</code> 
as <code>VectorProcessContainerPtr</code>.</li> 
<li>Derived classes in Python, for a subset of PYTHIA classes, can be 
passed back to PYTHIA. This is possible for all classes that can be 
passed to the <code>Pythia</code> class via the <code>setXPtr</code> 
functions and includes the following classes: <code>BeamShape</code>, 
<code>DecayHandler</code>, <code>LHAup</code>, 
<code>MergingHooks</code>, <code>PDF</code>, <code>PhaseSpace</code>, 
<code>ResonanceWidths</code>, <code>RndmEngine</code>, 
<code>SigmaProcess</code>, <code>SpaceShower</code>, <code>TimeShower</code>, 
and <code>UserHooks</code>. The protected functions and members of 
these classes are available through the Python 
interface. See <code>main10.py</code> for a direct Python translation 
of the C++ <code>main10.cc</code> example which uses a derived class 
from the <code>UserHooks</code> class to veto events.</li> 
</ul> 
 
<h3>Limitations</h3> 
 
A variety of methods to interface C++ code with Python exist, each 
with its own advantages and disadvantages, some of which are being 
rapidly developed. Currently, for the purposes of an interface to 
PYTHIA, SWIG provides the best option. However, this may not remain 
the case, and the technical interface may be changed to some other 
method, e.g. 
<a href="http://doc.pypy.org/en/latest/cppyy.html" target="_top">cppyy</a>, 
in the future. The Python interface to PYTHIA currently suffers the 
following limitations: 
<ul> 
<li>In the <code>CoupSUSY</code> class all public members that are 
3-by-3 arrays cannot be accessed, these 
include <code>LsddX</code>, <code>LsuuX</code>, <code>LsduX</code>, 
<code>LsudX</code>, <code>LsvvX</code>, <code>LsllX</code>, 
<code>LsvlX</code>, <code>LslvX</code>, as well as the 
equivalent <code>R</code> versions of these 
members. Additionally, <code>rvLLE</code>, <code>rvLQD</code>, 
and <code>rvUDD</code> cannot be accessed.</li> 
<li>In the <code>MergingHooks</code> class, the protected 
methods <code>orderHistories</code>, <code>allowCutonRecState</code>, 
and <code>doWeakClustering</code> with <code>bool</code> return values 
have been renamed 
as <code>getOrderHistories</code>, <code>getAllowCutonRecState</code>, 
and <code>getDoWeakClustering</code>, respectively, in the Python 
interface.</li> 
<li>The public <code>headerStream</code>, <code>initStream</code>, 
and <code>eventStream</code> members of the <code>Writer</code> class, 
used for writing LHEF output, cannot be accessed from the Python 
interface.</li> 
<li>For derived Python classes of the PYTHIA class <code>LHAup</code>, 
the protected member <code>osLHEF</code> cannot be accessed.</li> 
<li>The wrapper generated by SWIG is large (10 MB), and 
consequently the compile time can be significant. The only way to 
reduce the size of the wrapper is to remove functionality from the 
interface.</li> 
<li>Creating a derived Python class from a PYTHIA class, as described 
above in the features, is only possible for a subset of PYTHIA 
classes. However, if this feature is needed for specific classes, they 
can be added in the future upon request. This feature is not enabled 
by default for all classes to reduce the generated wrapper size.</li> 
<li>Python interfaces have not been generated for plugins 
within <code>include/Pythia8Plugins</code> which have direct external 
dependencies. This means there are no Python interfaces for any of the 
classes or functions defined 
in <code>EvtGen.h</code>, <code>FastJet3.h</code>, <code>HepMC2.h</code>, 
or <code>LHAFortran.h</code>. However, interfaces are available for 
all remaining plugins, including both <code>LHAMadgraph.h</code> 
and <code>PowhegProcs.h</code>.</li> 
</ul> 
 
<h3>Installation</h3> 
 
To install the Python interface, both the <code>python</code> command, 
as well as the Python system header <code>Python.h</code> must be 
available. The directory containing the <code>python</code> command 
can be passed to the PYTHIA configuration via the 
option <code>--with-python-bin</code>, while the directory 
containing <code>Python.h</code> can be set with the 
option <code>--with-python-include</code>. An example configuration 
could be as follows, 
<pre> 
    ./configure --with-python-include=/usr/include/python2.7 \ 
    --with-python-bin=/usr/bin 
</pre> where the paths must be changed accordingly for the local 
system. If the location of <code>Python.h</code> is unknown, 
oftentimes the command <code>python-config --includes</code> will 
supply the correct path. Please note that the Python versions for 
the <code>python</code> command and <code>Python.h</code> header must 
match. This is of particular importance when compiling against Python 
3. Many systems will provide the Python 3 command 
via <code>python3</code> rather than <code>python</code>, so either a 
temporary alias should be made, or a soft link of 
the <code>python3</code> command to <code>python</code> could also be 
made. However, take care, as many systems rely on Python 2 for things 
such as package managers, etc. Also note that if one wishes to utilize 
GZIP support (needed for the <code>LHAupMadgraph</code> plugin) then 
the option <code>--with-gzip</code> should also be provided. 
 
<p/> 
After configuring the Python interface for PYTHIA to be built and 
running <code>make</code> as usual, the following files should be 
generated in the directory <code>lib</code>. 
<ul> 
<li><code>pythia8.py</code>: the Python code for the interface.</li> 
<li><code>pythia8.pyc</code>: the byte-compiled Python code for the 
interface.</li> 
<li><code>_pythia8.so</code>: the C++ library needed for the interface.</li> 
<li><code>libpythia8.[so,dylib]</code>: the standard shared PYTHIA 
library.</li> 
</ul> 
 
<p/> 
To ensure that the <code>pythia8.py</code> module is available to 
Python, the system variable <code>PYTHONPATH</code> should be set similar to 
<pre> 
    export PYTHONPATH=$(PREFIX_LIB):$PYTHONPATH 
</pre> 
where <code>PREFIX_LIB</code> is the directory <code>lib</code> which 
contains <code>pythia8.py</code>. Generally, the library paths should 
be set correctly, but it also may be necessary to set 
<pre> 
    export LD_LIBRARY_PATH=$(PREFIX_LIB):$LD_LIBRARY_PATH 
</pre> 
where <code>DYLD</code> should be substituted for <code>LD</code> in 
OS X. Alternatively, it is also possible to define the Python path from 
within Python, as is done within the provided examples. 
 
<h3>Examples</h3> 
 
To use the Python interface for PYTHIA, start Python 
and <code>import pythia8</code>. The provided examples can be run 
by <code>python mainXX.py</code> where <code>XX</code> is the number 
of the example. 
 
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
