 
<html>
<head>
<title>CKKW-L Merging</title>
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

<form method='post' action='CKKWLMerging.php'>
 
<h2>CKKW-L Merging</h2> 
 
CKKW-L merging [<a href="Bibliography.php#refLon01" target="page">Lon01</a>] allows for a consistent combination 
of tree-level matrix elements containing multiple well-separated partons 
with each other and with parton showers. The result is a calculation that 
contains a mix of processes with different number of well-separated jets 
with fixed-order accuracy, improved by all-order resummation. The 
algorithm implemented  in PYTHIA is described in [<a href="Bibliography.php#refLon11" target="page">Lon11</a>]. To 
perform matrix element merging,  the user has to supply LHE 
files [<a href="Bibliography.php#refAlw07" target="page">Alw07</a>] for the hard process and the corresponding 
process with up to N additional jets. This mix of processes is then 
internally disentangled to ensure that the inclusive fixed-order inputs 
can be converted to exclusive cross sections that no longer overlap. 
Please note that subtleties (and setting scheme) for the EW-improved way of 
disentangling processes presented in [<a href="Bibliography.php#refChr15a" target="page">Chr15a</a>] is discussed in 
the section Electroweak Merging below. 
 
<p/> The usage of the merging procedure is illustrated in a few 
example main  programs 
(<code>main81.cc</code>, <code>main82.cc</code>, 
<code>main83.cc</code>, <code>main84.cc</code> and <code>main85.cc</code>, 
together with the input files <code>main81.cmnd</code>, 
<code>main82.cmnd</code>, <code>main84.cmnd</code> and 
<code>main85.cmnd</code>). These examples should of course only serve 
as  an illustration, and as such will not make use of the merging in 
all  possible ways. For full generality, the example programs link to 
LHAPDF,  FastJet and HepMC. Of course the user is welcome to  remove 
these dependencies. To remove the FastJet dependence, the functions 
calculating example observables have to be deleted. Removing the 
LHAPDF  dependence requires changing the cmnd input files to choose an 
inbuilt PDF,  as outlined in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>PDF documentation</a>.  The 
HepMC dependence can be removed by erasing the code allowing for HepMC 
output. 
 
<p/> Please note that a detailed tutorial on merging in Pythia is available 
from 
<a href="http://home.thep.lu.se/Pythia/pythia8/mergingworksheet8160.pdf"> 
http://home.thep.lu.se/Pythia/pythia8/mergingworksheet8160.pdf</a>. 
 
<p/> Three very short LHE files (<code>w+_production_lhc_0.lhe</code>, 
<code>w+_production_lhc_1.lhe</code>, <code>w+_production_lhc_2.lhe</code>) 
are included in the distribution. These files are not intended for 
physics  studies, but only serve as input for the example main 
programs. For  realistic studies, the user has to supply LHE files. 
 
<p/> In the generation of LHE files, the value of the factorisation 
scale used in  the PDFs is not important, since the cross section will 
be multiplied by ratios  of PDFs to adjust to the PYTHIA starting 
scales. The same is true for the  renormalisation scale (and starting 
value <i>&alpha;<sub>s</sub>(M<sub>Z</sub>)</i>)  used to evaluate 
<i>&alpha;<sub>s</sub></i>. Coupling and scale choices by the user 
will be transferred to the merging routines. 
 
<p/> Multi-jet events can suffer from infrared divergences in the 
calculation. Sensible matrix element generator (MEG) outputs should not 
contain phase space points in which the calculation is ill-defined, meaning 
infrared regions need to be removed by cuts. This is most conveniently done 
by disallowing the MEG to produce partons below a 
minimal parton-parton separation in a certain jet algorithm. Using 
dedicated cuts to regularise MEG output is of course possible as well. Any 
regularisation criterion defines the matrix element region: The parts of 
phase space in which the fixed order calculation is considered valid and 
preferable to the parton shower. Matrix element merging is combining 
MEG events in the matrix element region with parton shower events in regions 
outside the regularisation cut (often called parton shower region). Because 
the regularisation cut defines a boundary between the matrix element 
and parton shower regions, i.e. the regions to be merged into one inclusive 
sample, it is usually called <i> merging scale </i>. Since many different 
cut choices may regularise the MEG calculation, many different merging scale 
definitions are possible. A few standard choices are listed below, as well as 
documentation on how to use a user-defined cut criterion. In combining matrix 
element and parton shower regions, the CKKW-L prescription tries to minimise 
the dependence on the merging scale. This can only be achieved if the 
combination of MEG events and parton shower populates the whole phase space. 
Additional cuts on the partons in the LHEF generation should hence be 
avoided as much as possible, meaning that the merging scale cut should always 
pose a more stringent cut than all other cuts on the partons. Of course, if 
the hard process itself is divergent, cuts need to be made. However, this 
should be chosen in such a way as to not exclude regions that will be 
available to the matrix elements with additional jets. An example is QCD 
di-jet production with additional jets: Say the <i>2 &rarr; 2</i> process is 
regularised with a <i>pTmin</i> cut of pTminCut = 100 GeV, and 
the <i>2 - >3</i> sample is regularised with a <i>kTmin</i>-cut of 
kTminCut = 50 GeV. This would mean that when reclustering 
the  emission in the 2 &rarr; 3 sample, we could end up with a 
<i>pT</i> value <i>pTminNow</i> of the 2 &rarr; 2 configuration with 
<i>pTminCut > pTminNow</i>, which is excluded in the 
2 &rarr; 2 sample. Thus, the 2 &rarr; 3 sample will include a 
Sudakov factor  not included in the 2 &rarr; 2 sample, resulting 
in merging scale  dependencies. Such dependencies can be avoided if 
the additional cuts on the hard process are minimal. 
 
<p/> Of course, additional cuts on electroweak particles are 
allowed. These  should be the same for all samples with 
 <i>0 &lt;= n &lt;= N</i>  additional partons. 
 
<p/> If it is not possible to generate LHE files with minimal cuts, 
the user can choose to use the <code>MergingHooks</code> structures in 
order to decide how much influence to attribute to parton shower 
histories in which the reclustered lowest multiplicity process does 
not pass the matrix element cuts. This is  described below. 
 
<p/> When generating LHE files, we advise against putting 
unstable  particles (e.g. massive gauge bosons) in the final state. 
Rather, specify a  resonance by its decay products, e.g. if Les Houches 
events for the <i>pp &rarr; Z + jets &rarr; e+e- + jets</i> process 
are desired, generate the matrix element events with the <i>Z</i> decay 
included. From a physical  point of view, on-shell final massive gauge 
bosons should not be considered  part of a hard process, since only 
the boson decay products will be detectable.  Furthermore, 
non-narrow-width approximation contributions are not present if  the 
ME generator only produces on-shell bosons. Interference effects 
between  different production channels for the decay products would 
also be neglected.  These points seem an unnecessary restriction on 
the accuracy of the ME  calculation.  In addition, there is a 
technical reason for this strategy. Since  some matrix element 
generators choose to put additional information on  intermediate 
bosons into Les Houches events, depending on if they pass a certain 
criterion (e.g. being close to the mass shell), without exact 
knowledge of this  criterion, the only feasible way of bookkeeping the 
hard process is by  identifying outgoing decay products. 
 
<p/> Despite these considerations, (massive) gauge bosons in the final state 
are allowed in the hard process definition. This is useful particularly for 
Higgs physics, when different decays of the Higgs boson need to be simulated 
after the LHEF generation. 
 
<p/> For all merging purposes, processes with different charge of 
outgoing leptons are considered different processes. That means 
e.g. that <i>e+&nu;<sub>e</sub>+ jets</i> and 
<i>e-&nu;&#772;<sub>e</sub> + jets</i> 
are considered independent processes. If the user wishes to generate 
distributions including effects of more than one  process, merged 
samples for all independent processes should be generated  separately 
and added afterwards. Alternatively, to combine simple processes, 
combined LHE files can be used in conjunction with flexible containers (see 
below). 
 
<p/> When the matrix element merging is used to produce HepMC 
[<a href="Bibliography.php#refDob01" target="page">Dob01</a>] files to be analysed  with RIVET [<a href="Bibliography.php#refBuc10" target="page">Buc10</a>], 
special care  needs to taken in how the cross section is read by RIVET 
(see below). 
 
<p/> To specify the merging conditions, additionally information on 
the merging scale value and the functional definition of the merging 
scale is needed. A few  standard definitions of merging scales are 
available. We hope this makes the user interface less cumbersome. 
 
<p/> Different choices intrinsic to the CKKW-L merging procedure might 
be relevant for the user as well. The base 
class <code>MergingHooks</code> gives the user the opportunity to 
define the functional form of the merging scale.  In the following, 
the usage of the merging machinery to consistently include LHE files 
with additional jets into PYTHIA  will be discussed. 
 
<br/><br/><hr/> 
<h3>Merging scale definitions</h3> 
 
<p/> The quickest way to include processes with additional jets is to 
produce LHE files with one of the standard ways to define the merging 
scale. Three standard ways to define a merging scale (minimal <i>kT</i>, 
minimal evolution <i>pT</i> and by three cuts) are implemented. All of these 
prescriptions are equivalent - different definitions have only been introduced 
for the convenience of users, who might be limited by which cuts can be used 
in the generation of LHE files. Below, we describe how to switch on and use 
these different merging scale definitions. 
 
<h4>Merging with merging scale defined in kT:</h4> 
 
<br/><br/><strong>Merging:doKTMerging</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If the additional jets in the LHE files have been regulated by 
a <i>kT</i> cut, the user can supply the merging scale definition by 
setting this flag to on. <i>kT</i> here and below means cutting on 
Durham <i>kT</i> for <i>e+e-</i>  collisions, and cutting on 
longitudinally invariant <i>kT</i> for hadronic  collisions. Please note 
that this particular merging scale definition will check <i>kT</i> between 
all pairs of <i>u,d,c,s,b,g</i> partons. 
   
 
<p/> 
Currently, the name longitudinally invariant <i>kT</i> is used 
for a few jet recombination algorithms with slightly different 
jet measures. A specific form can be chosen by setting the switch 
 
<br/><br/><table><tr><td><strong>Merging:ktType  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 3</code>)</td></tr></table>
Precise functional definition of longitudinally 
invariant <ei>kT</ei>. For e+e- collisions, <ei>Durham kT</ei> is 
always defined by the square root of <ei>min{ 2*min[ 
E<sub>i</sub><sup>2</sup>, E<sub>j</sub><sup>2</sup>] * [ 1 - 
cos&theta;<sub>ij</sub>] }</ei>, so that this switch will have no effect. 
<br/>
<input type="radio" name="2" value="1" checked="checked"><strong>1 </strong>:  Longitudinally invariant <ei>kT</ei> is defined by  the  square root of the minimum of minimal jet kinematic <ei>pT</ei>  (<ei>p<sub>Tkin,min</sub><sup>2</sup> = min{ p<sub>T,i</sub><sup>2</sup> }   </ei>) and <ei>p<sub>Tlon,min</sub><sup>2</sup> =  min{ min[ p<sub>T,i</sub><sup>2</sup>, p<sub>T,j</sub><sup>2</sup>] *   [ (&Delta;y<sub>ij</sub>)<sup>2</sup> +  (&Delta;&phi;<sub>ij</sub>)<sup>2</sup> ] / D<sup>2</sup> }</ei> , i.e.  <ei>kT = min{ &radic;p<sub>Tkin,min</sub><sup>2</sup>,   &radic;p<sub>Tlon,min</sub><sup>2</sup> }</ei> for hadronic collisions. Note  that the true rapidity of partons is used.  <br/>
<input type="radio" name="2" value="2"><strong>2 </strong>: Longitudinally invariant <ei>kT</ei> is defined by  the  square root of the minimum of minimal jet kinematic <ei>pT</ei>  (<ei>p<sub>Tkin,min</sub><sup>2</sup> = min{ p<sub>T,i</sub><sup>2</sup>   } </ei>) and <ei>p<sub>Tlon,min</sub><sup>2</sup> =  min{ min[ p<sub>T,i</sub><sup>2</sup>,  p<sub>T,j</sub><sup>2</sup>] * [  (&Delta;&eta;<sub>ij</sub>)<sup>2</sup> +  (&Delta;&phi;<sub>ij</sub>)<sup>2</sup> ] / D<sup>2</sup> }</ei>, i.e.  <ei>kT = min{ &radic;p<sub>Tkin,min</sub><sup>2</sup>,   &radic;p<sub>Tlon,min</sub><sup>2</sup> }</ei>  for hadronic collisions. Note that the pseudorapidity of partons is used.  <br/>
<input type="radio" name="2" value="3"><strong>3 </strong>:  Longitudinally invariant <ei>kT</ei> is defined by  the  square root of the minimum of minimal jet kinematic <ei>pT</ei>  (<ei>p<sub>Tkin,min</sub><sup>2</sup> = min{ p<sub>T,i</sub><sup>2</sup>   } </ei>) and  <ei>p<sub>Tlon,min</sub><sup>2</sup> =  min{ min[ p<sub>T,i</sub><sup>2</sup>,  p<sub>T,j</sub><sup>2</sup>] * [ cosh(&Delta;&eta;<sub>ij</sub>) -  cos(&Delta;&phi;<sub>ij</sub>) ] / D<sup>2</sup> } </ei>,  i.e.  <ei>kT = min{ &radic;p<sub>Tkin,min</sub><sup>2</sup>  , &radic;p<sub>Tlon,min</sub><sup>2</sup> }</ei>  for hadronic collisions.  <br/>
 
<br/><br/><table><tr><td><strong>Merging:Dparameter </td><td></td><td> <input type="text" name="3" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
The value of the <i>D</i> parameter needed in the definition of 
longitudinally invariant <i>kT</i> separation. 
   
 
<br/><br/><table><tr><td><strong>Merging:nJetMax  </td><td></td><td> <input type="text" name="4" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Maximal number of additional jets in the matrix element. Note that the 
EW-improved "merging of mergings" strategy presented in [<a href="Bibliography.php#refChr15a" target="page">Chr15a</a>] 
requires a different meaning of "additional", as explained in the 
"Electroweak Merging" section below. 
   
 
<br/><br/><table><tr><td><strong>Merging:TMS </td><td></td><td> <input type="text" name="5" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>)</td></tr></table>
The value of the merging scale. The name is inspired by the scale in 
evolution equations, which is often called 't', and the suffix 'MS' stands 
for merging  scale.  In the particular case of <i>kT</i>-merging, this 
would be the value of the <i>kT</i>-cut  in GeV. For any merging scale 
definition, this input is considered the actual value of the merging 
scale. 
   
 
<br/><br/><table><tr><td><strong>Merging:Process  </td><td></td><td> <input type="text" name="6" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
The string specifying the hard core process, in MG4/ME notation. 
   
 
<p/> 
If e.g. <i>W + jets</i> merging should be performed, set this to 
<code>pp>e+ve</code> (<i>without white spaces or  quotation marks</i>). 
This string may contain resonances in the MG/ME notation, e.g. for merging 
<i>pp&rarr;Z W<sup>+</sup>&rarr;q q&#772; e+&nu;<sub>e</sub> + jets</i>, 
the string <code>pp>(z>jj)(w+>e+ve)</code> would be applicable. 
 
<p/> 
A lot more flexible hard process definitions are possible. To not dwell too 
much on these details here, we will come back to the process string at the end 
of this section. 
 
 
<br/><br/><strong>Merging:doMGMerging</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Even easier, but highly non-general, is to perform the merging with 
MadGraph/MadEvent-produced LHE files, with a merging scale defined by 
a <i>kT</i> cut.  For this, set this switch to on. The merging scale 
value will be read from  the +1 jet LHE file by searching for the 
string <code>ktdurham</code>, and  extracting the value from <code> 
value  = ktdurham</code>. Also, the hard  process will be read from 
the +0 jet LHE file, from the line containing  the string <code>@1</code> 
(the tag specifying the first process in the  MadGraph process card). 
For this to work, PYTHIA should be initialised on LHE files called 
<code>NameOfYourLesHouchesFile_0.lhe</code> (+0 jet sample) and 
<code>NameOfYourLesHouchesFile_1.lhe</code> (+1 jet sample) and the 
same naming convention for LHE files with two or more additional jets. 
Since for this option, the merging scale value is read from the 
LHEF, no merging scale value needs to be supplied by setting <code> 
Merging:TMS </code>.  Also, the hard process is read from LHEF, the 
input <code>Merging::Process</code> does not have to be defined. 
However, the maximal number of merged jets still has to be supplied by 
setting <code>Merging:nJetMax</code>. 
   
 
 
<h4>Merging with merging scale defined in Pythia evolution <i>pT</i>:</h4> 
 
If the LHE file has been regularised by cutting on the minimal Pythia 
evolution <i>pT</i> in the state, this can also be used as a merging scale 
right away. For this, change the switch 
 
<br/><br/><strong>Merging:doPTLundMerging</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The merging scale is then defined by finding the minimal Pythia evolution 
<i>pT</i> between sets of radiator, emitted and recoiler partons. For this 
particular merging scale definition, <i>u,d,c,s,b,g</i> are considered 
partons. The Pythia evolution <i>pT</i> of a single three-parton set is 
defined by 
<p/> 
<i>pT<sub>evol</sub> = z<sub>ijk</sub>(1-z<sub>ijk</sub>) 
   Q<sub>ij</sub><sup>2</sup></i> for FSR, where <i>i</i> is the radiating 
   parton, <i>j</i> is the emitted parton and <i>k</i> is the recoiler, 
   and 
   <i> Q<sub>ij</sub><sup>2</sup> = 
      (p<sub>i</sub> + p<sub>j</sub>)<sup>2</sup> </i>, and 
   <i>z<sub>ijk</sub> = 
      x<sub>i,jk</sub> / (x<sub>i,jk</sub> + x<sub>j,ik</sub>)</i> with 
   <i>x<sub>i,jk</sub> = 
      2 p<sub>i</sub> (p<sub>i</sub> + p<sub>j</sub> + p<sub>k</sub>) 
        / (p<sub>i</sub> + p<sub>j</sub> + p<sub>k</sub>)<sup>2</sup> </i> 
<p/> 
<i>pT<sub>evol</sub> = (1-z<sub>ijk</sub>) 
   Q<sub>ij</sub><sup>2</sup></i> for ISR, where <i>i</i> is the radiating 
   parton, <i>j</i> is the emitted parton and <i>k</i> is the second 
   initial state parton, and 
   <i> Q<sub>ij</sub><sup>2</sup> = 
     -(p<sub>i</sub> - p<sub>j</sub>)<sup>2</sup> </i>, and 
   <i>z<sub>ijk</sub> = 
     (p<sub>i</sub> - p<sub>j</sub> + p<sub>k</sub>)<sup>2</sup> 
   / (p<sub>i</sub> + p<sub>k</sub>)<sup>2</sup> </i>. 
 
<p/> 
When using this option, the merging scale is defined by the minimum 
<i>pT<sub>evol</sub></i> for all combinations of three partons in the event, 
irrespective of flavour or colour-connections. The merging scale value will 
be read from the <code>Merging:TMS</code> parameter, so that this needs to be 
set just as in the case of the <i>kT</i>-merging prescription. Of course you 
will also need to set <code>Merging:Process</code> and the maximal number of 
additional matrix element jets <code>Merging:nJetMax</code>. 
   
 
<h4>Merging with merging scale defined by a combination of cuts:</h4> 
 
It is possible to regularise QCD divergences in a LHE file by applying cuts 
to the kinematical <i>pT</i> of jets (<i>pT<sub>i</sub></i>), combined 
with a cut on <i> &Delta;R<sub>ij</sub></i> between jets and a cut on 
invariant mass <i> Q<sub>ij</sub></i> of jet pairs. The combination of 
these standard cuts can also serve as a merging scale. For this, use this 
setting 
 
<br/><br/><strong>Merging:doCutBasedMerging</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
This switch will use cuts on (<i>pT<sub>i</sub></i>), 
<i> &Delta;R<sub>ij</sub></i>  and 
<i> Q<sub>ij</sub></i> to define when parton shower emissions are allowed. 
Please note for this particular merging scale definition, only light jets 
(<i>u,d,c,s,g</i>) are checked. 
   
 
<p/> 
The values of the cuts will then be read from 
 
<br/><br/><table><tr><td><strong>Merging:QijMS </td><td></td><td> <input type="text" name="10" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>)</td></tr></table>
The value of the invariant mass cut <i> Q<sub>ij</sub></i> of pairs of final 
state partons used in the matrix element generation. 
   
 
<br/><br/><table><tr><td><strong>Merging:pTiMS </td><td></td><td> <input type="text" name="11" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>)</td></tr></table>
The value of the minimal transverse momentum cut <i> pT<sub>i</sub></i> on 
final state partons, as used in the matrix element generation. 
   
 
<br/><br/><table><tr><td><strong>Merging:dRijMS </td><td></td><td> <input type="text" name="12" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>)</td></tr></table>
The value of the minimal <i> &Delta;R<sub>ij</sub></i> separation between 
pairs of final state partons used in the matrix element generation, where 
<i> &Delta;R<sub>ij</sub><sup>2</sup> = (&Delta;y<sub>ij</sub>)<sup>2</sup> + 
(&Delta;&phi;<sub>ij</sub>)<sup>2</sup></i>. 
   
 
<p/> 
With knowledge of these values, and <code>Merging:doCutBasedMerging</code>, 
Pythia will use these cuts as a separation between matrix element phase space 
and parton shower region. If e.g. the Les Houches Events have been generated 
with the cuts <i> &Delta;R<sub>ij</sub> = 0.1 </i>, 
<i>pT<sub>i</sub>= 20 GeV</i> and <i> Q<sub>ij</sub> = 40 GeV</i>, set 
<code>Merging:QijMS=40.</code>, 
<code>Merging:pTjMS=20.</code>, 
<code>Merging:dRijMS=0.1</code> to perform a cut-based merging. Of course 
you will also need to set <code>Merging:Process</code> and the 
maximal number of additional matrix element jets 
<code>Merging:nJetMax</code>. 
 
<h4>Les Houches events outside the matrix element region</h4> 
 
<p/> 
Before continuing, we would like to point out that in general, the user 
should make sure that the events in the Les Houches file are actually 
calculated using the regularisation cut definition and value(s) supplied to 
Pythia as merging scale definition and value(s). However, if LH files with 
a large number of events and loose merging scale cuts are available, it 
might be useful to choose a higher merging scale value, e.g. for merging 
scale variations as part of uncertainty assessments. If CKKW-L merging is 
enabled, Pythia will by default check if events read from Les Houches file 
are in the matrix element region as defined by the merging scale definition 
and merging scale value. Events outside the matrix element region will be 
discarded. This will lead to warnings of the form "<code>Les Houches Event 
fails merging scale cut. Cut by rejecting event</code>". These warnings 
should, in this case, rather be regarded as information. 
To change the default behaviour, use the flag 
 
<br/><br/><strong>Merging:enforceCutOnLHE</strong>  <input type="radio" name="13" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="13" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This will check if the events read from LHE file are in the matrix element 
region as defined by the merging scale definition and value(s). If on, LHE 
input outside the matrix element region will be rejected. If off, every 
event is assumed to pass the merging scale cut. 
   
 
<h4>Defining the hard process</h4> 
 
To perform CKKW-L matrix element merging, the user has to decide on a hard 
process, and communicate this choice to Pythia. This is done by setting the 
input <code>Merging:Process</code>. Note that the EW-improved 
"merging of mergings" strategy presented in [<a href="Bibliography.php#refChr15a" target="page">Chr15a</a>] requires 
a very loose process definition that is discussed in the section Electroweak 
Merging below. 
 
<p/> 
For single processes in the Standard Model or the MSSM, MG4/ME notation is 
applicable. However, for some purposes, using a single simple process 
string is not satisfactory. Mixed W<sup>+</sup> and W<sup>-</sup> events 
in a single LHE file is a common example. For this case, it would of course 
be perfectly allowed to perform twice, once for W<sup>+</sup> production and 
once for W<sup>-</sup> production, and then add the results. Nevertheless, it 
seems reasonable to alleviate difficulties by allowing for less restrictive 
hard process definitions. Two generalisations of the process tag are 
available: Containers and user-defined particle tags. The syntax of these 
settings is described below. 
 
<p/> 
In case you want multiple processes in a LHE file to be treated on equal 
footing (e.g. <i>W<sup>+</sup> + jets</i> and 
<i>W<sup>-</sup> + jets</i>), you should use flexible containers do 
specify the hard process. So far, we allow the use of the containers 
<code>LEPTONS</code>, <code>NEUTRINOS</code>, <code>BQUARKS</code>. If you 
use these containers, the hard process definition is relatively flexible, 
meaning that Pythia will attempt a merging of QCD jets for each event in 
the LHE file, and assume that all particles matching one of the containers 
are products of the hard process. This is best explained by examples. 
If you want to have both <i>pp &rarr; e+ &nu;<sub>e</sub> + jets</i> 
and <i>pp &rarr; e- &nu;&#772;<sub>e</sub> + jets</i> events in a single 
file, you can set <code>Merging:Process=pp>LEPTONS,NEUTRINOS</code> as hard 
process (note that for readability, you are allowed to use commata to separate 
container names). Combining e.g. the processes 
<i>pp &rarr; e+ &nu;<sub>e</sub></i> and 
<i>pp &rarr; &mu;+ &nu;<sub>&mu;</sub></i> is possible with the hard process 
definition <code>pp>LEPTONS,NEUTRINOS</code>. 
 
<p/> 
For maximal flexibility, the definition of the hard process by these 
containers does not mean that each Les Houches event needs to contain 
particles to match each container. It is enough if one container is matched. 
This means that with the string <code>pp>LEPTONS,NEUTRINOS</code>, you can 
immediately process <i>pp &rarr; e+ e- </i> events mixed with 
<i>pp &rarr; e+ &nu;<sub>e</sub></i> events, since particles matching at 
least one container can be found in both cases. Another example for the usage 
of containers is mixing <i>pp &rarr; e+ &nu;<sub>e</sub></i> and 
<i>pp &rarr; tt&#772; &rarr; e+ &nu;<sub>e</sub> e- &nu;&#772;<sub>e</sub> 
  bb&#772;</i>. This can be accommodated by the hard process string 
<code>Merging:Process=pp>LEPTONS,NEUTRINOS,BQUARKS</code>. 
 
<p/> 
There is however a conceptual limitation to containers: The hard process 
definition is necessary to ensure that when constructing lower multiplicity 
states (that will be used to calculate the correct merging weight), the 
structure of the hard process will be preserved. If e.g. we want the hard 
process to be <i>pp &rarr; Z &rarr; bb&#772; </i>, we should ensure that 
the lowest multiplicity state contains a colour-singlet bb&#772;-pair. When 
reconstructing intermediate lower multiplicity states from multi-jet matrix 
elements, we should thus always be able to find at least one bb&#772;-pair. By 
mixing different processes in a LHE file, this requirement might already be 
violated at the level of Les Houches events. Flexible containers cannot give 
strong conditions which flavours should be preserved in the construction of 
the hard process. In order to avoid non-sensible results, it is hence <i> 
assumed that all</i> particles matching any of the containers will be part 
of the lowest multiplicity process. This implies that if you decide to use the 
<code>BQUARKS</code> container, all b-quarks in the LHE file will be 
interpreted as hard process particles, and never as additional radiation. 
 
<p/> 
Another way to specify the hard process particles is to explicitly define the 
particle names and identifiers. This is necessary if the matrix element 
merging in Pythia does not contain the particles of interest. To make sure 
that the hard process is still treated correctly, it is possible to define 
particles in the process string. If you e.g. want the hard process to contain 
a particle "zeta~" with PDG identifier "12345", produced in proton collisions, 
you have to include a user-defined particle tag by setting the process string 
to <code>pp>{zeta~,12345}</code>. The  user-defined particle is enclosed in 
curly brackets, with syntax 
<code>{particle_name,particle_identifier}</code>, where "particle_name" 
and "particle_identifier" are the particle name and particle identifier used 
for this particle in the input LHE file. User-defined particles are only 
allowed in the final state. You are free to fix user-defined particles with 
more common ones, as long as user-defined particles are put before more common 
particles in the process string. This means that if you e.g. wanted the hard 
process to contain a graviton in association with a positron and an 
electron-neutrino, you have to define the hard process as 
<code>pp>{G,39}e+ve</code>. 
 
<p/> 
Below you can find a list of particles predefined in the merging. If you wish 
to include a hard process with different final state particles, you may use 
the "curly bracket notation" outlined above. 
 
<p/> 
The set of incoming particles us limited to: 
<code>e-</code> (electron), <code>e+</code> (positron), <code>mu-</code> 
(muon), <code>mu+</code> (antimuon), <code>p</code> (proton, container to 
hold all initial state coloured particles),  <code>p~</code> (identical to 
<code>p</code> container). 
 
<p/> 
The following intermediate particles are allowed: 
<code>a</code> (photon), <code>z</code> (Z boson), 
<code>w-</code> (W<sup>-</sup> boson), <code>w+</code> (W<sup>+</sup> boson), 
<code>h</code> (scalar Higgs boson), <code>W</code> (container to hold both 
W<sup>-</sup> and W<sup>+</sup> boson), <code>t</code> (top quark), 
<code>t~</code> (anti-top), 
<code>dl</code>, <code>dl~</code>, <code>ul</code>, <code>ul~</code>, 
<code>sl</code>, <code>sl~</code>, <code>cl</code>, <code>cl~</code>, 
<code>b1</code>, <code>b1~</code>, <code>t1</code>, <code>t1~</code>, 
<code>dr</code>, <code>dr~</code>, <code>ur</code>, <code>ur~</code>, 
<code>sr</code>, <code>sr~</code>, <code>cr</code>, <code>cr~</code>, 
<code>b2</code>, <code>b2~</code>, <code>t2</code>, <code>t2~</code> 
(all MSSM squarks). 
 
<p/> 
We have pre-defined the outgoing particles: 
<code>e+</code>, <code>e-</code>, <code>ve~</code>, 
<code>ve</code>, <code>mu+</code>, <code>mu-</code>, 
<code>vm~</code>, <code>vm</code>, <code>ta+</code>, <code>ta-</code>, 
<code>vt~</code>, <code>vt</code> (all SM leptons and neutrinos), 
<code>j~</code> (container to hold all final state coloured particles), 
<code>j</code> (container to hold all final state coloured particles), 
<code>NEUTRINOS</code> (container to hold all final state neutrinos and 
anti-neutrinos), <code>LEPTONS</code> (container to hold all final state 
leptons and anti-leptons), <code>BQUARKS</code> (container to hold final 
state b-quarks), <code>d~</code>, <code>d</code>, <code>u~</code>, 
<code>u</code>, <code>s~</code>, <code>s</code>, <code>c~</code>, 
<code>c</code>, <code>b~</code>, <code>b</code>, <code>t~</code>, 
<code>t</code> (all SM quarks), <code>a</code>, <code>z</code>, 
<code>w-</code>, <code>w+</code> (all SM electro-weak bosons), 
<code>h</code> (scalar Higgs boson), <code>W</code> (container to hold both 
W<sup>-</sup> and W<sup>+</sup> boson), <code>n1</code> (MSSM neutralino), 
<code>dl~</code>, <code>dl</code>, <code>ul~</code>, <code>ul</code>, 
<code>sl~</code>, <code>sl</code>, <code>cl~</code>, <code>cl</code>, 
<code>b1~</code>, <code>b1</code>, <code>t1~</code>, <code>t1</code>, 
<code>dr~</code>, <code>dr</code>, <code>ur~</code>, <code>ur</code>, 
<code>sr~</code>, <code>sr</code>, <code>cr~</code>, <code>cr</code>, 
<code>b2~</code>, <code>b2</code>, <code>t2~</code>, <code>t2</code> 
(all MSSM squarks). Other outgoing particles are possible if you use the 
"curly bracket notation" described earlier. 
 
<br/><br/><hr/> 
<h3>Histogramming the events</h3> 
After the event has been processed, histograms for observables of interest 
need to be filled. In order to achieve good statistical accuracy for all jet 
multiplicities and all subprocesses contributing to one jet multiplicity, 
generally a fixed number of unit-weighted events is read from each Les 
Houches Event file. To then arrive at the correct prediction, for each of 
these events, histogram bins should be filled with the corresponding cross 
section, or weighted with unit weight and normalised at the end to 
the generated cross section for each jet multiplicity separately. 
 
<p/> Still another, even more important, event weight that has to 
applied on an  event-by-event basis is the CKKW-L-weight. This 
corrective weight is the main  outcome of the merging procedure and 
includes the correct no-emission  probabilities, PDF weights and 
coupling (&alpha;<sub>s</sub> or &alpha;<sub>em</sub>) factors. This means 
that the merging implementation will generate weighted events. The 
CKKW-L-weight can be accessed by the following function: 
 
<p/><strong> double Info::mergingWeight() &nbsp;</strong> <br/> 
Returns the CKKW-L weight for the current event. 
 
<p/> Note that to avoid confusion, this function does not include the 
the weight of a phase space point (given 
by <strong>Info::weight()</strong>). This weight will differ from 
unity when reading in weighted Les Houches events. In this case, the 
full weight with which to fill histogram bins is 
<strong>Info::mergingWeight() * Info::weight()</strong>. 
 
<p/> Finally, to arrive at a correct relative normalisation of the 
contributions from different number of additional jets in the matrix 
element, each histogram should be rescaled with the accepted cross 
section given by 
<strong>Info::sigmaGen()</strong>. The accepted cross section includes 
the  effect of vetoes generating Sudakov form factors for the matrix 
elements, and  is in general only known after the run. 
 
<p/> This final step can of course be skipped if the accepted cross 
section had been estimated before the histogramming run, and  histogram 
bins had instead been filled with the weight 
<strong>Info::mergingWeight() * &sigma;<sub>est</sub>(number of 
additional jets in current ME sample)</strong>. This is the way HepMC 
events should be weighted to produce correct relative weights of 
events (see below, and particularly examine the example programs 
<code>main84.cc</code> and <code>main85.cc</code>). 
 
<p/> Examples how to use these options are given in <code>main81.cc</code> 
(<i>kT</i> merging), <code>main84.cc</code> (automatic MG/ME merging 
for RIVET usage), and <code>main85.cc</code> (HepMC output for RIVET usage). 
 
<br/><br/><hr/> 
<h3>Merging with user-defined merging scale function</h3> 
 
<p/> For all other merging scale definitions, the procedure is 
slightly more  complicated, since the user has to write a small piece 
of code defining the  merging scale. To allow for a user defined 
procedure, set the input 
 
<br/><br/><strong>Merging:doUserMerging</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
General user defined merging on/off. 
   
 
<p/> 
Then, set 
the <strong>Merging:nJetMax</strong>, <strong>Merging:TMS</strong> 
and <strong>Merging:Process</strong> input as before. 
 
<p/> Since during execution, PYTHIA needs to evaluate the merging 
scale with the  definition of the user, the user interface is designed 
in a way similar to the 
<code>UserHooks</code> strategy. The class controlling the merging 
scale  definition is called <code>MergingHooks</code>. 
 
<h4>Initialisation</h4> 
 
<p/> To initialise the merging with user-defined merging scale, we 
should construct a class derived from <code>MergingHooks</code>, with 
a constructor and destructor 
 
<p/> 
<a name="anchor1"></a>
<p/><strong>MergingHooks::MergingHooks() &nbsp;</strong> <br/>
   
<a name="anchor2"></a>
<p/><strong>virtual MergingHooks::~MergingHooks() &nbsp;</strong> <br/>
   
The constructor and destructor do not need to do anything. 
 
<p/> For the class to be called during execution, a pointer to an 
object of the class should be handed in with the 
<br/><code><?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?> 
Pythia::setMergingHooksPtr( MergingHooks*)</a></code> method. 
 
An examples of this procedure are given in <code>main82.cc</code>. 
 
<h4>Defining a merging scale</h4> 
 
<p/> Then, in the spirit of the <code>UserHooks</code> class, the user 
needs to  supply the process to be merged by defining a methods to 
evaluate the merging scale variable. 
 
<a name="anchor3"></a>
<p/><strong>virtual double MergingHooks::tmsDefinition(const Event& event) &nbsp;</strong> <br/>
This method will have to calculate the value of the merging scale 
defined in  some variable from the input event record. An example of 
such a function is  given in <code>main82.cc</code>. 
   
 
<p/> The base class <code>MergingHooks</code> contains many functions 
giving  information on the hard process, to make the definition of the 
merging scale as easy as possible: 
 
<a name="anchor4"></a>
<p/><strong>int MergingHooks::nMaxJets() &nbsp;</strong> <br/>
Return the maximum number of additional jets to be merged. 
   
 
<a name="anchor5"></a>
<p/><strong>int MergingHooks::nHardOutPartons() &nbsp;</strong> <br/>
Returns the number of outgoing partons in the hard core process. 
   
 
<a name="anchor6"></a>
<p/><strong>int MergingHooks::nHardOutLeptons() &nbsp;</strong> <br/>
Returns the number of outgoing leptons in the hard core process. 
   
 
<a name="anchor7"></a>
<p/><strong>int MergingHooks::nHardInPartons() &nbsp;</strong> <br/>
Returns the number of incoming partons in the hard core process. 
   
 
<a name="anchor8"></a>
<p/><strong>int MergingHooks::nHardInLeptons() &nbsp;</strong> <br/>
Returns the number of incoming leptons in the hard core process. 
   
 
<a name="anchor9"></a>
<p/><strong>int MergingHooks::nResInCurrent() &nbsp;</strong> <br/>
The number of resonances in the hard process reconstructed from the 
current event. If e.g. the ME configuration was 
<i>pp &rarr; (w+&rarr;e+ve)(z &rarr; mu+mu-)jj</i>, and the ME generator put 
both intermediate bosons into the LHE file, this will return 2. 
   
 
<a name="anchor10"></a>
<p/><strong>double MergingHooks::tms() &nbsp;</strong> <br/>
 Returns the value used as the merging scale. 
   
 
<p/> Filling output histograms for the event then proceeds along the 
lines described above in "Histogramming the events". 
 
<p/> The full procedure is outlined in <code>main82.cc</code>. Special 
care needs to be  taken when the output is stored in the form of HepMC 
files for RIVET usage. 
 
<h4>Defining a cut on lowest jet multiplicity events</h4> 
 
<p/> It can sometimes happen that when generating LHE files, a fairly 
restrictive cut has been used when generating the lowest multiplicity 
matrix element  configurations. Then, it can happen that states that 
are (in the generation of a parton shower history) constructed by 
reclustering from higher multiplicity  configurations, do not pass 
this matrix element cut. 
 
<p/> Consider as an example  pure QCD dijet merging, when up to one 
additional jet should be merged.  Three-jet matrix element 
configurations for which the reclustered two-jet state does not pass 
the cuts applied to the two-jet matrix element would never have  been 
produced by showering the two-jet matrix element. This means that the 
three-jet matrix element includes regions of phase space that would 
never have  been populated by the parton shower. Thus, since the 
matrix element phase space is larger than the shower phase space, 
merging scale dependencies are expected.  A priori, this is not 
troublesome, since the aim of matrix element merging is  to include 
regions of phase space outside the range of the parton shower 
approximation into the shower. An example is the inclusion of 
configurations  with only unordered histories. 
 
<p/> Clearly, if the parton shower phase space is very constrained by 
applying  stringent cuts to the two-jet matrix element, merging scale 
dependencies can  become sizable, as was e.g. seen in [<a href="Bibliography.php#refLon11" target="page">Lon11</a>] 
when forcing shower emissions to be ordered both in the evolution 
variable and in rapidity. To  influence the effect of large phase 
space differences for shower emissions and matrix element 
configurations due to LHEF generation cuts, the user has to  write a 
small piece of code overwriting method 
 
<a name="anchor11"></a>
<p/><strong>virtual double MergingHooks::dampenIfFailCuts(const Event& event) &nbsp;</strong> <br/>
multiplicity  reclustered state as an input Event. From this input 
event, the user can then check if matrix element cuts are 
fulfilled. The return value will be internally multiplied to the 
CKKW-L weight of the current event. Thus, if the user wishes  to 
suppress contributions not passing particular cuts, a number smaller 
than  unity can be returned. 
   
 
<p/> Note that this method gives the user access to the lowest 
multiplicity state,  which ( e.g. in the case of incomplete histories) 
does not have to be a <i>2 &rarr; 2</i> configuration. Also, changing the 
weight of the current event by  hand is of course a major intervention 
in the algorithm, and should be  considered very carefully. Generally, 
if this facility would have to be used extensively, it is certainly 
preferable to be less restrictive when applying  additional, 
non-merging-scale-related cuts to the matrix element. 
 
<p/> An example how to force a cut on lowest multiplicity reclustered 
states for pure QCD matrix element configurations is given by 
<code>main83.cc</code> (to be used with e.g. <code>main82.cmnd</code>). 
 
<h4>Influencing the construction of all possible histories</h4> 
 
<p/> Even more powerful - and dangerous - is influencing the construction 
of histories directly. This should only be attempted by expert users. If you 
believe manipulations completely unavoidable, we advise you to take great care 
when redefining the following functions. 
 
<a name="anchor12"></a>
<p/><strong>virtual bool MergingHooks::canCutOnRecState() &nbsp;</strong> <br/>
In the base class this method returns false. If you redefine it 
to return true then the method <code>doCutOnRecState(...)</code> 
will be called for each reclustered state encountered in the generation of 
all possible histories of the matrix element state. 
   
 
<a name="anchor13"></a>
<p/><strong>virtual bool MergingHooks::doCutOnRecState(const Event& event) &nbsp;</strong> <br/>
This routine will be supplied internally with every possible reclustered 
event that can be reached by reclustering any number of partons in 
the matrix element input state. The new, reclustered, states can then be 
analysed. If the method returns false, the history to which the state belongs 
will be treated as if it were unordered, i.e. this path will only be chosen 
if no other histories are available. In this way, the number of histories 
not fulfilling the user criterion will be minimised. 
   
 
<p/> 
Clearly, these methods are highly intrusive. It could e.g. happen that no 
history is allowed, which would make merging impossible. One example where 
this method could be useful is if cuts on the core <i>2 &rarr; 2</i> 
processes have to be checked, and the method 
<code>MergingHooks::dampenIfFailCuts(const Event& event)</code> is not 
sufficiently effective. 
 
<h4>Defining the hard process matrix element</h4> 
 
<p/> The MergingHooks class also allows the expert user to define the matrix 
element of the hard process, by defining the method 
 
<a name="anchor14"></a>
<p/><strong>virtual double MergingHooks::hardProcessME(const Event& inEvent) &nbsp;</strong> <br/>
This routine will be supplied internally with the reconstructed 
lowest-multiplicity event. From this, it is possible to calculate the squared 
matrix element of the hard process, by using the information stored in the 
event record. The function should return a <code>double</code> value that 
corresponds to the matrix element at the phase space point given by the input 
event record. This number will then be multiplied to the product of splitting 
functions that define the probability of the current path of the parton 
shower history. In this way, the hard process configuration can be taken into 
account when choosing the parton shower history, which is, internally, used 
to generate the "merging weight". 
   
 
<p/> The inclusion of the hard process matrix element into the choice 
of histories becomes relevant when the hard process matrix element has very 
strong phase space dependencies. QCD dijet cross sections for example strongly 
depend on the transverse momentum of the jets. So far, the authors have not 
encountered any changes upon inclusion of the full hard process matrix 
element, even for the QCD dijet case. 
 
<br/><br/><hr/> 
<h3>Matrix element merging and HepMC output for RIVET</h3> 
 
Examples how to produce matrix element merged events to be analysed 
with RIVET are given by <code>main84.cc</code> and <code>main85.cc</code>. 
 
<p/> The main issue is that the output of separate RIVET runs can not 
in general be combined. To perform a matrix element merging, we 
however need to runs over  different LHE files. The solution to this 
problem (so far) is to only perform  one RIVET run for all matrix 
elements, i.e. print the events for all ME parton  multiplicities, 
with the correct weights, to a single HepMC file. Since the correct 
weight includes the cross section of the different samples after 
Sudakov vetoes --- which is not a priori known --- the cross sections 
have to be  estimated in a test run, before the actual production run 
is performed. Finally, the cross section of the last event in the 
HepMC file has to be taken as the  full merged cross section 
<i>sigma_merge = Sum_{i=0}^N Sum_{j=0}*^{nEvents} 
sigma_est(i)*wckkwl(j)</i>. 
 
<p/> This procedure is outlined in <code>main84.cc</code>. 
 
Input LHE files with only very inclusive cuts pose further difficulties. For 
such files (which were already addressed under the heading <strong>Les Houches 
events outside the matrix element region</strong>), the cross section after 
the merging scale cut is not known before the cut is performed. Using Pythia's 
<code>UserHooks</code> facilities, it is possible to produce a valid estimate 
of the cross section after cuts. This however entails a careful cut definition 
by the user, which might become cumbersome for some in-built merging scale 
definitions. A reasonable alternative is using the switch 
 
<br/><br/><strong>Merging:doXSectionEstimate</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, estimate cross section after merging scale cut. This switch has to be 
used in conjunction with a merging scale definition (e.g. 
<code>Merging:doPTLundMerging = on</code>). Then, this merging scale 
definition will be used as a cut on the input events. After the requested 
number of Monte Carlo events, the cross section after the cut can be extracted 
by inferring the <code>Info::sigmaGen()</code> method, and the number of 
accepted events by using <code>Info::nAccepted()</code> 
   
 
<p/> 
This switch also relies on knowledge on how many partons a LHE file should 
contain. This is important for real-emission kinematics in the case of 
NLO merging. The number of (additional) partons in a LHE file can be set with 
 
<br/><br/><table><tr><td><strong>Merging:nRequested  </td><td></td><td> <input type="text" name="16" value="-1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>; <code>minimum = -1</code>)</td></tr></table>
Exact number of additional jets requested for a particular LHE file. If 
a file should for example only contain <i>W<sup>+</sup> g g </i> events, 
this switch should be set to "2" for this LHE file. For NLO merging 
schemes (see <?php $filepath = $_GET["filepath"];
echo "<a href='NLOMerging.php?filepath=".$filepath."' target='page'>";?>NLO Merging</a>), this number has to 
be set. 
   
 
<p/> 
The usage of these switches to obtain the necessary cross section estimate is 
illustrated in <code>main85.cc</code>. The example <code>main85.cc</code> 
program is intended as a "front-end" for CKKW-L merging in Pythia8, so we will 
discuss the program briefly. <code>main85.cc</code> should be used together 
with an input file (like <code>main85.cmnd</code>). The executable should be 
invoked with three arguments: the input file, the "name" of the input LHE 
files, and the name of the output HepMC event file. To use the LHE files that 
are shipped with the Pythia distribution, a valid usage would be 
<p/> 
<code>./main85.exe ./main85.cmnd ./w_production ./myhepmc.hepmc</code> 
<p/> 
If you want to use other input LHE files, note 
that <code>main85.cc</code> assumes the naming 
convention <i>name_tree_#nAdditionalJets.lhe</i>. All settings can 
be included though the input file, so that <code>main85.cc</code> does not 
have to be changed explicitly. <code>main85.cc</code> first switches off 
showers, multiparton interactions and hadronisation, and estimates the cross 
sections (after enforcing the merging scale cut) of samples for different 
numbers of additional jets. Then, showers, MPI and hadronisation are switched 
on again, and the CKKW-L merging procedure is performed. Events will be read 
in a decreasing sequence of jet multiplicities, which also means that e.g. 
events with two additional partons in the LHE file will be printed to the 
HepMC file before events with one additional parton. 
 
<br/><br/><hr/> 
<h3>Electroweak Merging</h3> 
 
Merging strategies like CKKW-L usually assume that the description of a 
(relatively simple) underlying process should be improved by combining with 
states that contain additional well-separated partons - with "additional" 
measured with respect to the underlying process. As discussed 
in [<a href="Bibliography.php#refChr15a" target="page">Chr15a</a>], this philosophy is not always sensible, and may lead to 
an unconvincing physics model. The bias can be greatly reduced by considering 
that in perturbation theory, corrections to seemingly very different 
underlying processes mix, so that there is no justification to classify 
some states as corrections to only one underlying process. Interactions 
that only contain vertices of only one theory (e.g. QCD) will mix with 
processes that contain only vertices of another interactions (e.g. QED). 
Underlying processes with very different coupling structures should thus 
be considered. This is the main aim of the electroweak merging scheme. 
The process <code> p p &rarr; W jet jet</code> provides a good example, since 
it can be interpreted either as double-real-QCD-emission correction to 
<code> p p &rarr; W</code> or as real-electroweak-emission 
correction to <code> p p &rarr; jet jet</code>. The distinction is artificial, 
but the all-order resummation is very different in either case, leading to 
distinctly different predictions. Thus, a minimally biased method for 
assigning an underlying process has to be found. 
 
<p/> 
The method of [<a href="Bibliography.php#refChr15a" target="page">Chr15a</a>] chooses the underlying process 
probabilistically for each phase space point based on a product of splitting 
kernels and full hard process matrix elements, and includes the correct 
all-order factors after this choice. In our previous example 
(<code> p p &rarr; W jet jet</code>), this would mean that the 
hard scattering could be <code> p p &rarr; W</code> or 
<code> p p &rarr; jet jet</code> (or, depending on phase space considerations, 
also <code> p p &rarr; W jet</code>). The availability of 
<?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>electro-weak showers</a> in Pythia 8 is crucial 
for a consistent scheme that removes of double-counting, since 
double-counting is (in part) removed by applying no-emission probabilities. 
In order to e.g. ensure that a state in the <code> p p &rarr; jet jet</code> 
does not, via W-boson emission, evolve into a state overlapping 
with <code> p p &rarr; W</code> + two QCD emissions, the former 
has to be reweighted with an all-order no-electroweak-emission probability. 
As always, a small merging scale dependence is facilitated by accounting for 
for dynamical PDF evaluation and running couplings. An electroweak merging 
thus includes a reweighting with &alpha;<sub>em</sub> ratios that are 
automatically included in the "merging weight". 
 
<p/> 
As another consequence of probabilistically assigning the underlying process 
is the "merging of mergings": Since purely partonic final states can 
evolve into jets + electroweak bosons states, it is necessary to treat bosons 
and partons on equal footing, meaning that well-separated boson states should 
be corrected with fixed-order inputs, while soft/collinear bosons should be 
associated with parton showering. Shortly, bosons and partons are treated 
identically. For our previous example, this means that the "correct" 
set of fixed-order corrections that include up to three final state particles 
is <code> p p &rarr; n Ws + m jets</code>, where any combination of 
n and m subject to <code> n + m &lt;= 3</code> has to be included 
(<code> p p &rarr; W</code>, <code> p p &rarr; W jet</code>, 
<code> p p &rarr; W jet jet</code>, <code> p p &rarr; W W</code>, 
<code> p p &rarr; W W W</code>, 
<code> p p &rarr; W W jet</code>, <code> p p &rarr; jet jet</code>, 
<code> p p &rarr; jet jet jet</code>) For practical purposes, it is sometimes 
permitted to not combine a complete set of processes. Only the single 
W-state has been explicitly validated and in addition the weak PS does 
not include all possible splittings for multiple W emissions, therefore 
caution has to be taken if using this for multiple W states. 
 
<p/> 
The "merging of mergings" has important consequences for the (interested) 
user. Below, we give instructions on the usage of the electroweak merging. 
 
<h4>Fixed-order inputs for electroweak merging</h4> 
 
The electroweak merging leads to the idea of a "merging of mergings". This has 
to be enforced also at the fixed-order sample generation stage, with two main 
requirements. 
 
<p/> A fully consistent treatment requires the generation of samples 
containing all states with a number of emissions that is less than 
or equal to the "maximal possible number of emissions of any type" that 
should be corrected, cf. the <code>W jet jet</code> example above. 
 In practise, it is often permitted 
to disregard some (set of) samples since their impact on an analysis is 
negligible. If e.g. an analyis always requires missing transverse momentum 
and a single lepton, and we assume perfect lepton acceptances, then it would 
be permitted to disregard the multi-W samples in the example. If the 
collision energy is in addition low (~up to fews of TeVs), then the 
probability for a pure QCD state to emit W-bosons is often low enough so 
that the pure jet samples can be neglected. However, you should think very 
carefully before settling on any shortcuts. 
 
<p/> Any particle that could count as an emission has to be included in the 
calculation of the particle separations that define the merging scale. If e.g. 
W-bosons are considered emissions, then any state with W-bosons that are 
collinear with another emission should be removed from the fixed-order sample. 
Such configurations will instead be produced through parton showering. This 
requirement means that you might have to define your own cut including this 
condition in your favourite fixed-order matrix element generator. At the 
risk of losses in efficiency, you can also use samples with very loose cuts 
and have Pythia enforce the merging 
scale cut when reading your input events. The latter is only possible if you 
use the merging scale definition <code>Merging:doPTLundMerging = on</code>. 
 
<h4>Enabling the electroweak merging</h4> 
 
The electroweak merging is currently only tested for processes containing 
W-bosons and jets. For a consistent merging, it is necessary to enable 
W-boson emissions by using 
<code>TimeShower:weakShower      = on</code>, 
<code>TimeShower:weakShowerMode  = 1</code>, 
<code>SpaceShower:weakShower     = on</code>, 
<code>SpaceShower:weakShowerMode = 1</code>, and 
<code>WeakShower:externalSetup   = on</code>. 
 
<p/> To enable the electroweak merging, use the following switch. 
 
<br/><br/><strong>Merging:allowWeakClustering</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow clustering of weak bosons, as necessary if a merging of matrix elements 
with QCD and weak showering is attempted. Currently, only emissions of 
W-bosons are accounted for. This switch should only be used if weak 
showering is turned on (see 
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>Timelike showers</a>, 
<?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>Spacelike showers</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>Weak showers</a> for details on weak 
showering). 
   
 
<p/> 
As explained above, the electroweak (EW) merging treats all shower-producible 
particles equally. This also means that the merging scale definition must 
include a separation of W-bosons and partons to define if a state is in the 
(well-separated) fixed-order region or if it is in the (soft/collinear) parton 
shower region. Such a cut can be implemented inside the (external) matrix 
element generator. On the other hand, Pythia 8 allows the usage of samples 
with very loose cuts and can enforce the correct merging scale cut by 
rejecting input events that do not pass the cut. This is also possible for 
the EW merging, albeit only for the merging scale definition that is enabled 
by using <code>Merging:doPTLundMerging = on</code>. We recommend using this 
strategy for users that do not wish to implement the cut directly into the ME 
generation, and who are prepared to accept a loss of efficiency because 
of Pythia's a-posteriori rejection. 
 
 
<h4>Defining the "inclusive" merging process</h4> 
 
Since the concept of a single "hard process" is not suitable for the EW 
merging, the process should be defined in a rather loose manner. 
This loose definition is still done by setting the 
input <code>Merging:Process</code>. 
 
Processes for EW merging should use the containers <code>Jinc</code>, 
<code>Winc</code>, 
<code>Ainc</code> and <code>Zinc</code>, which tell the code which particles 
could be possible "additional" emissions. No other particle defnitions are 
allowed, and none of the settings discussed in the "Defining the hard process" 
section are relevant here. 
 
Examples of allowed process definitions are 
 
<p/> <code>Merging:Process = pp > Jinc,Winc</code> meaning that W-bosons and 
partons are treated on equal footing (i.e. this is the setting applicable to 
the example used earlier). The merging will then include pure QCD multijet 
events, W+jets events, multi-W+jets events and pure multi-W events; 
 
<p/> <code>Merging:Process = pp > Jinc,Zinc</code> meaning that Z-bosons and 
partons are treated on equal footing; 
 
<p/> <code>Merging:Process = pp > Jinc,Winc,Zinc</code> meaning that W-bosons, 
Z-bosons and partons are treated on equal footing. 
 
<h4>Setting the number of additional particles</h4> 
 
Since the EW merging probabilistically decides on the "underlying process", it 
is a priori not possible to set the maximal number of additional emissions on 
top of this underlying process. A <code> p p &rarr; W jet jet</code> 
state would e.g. contain two additional emissions if interpreted as correction 
to <code> p p &rarr; W</code>, and only one additional emission 
if interpreted as correction to <code> p p &rarr; jet jet</code>. Pythia 8 
consequently decides dynamically how to set the additional number of 
emissions. 
 
<p/> 
The maximal number of emissions, set by using the <code>Merging:nJetMax</code> 
setting, still has to be defined to allow a sensible treatment of the 
"highest-multiplicity" states. We thus redefine the meaning of 
<code>Merging:nJetMax</code> to "maximal possible number of emissions of any 
type". As an example, <code>Merging:nJetMax = 3</code> if you want to perform 
a "merging of mergins" containing states with up to three partons, or up to 
two partons and one W-boson, or up to one parton and two W-bosons, or up to 
three W-bosons. 
 
 
 
<br/><br/><hr/> 
<h3>Further variables</h3> 
 
For more advanced manipulations of the merging machinery, all 
parameter  changes that were investigated in [<a href="Bibliography.php#refLon11" target="page">Lon11</a>] are 
supplied. Please  check [<a href="Bibliography.php#refLon11" target="page">Lon11</a>] for a detailed discussion of 
the switches. 
 
<p/> These switches allow enthusiastic users to perform a systematic 
assessment of the merging prescription. Apart from this, we advise the 
non-expert user to keep the default values. 
 
<br/><br/><table><tr><td><strong>Merging:nQuarksMerge  </td><td></td><td> <input type="text" name="18" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 2</code>; <code>maximum = 5</code>)</td></tr></table>
This switch controls which quarks flavours (labelled by PDG id's) are 
considered additional partons. If e.g. set to 4, then u-, d-, c- and s-quarks 
will be merged, while b-quarks will not be considered in the merging 
(corresponding to a 4-flavour merging scheme). We advise caution when 
changing this number. In particular, please ensure that the allowed flavour 
for additional partons in the input LHE file does not exceed this value, since 
unnecessary double-counting might occur otherwise. 
   
 
<br/><br/><strong>Merging:includeMassive</strong>  <input type="radio" name="19" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="19" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, use the correct massive evolution variable and massive 
splitting kernels in the reconstruction and picking of parton shower 
histories of the matrix  element. If off, reconstruct evolution 
scales, kinematics and splitting kernels  as if all partons were 
massless. 
   
 
<br/><br/><strong>Merging:enforceStrongOrdering</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, preferably pick parton shower histories of the matrix element 
which  have strongly ordered consecutive splittings, i.e. paths in 
which consecutive reclustered evolution scales are separated by a 
user-defined factor. 
   
 
<br/><br/><table><tr><td><strong>Merging:scaleSeparationFactor </td><td></td><td> <input type="text" name="21" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 1.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
The factor by which scales should differ to be classified as strongly 
ordered. 
   
 
<br/><br/><strong>Merging:orderInRapidity</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, preferably pick parton shower histories of the matrix element 
with  consecutive splittings ordered in rapidity and <i>pT</i>. 
   
 
<br/><br/><strong>Merging:pickByFullP</strong>  <input type="radio" name="23" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="23" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, pick parton shower histories of the matrix element by the full 
shower  splitting kernels, including potential ME corrections and 
Jacobians from joined evolution measures. 
   
 
<br/><br/><strong>Merging:pickByPoPT2</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, pick parton shower histories of the matrix element by the 
shower  splitting kernels divided by the evolution <i>pT</i>. 
   
 
<br/><br/><strong>Merging:pickBySumPT</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, exclusively pick parton shower histories of the matrix element 
for which have the smallest sum of scalar evolution <i>pT</i> for 
consecutive splittings has been calculated. 
   
 
<br/><br/><strong>Merging:includeRedundant</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, then also include PDF ratios and <i>&alpha;<sub>s</sub></i> 
factors in the  splitting probabilities used for picking a parton shower 
history of the matrix  element, when picking histories by the full shower 
splitting probability. As argued in  [<a href="Bibliography.php#refLon11" target="page">Lon11</a>], this should not 
be done since a reweighting with PDF ratios and <i>&alpha;<sub>s</sub></i> 
factors will be performed. However, it can give useful insight in how 
sensitive the results  are to the prescription on how to choose PS 
histories. 
   
 
<br/><br/><table><tr><td><strong>Merging:nonJoinedNorm </td><td></td><td> <input type="text" name="27" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
Normalisation factor with which to multiply splitting probability for 
splittings without joined evolution equation. 
   
 
<br/><br/><table><tr><td><strong>Merging:fsrInRecNorm </td><td></td><td> <input type="text" name="28" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
Normalisation factor with which to multiply splitting probability for 
final state splittings with an initial state recoiler. 
   
 
<br/><br/><table><tr><td><strong>Merging:aCollFSR </td><td></td><td> <input type="text" name="29" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
Factor with which to multiply the scalar <i>pT</i> of a final state 
splitting, when choosing the history by the smallest sum of scalar 
<i>pT</i>. Default value taken from Herwig++ [<a href="Bibliography.php#refTul09" target="page">Tul09</a>]. 
   
 
<br/><br/><table><tr><td><strong>Merging:aCollISR </td><td></td><td> <input type="text" name="30" value="0.9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
Factor with which to multiply the scalar <i>pT</i> of an initial state 
splitting, when choosing the history by the smallest sum of scalar 
<i>pT</i>. Default value taken from Herwig++ [<a href="Bibliography.php#refTul09" target="page">Tul09</a>]. 
   
 
<br/><br/><table><tr><td><strong>Merging:unorderedScalePrescrip  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
When the parton shower history of the matrix element contains a 
sequence of splittings which are not ordered in evolution <ei>pT</ei> 
(called an unordered history), this sequence is interpreted as a combined 
emission. Then, a decision on which starting scale for trial emissions 
off reconstructed states in this sequence of unordered splittings has 
to be made. Two options are available: 
<br/>
<input type="radio" name="31" value="0" checked="checked"><strong>0 </strong>:  Use larger of the two reconstructed (unordered)  scales as  starting scale.  <br/>
<input type="radio" name="31" value="1"><strong>1 </strong>:  Use smaller of the two reconstructed (unordered)  scales as  starting scale.  <br/>
 
<br/><br/><table><tr><td><strong>Merging:unorderedASscalePrescrip  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Prescription which scale to use to evaluate <ei>&alpha;<sub>s</sub></ei> 
weight for  splittings in a sequence of splittings which are not ordered 
in evolution <ei>pT</ei>. 
<br/>
<input type="radio" name="32" value="0"><strong>0 </strong>:  Use the combined splitting scale as argument in  <ei>&alpha;<sub>s</sub></ei>, for both splittings.  <br/>
<input type="radio" name="32" value="1" checked="checked"><strong>1 </strong>:  Use the true reconstructed scale  as as argument in  <ei>&alpha;<sub>s</sub></ei>, for each splitting separately.  <br/>
 
<br/><br/><table><tr><td><strong>Merging:unorderedPDFscalePrescrip  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Prescription which scale to use to evaluate ratios of parton distributions 
for splittings in a sequence of splittings which are not ordered 
in evolution <ei>pT</ei>. 
<br/>
<input type="radio" name="33" value="0" checked="checked"><strong>0 </strong>:  Use the combined splitting scale as argument in PDF ratios,  for both splittings.  <br/>
<input type="radio" name="33" value="1"><strong>1 </strong>:  Use the true reconstructed scale as argument in PDF  ratios, for each splitting separately.  <br/>
 
<br/><br/><table><tr><td><strong>Merging:incompleteScalePrescrip  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
When no complete parton shower history (i.e. starting from a 
<ei>2 &rarr; 2</ei> process)  for a matrix element with additional jets 
can be found, such a configuration is said to have an incomplete history. 
Since in incomplete histories, not all  shower starting scales are 
determined by clusterings, a prescription for setting the starting scale 
of trial showers in incomplete histories is needed. Three options are 
provided. 
<br/>
<input type="radio" name="34" value="0" checked="checked"><strong>0 </strong>:  Use factorisation scale as shower starting scale  for  incomplete histories.  <br/>
<input type="radio" name="34" value="1"><strong>1 </strong>:  Use <ei>sHat</ei> as shower starting scale for  incomplete histories.  <br/>
<input type="radio" name="34" value="2"><strong>2 </strong>:  Use <ei>s</ei> as shower starting scale for  incomplete histories.  <br/>
 
<br/><br/><strong>Merging:allowColourShuffling</strong>  <input type="radio" name="35" value="on"><strong>On</strong>
<input type="radio" name="35" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, this will allow the algorithm to swap one colour index in the state, 
when trying to find all possible clusterings, if no clustering has been 
found, but more clusterings had been requested. In this way, some incomplete 
histories can be avoided. Generally, we advise the non-expert user to not 
touch this switch, because a slight change in the colour structure can change 
the radiation pattern. To however study the sensitivity of the predictions on 
these effects, allowing for colour reshuffling can be useful. 
   
 
<br/><br/><strong>Merging:usePythiaQRenHard</strong>  <input type="radio" name="36" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="36" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, this will allow the algorithm to use a dynamical renormalisation scale 
to evaluate the strong couplings of the core hard process in dijet and 
prompt photon events. 
This means that the value of <i>&alpha;<sub>s</sub></i> used as coupling 
of the hard process in the matrix element generation will be replaced with 
a running coupling evaluated at the geometric mean of the squared transverse 
masses of the two outgoing particles, as is the default prescription in 
Pythia. 
   
 
<br/><br/><strong>Merging:usePythiaQFacHard</strong>  <input type="radio" name="37" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="37" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, this will allow the algorithm to use a dynamical factorisation scale 
to evaluate parton distributions associated with the hadronic cross section 
of the core hard process in dijet and prompt photon events. 
In the calculation of PDF ratios as part of the CKKW-L weight of an event, 
parton distributions that should be evaluated at the scale of the core 
2 - >2 process will be evaluated using the dynamical factorisation scale 
Pythia would attribute to this process. This means that the hard process 
factorisation scale is set to the smaller of the squared transverse masses 
of the two outgoing particles. 
   
 
<br/><br/><strong>Merging:mayRemoveDecayProducts</strong>  <input type="radio" name="38" value="on"><strong>On</strong>
<input type="radio" name="38" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Remove products of resonances in the hard process, in case Pythia generates 
decay products before merging. This makes merging possible even for an 
indeterminate final state, if Pythia itself has produced the decay products. 
The merging methods will instead be invoked on the "non-decayed" event, 
thus removing the limitation to only one decay channel when performing 
the merging. 
This switch is necessary e.g. for slepton pair production in association with 
additional QCD jets, if the input LHE file contains the resonant sleptons, 
and Pythia decides on a decay according to the branching fractions read from 
SLHA input. 
   
 
<br/><br/><strong>Merging:allowSQCDClustering</strong>  <input type="radio" name="39" value="on"><strong>On</strong>
<input type="radio" name="39" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow clustering of gluon emission off squarks. 
   
 
<br/><br/><strong>Merging:useShowerPlugin</strong>  <input type="radio" name="40" value="on"><strong>On</strong>
<input type="radio" name="40" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Use the splitting probabilities, evolution variables and phase space mappings 
of an external shower plugin. This will become possible as soon as new 
showers containing the necessary ingredients are available in Pythia. 
   
 
<br/><br/><strong>Merging:applyVeto</strong>  <input type="radio" name="41" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="41" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If off, no event veto based on the merging scale is applied in CKKW-L merging. 
This means that the user has to implement the veto by hand in the Pythia main 
program. It can be useful to postpone event vetoes for the purpose of merging 
scale variations. 
   
 
<br/><br/><strong>Merging:includeWeightInXsection</strong>  <input type="radio" name="42" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="42" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, then the reweighting of events in the CKKW-L scheme is included in 
the event weight <code>Info::weight()</code>, the merging weight 
<code>Info:mergingWeight()</code> is unity, and the cross section printed 
by <code>Info::sigmaGen()</code> includes the effect of CKKW-L merging. 
   
 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "off")
{
$data = "Merging:doKTMerging = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1")
{
$data = "Merging:ktType = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1.0")
{
$data = "Merging:Dparameter = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "Merging:nJetMax = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.0")
{
$data = "Merging:TMS = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "void")
{
$data = "Merging:Process = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "Merging:doMGMerging = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "Merging:doPTLundMerging = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "Merging:doCutBasedMerging = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.0")
{
$data = "Merging:QijMS = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.0")
{
$data = "Merging:pTiMS = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.0")
{
$data = "Merging:dRijMS = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "on")
{
$data = "Merging:enforceCutOnLHE = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Merging:doUserMerging = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "Merging:doXSectionEstimate = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "-1")
{
$data = "Merging:nRequested = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Merging:allowWeakClustering = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "5")
{
$data = "Merging:nQuarksMerge = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "on")
{
$data = "Merging:includeMassive = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "Merging:enforceStrongOrdering = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "1.0")
{
$data = "Merging:scaleSeparationFactor = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "Merging:orderInRapidity = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "on")
{
$data = "Merging:pickByFullP = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "Merging:pickByPoPT2 = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "Merging:pickBySumPT = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "Merging:includeRedundant = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "1.0")
{
$data = "Merging:nonJoinedNorm = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "1.0")
{
$data = "Merging:fsrInRecNorm = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "1.0")
{
$data = "Merging:aCollFSR = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "0.9")
{
$data = "Merging:aCollISR = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "0")
{
$data = "Merging:unorderedScalePrescrip = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "1")
{
$data = "Merging:unorderedASscalePrescrip = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "0")
{
$data = "Merging:unorderedPDFscalePrescrip = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "0")
{
$data = "Merging:incompleteScalePrescrip = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "off")
{
$data = "Merging:allowColourShuffling = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "on")
{
$data = "Merging:usePythiaQRenHard = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "on")
{
$data = "Merging:usePythiaQFacHard = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "Merging:mayRemoveDecayProducts = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "off")
{
$data = "Merging:allowSQCDClustering = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "off")
{
$data = "Merging:useShowerPlugin = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "on")
{
$data = "Merging:applyVeto = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "on")
{
$data = "Merging:includeWeightInXsection = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
