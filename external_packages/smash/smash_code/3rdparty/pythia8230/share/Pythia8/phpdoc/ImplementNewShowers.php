<html>
<head>
<title>Implement New Showers</title>
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

<form method='post' action='ImplementNewShowers.php'>
 
<h2>Implement New Showers</h2> 
 
In case you want to replace the PYTHIA initial- and final-state 
showers by your own, it is possible but not trivial. The point is 
that multiparton interactions (MPI), initial-state radiation (ISR) and 
final-state radiation (FSR) in general appear in one single 
interleaved sequence of decreasing <i>pT</i> values. Therefore 
shower replacements would have to be able to play the game by such 
rules, as we will outline further below. Of course, this still 
leaves the field open exactly how to define what to mean by 
<i>pT</i>, how to handle recoil effects, how the colour flow is 
affected, and so on, so there is certainly room for alternative 
showers. A first example of a shower implemented within the PYTHIA 
context is <a href="http://home.fnal.gov/~skands/vincia/">VINCIA</a>. 
 
<p/> 
For the moment we assume you want to keep the MPI part of the story 
unchanged, and make use of the existing beam-remnants (BR) machinery. 
If you want to replace both MPI, ISR, FSR and BR then you had better 
replace the whole <code>PartonLevel</code> module of the code. 
If, in addition, you want to produce your own hard processes, 
then you only need the 
<?php $filepath = $_GET["filepath"];
echo "<a href='HadronLevelStandalone.php?filepath=".$filepath."' target='page'>";?>hadron-level standalone</a> 
part of the machinery. 
 
<p/> 
In order to write replacement codes for ISR and/or FSR it is useful 
to be aware of which information has to be shared between the 
different components, and which input/output structure is required 
of the relevant methods. For details, nothing beats studying the 
existing code. However, here we provide an overview, that should 
serve as a useful introduction. 
 
<p/> 
It should be noted that we here primarily address the problem in 
its full generality, with interleaved MPI, ISR and FSR. There exists 
an option <code>TimeShower:interleave = off</code> where only 
MPI and ISR would be interleaved and FSR be considered after these 
two, but still before BR. Most of the aspects described here would 
apply also for that case. By contrast, resonance decays are only 
considered after all the four above components, and timelike 
showers in those decays would never be interleaved with anything 
else, so are much simpler to administrate. 
 
<p/> 
Therefore the <code><?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?> 
pythia.setShowerPtr( timesDecPtr, timesPtr, spacePtr)</a></code> 
method allows two separate pointers to be set to instances of 
derived <code>TimeShower</code> classes. The first is only required 
to handle decays, say of <i>Z^0</i> or <i>Upsilon</i>, with no 
dependence on beam remnants or ISR. The second, as well as 
<code>spacePtr</code>, has to handle the interleaved evolution of MPI, 
ISR and FSR. Therefore you are free to implement only the first, and 
let the PYTHIA default showers take care of the latter two. But, if 
you wanted to, you could also set <code>timesDecPtr = 0</code> and 
only provide a <code>timesPtr</code>, or only a <code>spacePtr</code>. 
If your timelike shower does both cases, the first two pointers 
can agree. The only tiny point to take into account then is that 
<code>init( beamAPtr, beamBPtr)</code> is called twice, a first time 
to <code>timesDecPtr</code> with beam pointers 0, and a second time 
to <code>timesPtr</code> with nonvanishing beam pointers. 
 
<h3>The event record and associated information</h3> 
 
Obviously the main place for sharing information is the event 
record, specifically the <code>Event event</code> member of 
<code>Pythia</code>, passed around as a reference. It is 
assumed you already studied how it works, so here we only draw 
attention to a few aspects of special relevance. 
 
<p/> 
One basic principle is that existing partons should not be 
overwritten. Instead new partons should be created, even when a 
parton only receives a slightly shifted momentum and for the rest 
stays the same. Such "carbon copies" by final-state branchings 
should be denoted by both daughter indices of the original parton 
pointing to the copy, and both mother indices of the copy to the 
original. If the copy instead is intended to represent an earlier 
step, e.g. in ISR backwards evolution, the role of mothers and 
daughters is interchanged. The 
<code>event.copy( iCopy, newStatus)</code> 
routine can take care of this tedious task; the sign of 
<code>newStatus</code> tells the program which case to assume. 
 
<p/> 
To make the event record legible it is essential that the 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>status codes</a> 
are selected appropriately to represent the reason why each new 
parton is added to the record. Also remember to change the 
status of a parton to be negative whenever an existing parton 
is replaced by a set of new daughter partons. 
 
<p/> 
Another important parton property is <code>scale()</code>, 
which does not appear in the normal event listing, but only 
if you use the extended 
<code>Event:listScaleAndVertex = on</code> option. 
This property is supposed to represent the production scale 
(in GeV) of a parton. In the current FSR and ISR algorithms 
it is used to restrict from above the allowed <i>pT</i> 
values for branchings of this particular parton. 
Beam remnants and other partons that should not radiate are 
assigned scale 0. 
 
<p/> 
Auxiliary to the event record proper is the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='AdvancedUsage.php?filepath=".$filepath."' target='page'>";?>PartonSystems</a></code> 
class, that keep track of which partons belong together in the 
same scattering subsystem. This information must be kept up-to-date 
during the shower evolution. 
 
<p/> 
For initial-state showers it is also necessary to keep track of 
the partonic content extracted from the beams. This information 
is stored in the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='AdvancedUsage.php?filepath=".$filepath."' target='page'>";?>BeamParticle</a></code> 
class. 
 
<h3>The TimeShower interface</h3> 
 
If you want to replace the <code>TimeShower</code> class this would 
involve replacing the virtual methods among the following ones. 
 
<a name="anchor1"></a>
<p/><strong>TimeShower::TimeShower() &nbsp;</strong> <br/>
The constructor does not need to do anything. 
   
 
<a name="anchor2"></a>
<p/><strong>virtual TimeShower::~TimeShower() &nbsp;</strong> <br/>
The destructor does not need to do anything. 
   
 
<a name="anchor3"></a>
<p/><strong>void TimeShower::initPtr(Info* infoPtr, Settings* settingsPtr, ParticleData* particleDataPtr, Rndm* rndmPtr, CoupSM* coupSMPtr, PartonSystems* partonSystemsPtr, UserHooks* userHooksPtr) &nbsp;</strong> <br/>
This method only imports pointers to standard facilities, 
and is not virtual. 
   
 
<a name="anchor4"></a>
<p/><strong>virtual void TimeShower::init( BeamParticle* beamAPtrIn = 0, BeamParticle* beamBPtrIn = 0) &nbsp;</strong> <br/>
You have to store your local copy of the pointers to these objects, 
since they have to be used during the generation, as explained above. 
The pointers could be zero; e.g. a local copy of <code>TimeShower</code> 
is created to handle showers in decays such as <i>Upsilon &rarr; q qbar</i> 
from inside the <code>ParticleDecays class</code>. This is also the 
place to do initialization of whatever parameters you plan to use, 
e.g. by reading in them from a user-accessible database like the 
<code>Settings</code> one. 
   
 
<a name="anchor5"></a>
<p/><strong>virtual bool TimeShower::limitPTmax( Event& event, double Q2Fac = 0.,  double Q2Ren = 0.) &nbsp;</strong> <br/>
The question is whether the FSR should be allowed to occur at larger 
scales than the hard process it surrounds. This is process-dependent, 
as illustrated below for the the analogous 
<code>SpaeShower::limitPTmax(...)</code> method, although the two 
kinds of radiation need not have to be modeled identically. 
The <code>TimeShower:pTmaxMatch</code> switch allows you to force the 
behaviour among three options, but you may replace by your own logic. 
<br/>The internal PYTHIA implementation also allows intermediate options, 
where emissions can go up to the kinematical limit but be dampened above 
the factorization or renormalization scale. Therefore the (square of the) 
latter two are provided as optional input parameters. 
   
 
<a name="anchor6"></a>
<p/><strong>double TimeShower::enhancePTmax() &nbsp;</strong> <br/>
Relative to the default <i>pT_max</i> evolution scale of the process, 
it may still be convenient to vary the matching slightly for the hardest 
interaction in an event, to probe the sensitivity to such details. The 
base-class implementation returns the value of the 
<code>TimeShower:pTmaxFudge</code> parameter. 
   
 
<a name="anchor7"></a>
<p/><strong>virtual int TimeShower::shower( int iBeg, int iEnd, Event& event, double pTmax, int nBranchMax = 0) &nbsp;</strong> <br/>
This is an all-in-one call for shower evolution, and as such cannot be 
used for the normal interleaved evolution, where only the routines below 
are used. It also cannot be used in resonance decays that form part of 
the hard process, since there the 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hooks</a> insert a potential 
veto step. Currently this routine is therefore only used in the 
hadron-level decays, e.g. <i>Upsilon &rarr; g g g</i>. 
<br/><code>iBeg</code> and <code>iEnd</code> is the position of the 
first and last parton of a separate system, typically produced by a 
resonance decay. Such a system only evolves in isolation, and in 
particular does not relate to the beams. 
<br/>The <code>pTmax</code> value sets the maximum scale for evolution, 
but normally you would restrict that further for each individual 
parton based on its respective scale value. 
<br/>The <code>nBranchMax</code> value, if positive, gives the maximum 
number of allowed branchings in the call, as useful for matching studies. 
<br/>The routine is expected to return the number of FSR branchings that 
were generated, but only for non-critical statistics purposes. 
<br/>Since the real action typically is delegated to the routines 
below, it may well be that the existing code need not be replaced. 
   
 
<a name="anchor8"></a>
<p/><strong>virtual int TimeShower::showerQED( int iBeg, int iEnd, Event& event, double pTmax) &nbsp;</strong> <br/>
This is a further simplified version of the <code>shower</code> 
method above. Currently it only handles the emission of photons 
in the decay of a hadron into a pair of leptons, either a charged 
lepton-antilepton or a lepton-neutrino pair. It is properly matched 
to the matrix element in the decay via a virtual photon or 
<i>W^+-</i>, respectively. It is called as part of such decays 
if <code>ParticleDecays:allowPhotonRadiation = on</code>, which is 
not the default value. 
   
 
<a name="anchor9"></a>
<p/><strong>double TimeShower::pTLastInShower() &nbsp;</strong> <br/>
Can be used to return the <i>pT</i> evolution scale of the last 
branching in the cascade generated with the above 
<code>shower(...)</code> method. Is to be set in the internal 
<code>pTLastInShower</code> variable, and should be 0 if there 
were no branchings. Can be useful for matching studies. 
   
 
<a name="anchor10"></a>
<p/><strong>virtual void TimeShower::prepareGlobal( Event& event) &nbsp;</strong> <br/>
This method resets some counters and extracts the locations of 
outgoing partons, in preparation of using the optional global recoil 
scheme. Unlike <code>prepare(...)</code> below it is only called 
once during the parton-level evolution, since it only relates to 
the hardest interaction. Is probably of no use to most people. 
   
 
<a name="anchor11"></a>
<p/><strong>virtual void TimeShower::prepare( int iSys, Event& event, bool limitPTmaxIn = true) &nbsp;</strong> <br/>
This method is called immediately after a new interaction (or the 
products of a resonance decay) has been added, and should then be used 
to prepare the subsystem of partons for subsequent evolution. In 
the current code this involves identifying all colour and charge 
dipole ends: the position of radiating and recoiling partons, maximum 
<i>pT</i> scales, possible higher-order matrix elements matchings 
to apply, and so on. 
<br/>The <code>iSys</code> parameter specifies which parton system 
is to be prepared. It is used to extract the set of partons to be 
treated, with rules as described in the above section on subsystems. 
Specifically, the first two partons represent the incoming state, 
or are 0 for resonance decays unrelated to the beams, while the 
rest are not required to be in any particular order. 
<br/>The <code>limitPTmaxIn</code> switch conveys the choice made 
on maximum scale for the dipole-ends evolution associated with the 
hard interaction of the event (while scales for subsequent MPIs by 
default are set to respect <i>pT</i> ordering). If <code>true</code> 
then this scale is set by the user choice options, 
see the <code>TimeShower::limitPTmax</code> method above, 
while if <code>false</code> emissions are allowed to go up to the 
kinematical limit. For the two-hard-interactions scenario, the two 
class variables <code>dopTlimit1</code> and <code>dopTlimit2</code> 
instead convey the choice made. 
   
 
<a name="anchor12"></a>
<p/><strong>virtual void TimeShower::rescatterUpdate( int iSys, Event& event) &nbsp;</strong> <br/>
This method is called immediately after rescattering in the description 
of multiparton interactions. Thus the information on one or several 
systems is out-of-date, while that of the others is unchanged. 
We do not provide the details here, since we presume few implementors 
of new showers will want to touch the technicalities involved 
in obtaining a description of rescattering. 
   
 
<a name="anchor13"></a>
<p/><strong>virtual void TimeShower::update( int iSys, Event& event, bool hasWeakRad = false) &nbsp;</strong> <br/>
This method is called immediately after a spacelike branching in the 
<code>iSys</code>'th subsystem. Thus the information for that system is 
out-of-date, while that of the others is unchanged. If you want, you are 
free to throw away all information for the affected subsystem and call 
<code>prepare( iSys, event)</code> to create new one. Alternatively 
you may choose only to update the information that has changed. The 
optional final argument is used to signal if a weak radiation has 
occured in the spacelike evolution, which might be used to switch off 
the timelike weak emissions. 
   
 
<a name="anchor14"></a>
<p/><strong>virtual double TimeShower::pTnext( Event& event, double pTbegAll, double pTendAll, bool isFirstTrial = false) &nbsp;</strong> <br/>
This is the main driver routine for the downwards evolution. A new 
<i>pT</i> is to be selected based on the current information set up 
by the routines above, and along with that a branching parton or dipole. 
The <code>pTbegAll</code> scale is the maximum scale allowed, from which 
the downwards evolution should be begun (usually respecting the maximum 
scale of each individual parton). If no emission is found above 
<code>pTendAll</code> (and above the respective shower cutoff scales) 
then <code>0.</code> should be returned and no emissions will be allowed. 
Both scales can vary from one event to the next: if a scale has 
already been selected for MPI or ISR it makes no sense to look for 
a scale smaller than that from FSR, since it would not be able to 
compete, so <code>pTendAll</code> is set correspondingly. As it happens, 
FSR is tried before ISR and MPI in the interleaved evolution, 
but this is an implementation detail that could well change. 
<br/>Typically the implementation of this routine would be to set 
up a loop over all possible radiating objects (dipoles, dipole ends, ...), 
for each pick its possible branching scale and then pick the one 
with largest scale as possible winner. At this stage no branching should 
actually be carried out, since MPI, ISR and FSR still have to be compared 
to assign the winner. 
<br/>The <code>isFirstTrial</code> is <code>true</code> only for the 
very first emission, and is currently used in the global recoil option. 
   
 
<a name="anchor15"></a>
<p/><strong>virtual bool TimeShower::branch( Event& event, bool isInterleaved = false) &nbsp;</strong> <br/>
This method will be called once FSR has won the competition with 
MPI and ISR to do the next branching. The candidate branching found 
in the previous step should here be carried out in full. The 
pre-branching partons should get a negative status code and new 
replacement ones added to the end of the event record. Also the 
subsystem information should be updated, and possibly also the 
beams. 
<br/>Should some problem be encountered in this procedure, e.g. if 
some not-previously-considered kinematics requirement fails, it is 
allowed to return <code>false</code> to indicate that no branching 
could be carried out. 
<br/>Normally the optional <code>isInterleaved</code> argument would 
not be of interest. It can be used to separate resonance decays, false, 
from the interleaved evolution together with MPI and ISR, true. 
More precisely, it separates calls to the <code>timesDecPtr</code> 
and the <code>timesPtr</code> instances. 
   
 
<a name="anchor16"></a>
<p/><strong>virtual bool TimeShower::rescatterPropogateRecoil( Event& event, Vec4& pNew) &nbsp;</strong> <br/>
This method is only called if rescattering is switched on in the 
description of multiparton interactions. It then propagates a recoil 
from a timelike branching to internal lines that connect systems. 
As for <code>rescatterUpdate</code> above, this is not likely to be 
of interest to most implementors of new showers. 
   
 
<a name="anchor17"></a>
<p/><strong>int TimeShower::system() &nbsp;</strong> <br/>
This method is not virtual. If a branching is constructed by the 
previous routine this tiny method should be able to return the number 
of the selected subsystem <code>iSysSel</code> where it occurred, 
so that the spacelike shower can be told which system to update, 
if necessary. Therefore <code>iSysSel</code> must be set in 
<code>branch</code> (or already in <code>pTnext</code>). 
   
 
<a name="anchor18"></a>
<p/><strong>bool TimeShower::getHasWeaklyRadiated() &nbsp;</strong> <br/>
This method is not virtual. It is used to tell whether a weak gauge 
boson has been emitted in the shower evolution. 
   
 
<a name="anchor19"></a>
<p/><strong>virtual void TimeShower::list() &nbsp;</strong> <br/>
This method is not at all required. In the current implementation it 
outputs a list of all the dipole ends, with information on the 
respective dipole. The routine is not called anywhere in the public 
code, but has been inserted at various places during the 
development/debug phase. 
   
 
<p/> 
One of the major complications when combining fixed-order calculations 
with parton showers comes from constructing a parton shower history. This 
crucial step requires very precise knowledge of the functional forms 
of the evolution variable (e.g. transverse momentum), the auxiliary 
variable (e.g. an energy sharing variable), as well as the parton 
shower splitting probabilities and the (inverse) parton shower momentum 
mapping. All of these are readily available within the parton shower, 
but can be difficult to reconstruct, particularly if non-general improvements 
are included in the parton shower. To minimise code proliferation and 
bug potential, the timelike showers contain easy-access functions that can 
be loaded directly from PYTHIA's merging machinery. Note that within a new 
shower, you do not need to use these functions or all of the inputs 
transferred to these functions. Any dummy definition is acceptable. 
 
<p/> 
The easy-access functions are listed in the following. They use some common 
terminology, wherein <code>iRad</code> is the position in the event 
record of the radiating parton, <code>iRec</code> is the position of the 
recoiling partner of the dipole, which ensures overall energy and momentum 
conservation in the splitting, and <code>iEmt</code> is the position of the 
radiated parton, all after the splitting. (The ones before are stored in 
<code>iRadBef</code> and <code>iRecBef</code>, not used here.) Obviously 
the distinctions between <code>iRad</code> and <code>iEmt</code>, or 
between <code>iRad</code> and <code>iRec</code>, are a matter of choice. 
 
<a name="anchor20"></a>
<p/><strong>virtual bool TimeShower::isTimelike( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return true if the splitting that produced the particles 
<code>event[iRad]</code>, <code>event[iRec]</code> and 
<code>event[iEmt]</code> should be classified as timelike splittings 
(i.e. is handled by final state showers). The identifier <code>name</code> 
can be used for additional flexibility, e.g. if multiple kernels with 
identical post-branching states exist. 
   
 
<a name="anchor21"></a>
<p/><strong>virtual Event TimeShower::clustered( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return a PYTHIA event record in which the emission 
of the particle with index <code>iEmt</code> in the input 
<code>event</code> (also changing the particles with index <code>iRad</code> 
and <code>iRec</code>) is undone. The identifier 
<code>name</code> can be used for additional flexibility, e.g. if multiple 
kernels with identical post-branching states exist. 
Reclustered events are crucial in setting up consistent parton shower 
histories. 
   
 
<a name="anchor22"></a>
<p/><strong>virtual map &lt;string,double&gt; TimeShower::getStateVariables( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return a map of variables related to the splitting 
that produced the particles <code>event[iRad]</code>, <code>event[iRec]</code> 
and <code>event[iEmt]</code>. Four entries are compulsory, 
<br/>- the key <code>t</code> gives the evolution variable associated 
with the splitting, 
<br/>- the key <code>tRS</code> gives the evolution variable at 
which the shower should be restarted after a branching, 
<br/>- the key <code>scaleAS</code> gives the argument of 
<i>&alpha;<sub>s</sub></i> used for the branching, and 
<br/>- the key <code>scalePDF</code> gives the argument of PDF factors 
used for the branching. 
<br/>All other entries are optional and may be auxiliary variables related 
to the splitting (e.g. an energy sharing variable, an azimuthal angle, 
kinematical invariants etc.). 
<br/>The identifier <code>name</code> can be used for additional flexibility, 
e.g. if multiple kernels with identical post-branching states exist. 
   
 
<a name="anchor23"></a>
<p/><strong>virtual vector&lt;string&gt; TimeShower::getSplittingName( const Event& event, int iRad, int iEmt, int iRec) &nbsp;</strong> <br/>
This function should return a vector of string identifiers of the 
splitting producing the particles with indices 
<code>iRad</code>, <code>iEmt</code> and <code>iRec</code>in the input 
<code>event</code>. 
If e.g. <code>iRad</code> is the index of a final state quark 
and <code>iEmt</code> the index of a final state antiquark, we could imagine 
returning the string "fsr:G2QQ" (which is the name of such a branching 
in PYTHIA's UserHooks facilities). The return value is a vector to allow 
for multiple string identifiers, e.g. if multiple splittings lead to identical 
final states. 
   
 
<a name="anchor24"></a>
<p/><strong>virtual double TimeShower::getSplittingProb( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return the probability of an emission 
of the particle with index <code>iEmt</code> from the particles 
with index <code>iRad</code> and <code>iRec</code> 
All indices are relative to the input <code>event</code>. The identifier 
<code>name</code> can be used for additional flexibility. 
   
 
<a name="anchor25"></a>
<p/><strong>virtual bool TimeShower::allowedSplitting( const Event& event, int iRad, int iEmt) &nbsp;</strong> <br/>
This function is not used in the Pythia core code, and can thus be omitted. 
The purpose of this function is to allow ME+PS plugins for Pythia to ask the 
timelike shower if a particular combination of radiator and emission could 
have resulted from a splitting. The function arguments are the 
<code>event</code> which contains the radiator and emission, and the position 
of the radiator (<code>iRad</code>) and emission (<code>iEmt</code>) in said 
event. The function should return true if the combination was possible, and 
false otherwise. 
   
 
<a name="anchor26"></a>
<p/><strong>virtual vector&lt;int&gt;TimeShower::getRecoilers(  const Event& event, int iRad, int iEmt, string name) &nbsp;</strong> <br/>
This function is not used in the Pythia core code, and can thus be omitted. 
The purpose of this function is to allow ME+PS plugins for Pythia to ask the 
timelike shower for all possible recoilers for a particular combination of 
radiator and emission. The function arguments are the 
<code>event</code> which contains the radiator and emission, and the position 
of the radiator (<code>iRad</code>) and emission (<code>iEmt</code>) in said 
event. Since this information might not be enough to uniquely determine the 
list of recoilers, a further string argument that encodes e.g. the splitting 
name is supplied. The function should return a vector of the positions of all 
possible allowed recoilers in the input event. 
   
 
<h3>The SpaceShower interface</h3> 
 
If you want to replace the <code>SpaceShower</code> class this would 
involve replacing the virtual methods in the following. You will find 
that much of the story reminds of <code>TimeShower</code> above, and 
actually some cut-and-paste of text is involved. In some respects the 
description is simpler, since there are no special cases for resonance 
decays and non-interleaved evolution. Thus there is no correspondence 
to the <code>TimeShower::shower(...)</code> routine. 
 
<a name="anchor27"></a>
<p/><strong>SpaceShower::SpaceShower() &nbsp;</strong> <br/>
The constructor does not need to do anything. 
   
 
<a name="anchor28"></a>
<p/><strong>virtual SpaceShower::~SpaceShower() &nbsp;</strong> <br/>
Also the destructor does not need to do anything. 
   
 
<a name="anchor29"></a>
<p/><strong>void SpaceShower::initPtr(Info* infoPtrIn, Settings* settingsPtrIn, ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, CoupSM* coupSMPtrIn, PartonSystems* partonSystemsPtrIn, UserHooks* userHooksPtrIn, MergingHooks* mergingHooksPtrIn = 0) &nbsp;</strong> <br/>
This method only imports pointers to standard facilities, 
and is not virtual. 
   
 
<a name="anchor30"></a>
<p/><strong>virtual void SpaceShower::init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn) &nbsp;</strong> <br/>
You have to store your local copy of the pointers to these objects, 
since they have to be used during the generation, as explained above. 
This is also the place to do initialization of whatever parameters 
you plan to use, e.g. by reading in them from a user-accessible 
database like the <code>Settings</code> one. 
   
 
<a name="anchor31"></a>
<p/><strong>virtual bool SpaceShower::limitPTmax( Event& event, double Q2Fac = 0.,  double Q2Ren = 0.) &nbsp;</strong> <br/>
The question is whether the ISR should be allowed to occur at larger 
scales than the hard process it surrounds. This is process-dependent. 
For instance, if the hard process is <i>Z^0</i> production we know 
that ISR should be allowed to go right up to the kinematical limit. 
If it is a <i>2 &rarr; 2</i> QCD process the ISR should not exceed 
the scale of the hard process, since if so one would double-count. 
The <code>SpaceShower:pTmaxMatch</code> switch allows you to force the 
behaviour, or else to program your own logic. The current default 
implementation limits <i>pT</i> whenever the final state contains 
a quark (except top), gluon or photon, since then the danger of 
double-counting is there. You may replace by your own logic, or leave as is. 
<br/>The internal PYTHIA implementation also allows intermediate options, 
where emissions can go up to the kinematical limit but be dampened above 
the factorization or renormalization scale. Therefore the (square of the) 
latter two are provided as optional input parameters. 
   
 
<a name="anchor32"></a>
<p/><strong>virtual double SpaceShower::enhancePTmax() &nbsp;</strong> <br/>
When the above method limits <i>pT_max</i> to the scale of the process, 
it may still be convenient to vary the matching slightly for the hardest 
interaction in an event, to probe the sensitivity to such details. The 
base-class implementation returns the value of the 
<code>SpaceShower:pTmaxFudge</code> parameter. 
   
 
<a name="anchor33"></a>
<p/><strong>virtual void SpaceShower::prepare( int iSys, Event& event, bool limitPTmaxIn = true) &nbsp;</strong> <br/>
This method is called immediately after a new interaction has been 
added, and should then be used to prepare the subsystem of partons 
for subsequent evolution. In the current code this involves identifying 
the colour and charge dipole ends: the position of radiating and recoiling 
partons, maximum <i>pT</i> scales, and possible higher-order matrix 
elements matchings to apply. Depending on what you have in mind you may 
choose to store slightly different quantities. You have to use the 
subsystem information described above to find the positions of the two 
incoming partons (and the outgoing ones) of the system, and from there 
the scales at which they were produced. 
<br/> The <code>limitPTmax</code> input agrees with the output of the 
previous method for the hardest process, and is always true for 
subsequent MPI, since there an unlimited <i>pT</i> for sure 
would lead to double-counting. 
   
 
<a name="anchor34"></a>
<p/><strong>virtual void SpaceShower::update( int iSys, Event& event, bool hasWeakRad = false) &nbsp;</strong> <br/>
This method is called immediately after a timelike branching in the 
<code>iSys</code>'th subsystem. Thus the information for that system may 
be out-of-date, and to be updated. For the standard PYTHIA showers 
this routine does not need to do anything, but that may be different 
in another implementation. The optional final argument is used to signal 
if a weak radiation has occured in the timelike evolution, which might be 
used to switch off the spacelike weak emissions. 
   
 
<a name="anchor35"></a>
<p/><strong>virtual double SpaceShower::pTnext( Event& event, double pTbegAll, double pTendAll, int nRadIn = -1) &nbsp;</strong> <br/>
This is the main driver routine for the downwards evolution. A new 
<i>pT</i> is to be selected based on the current information set up 
by the routines above, and along with that a branching parton or dipole. 
The <code>pTbegAll</code> scale is the maximum scale allowed, from which 
the downwards evolution should be begun (usually respecting the maximum 
scale of each individual parton). If no emission is found above 
<code>pTendAll</code> (and above the respective shower cutoff scales) 
then <code>0.</code> should be returned and no emissions will be allowed. 
Both scales can vary from one event to the next: if a scale has 
already been selected for MPI or ISR it makes no sense to look for 
a scale smaller than that from FSR, since it would not be able to 
compete, so <code>pTendAll</code> is set correspondingly. As it happens, 
FSR is tried before ISR and MPI in the interleaved evolution, 
but this is an implementation detail that could well change. 
<br/>Typically the implementation of this routine would be to set 
up a loop over all possible radiating objects (dipoles, dipole ends, ...), 
for each pick its possible branching scale and then pick the one 
with largest scale as possible winner. At this stage no branching should 
actually be carried out, since MPI, ISR and FSR still have to be compared 
to assign the winner. 
<br/>The input <code>nRadIn</code> provides the total number of 
ISR and FSR emissions already generated in the event, and so allows a 
special treatment for the very first emission, if desired. 
   
 
<a name="anchor36"></a>
<p/><strong>virtual bool SpaceShower::branch( Event& event) &nbsp;</strong> <br/>
This method will be called once ISR has won the competition with 
MPI and FSR to do the next branching. The candidate branching found 
in the previous step should here be carried out in full. The 
pre-branching partons should get a negative status code and new 
replacement ones added to the end of the event record. Also the 
subsystem information should be updated, and possibly also the 
beams. 
<br/>Should some problem be encountered in this procedure, e.g. if 
some not-previously-considered kinematics requirement fails, it is 
allowed to return <code>false</code> to indicate that no branching 
could be carried out. Also a complete restart of the parton-level 
description may be necessary, see <code>doRestart()</code> below. 
   
 
<a name="anchor37"></a>
<p/><strong>int SpaceShower::system() &nbsp;</strong> <br/>
This method is not virtual. If a branching is constructed by the 
previous routine this tiny method should be able to return the number 
of the selected subsystem <code>iSysSel</code> where it occurred, 
so that the spacelike shower can be told which system to update, 
if necessary. Therefore <code>iSysSel</code> must be set in 
<code>branch</code> (or already in <code>pTnext</code>). 
   
 
<a name="anchor38"></a>
<p/><strong>bool SpaceShower::doRestart() &nbsp;</strong> <br/>
This method is not virtual. If <code>branch(...)</code> above fails 
to construct a branching, and the conditions are such that the whole 
parton-level description should be restarted, then it should return 
true, else not. Currently only the rescattering description can give 
this kind of failures, and so the internal <code>rescatterFail</code> 
boolean must be set true when this should happen, and else false. 
   
 
<a name="anchor39"></a>
<p/><strong>bool SpaceShower::getHasWeaklyRadiated() &nbsp;</strong> <br/>
This method is not virtual. It is used to tell whether a weak gauge 
boson has been emitted in the shower evolution. 
   
 
<a name="anchor40"></a>
<p/><strong>virtual void SpaceShower::list() &nbsp;</strong> <br/>
This method is not at all required. In the current implementation it 
outputs a list of all the dipole ends, with information on the 
respective dipole. The routine is not called anywhere in the public 
code, but has been inserted at various places during the 
development/debug phase. 
   
 
<p/> 
As mentioned at the end of the section on timelike showers, it can 
be beneficial to merging and matching machinery of PYTHIA if some 
additional functionality is available. This in particular includes 
easy access to a) the evolution variable, b) the auxiliary (energy-sharing) 
variable, c) to parton shower splitting probabilities and d) to the 
(inverse) parton shower momentum mapping. Thus, as in the timelike case, 
it can be beneficial to define the functions 
 
<a name="anchor41"></a>
<p/><strong>virtual bool SpaceShower::isSpacelike( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return true if the splitting that produced the particles 
<code>event[iRad]</code>, <code>event[iRec]</code> and 
<code>event[iEmt]</code> should be classified as spacelike splittings 
(i.e. is handled by initial state showers). The identifier <code>name</code> 
can be used for additional flexibility, e.g. if multiple kernels with 
identical post-branching states exist. 
   
 
<a name="anchor42"></a>
<p/><strong>virtual Event SpaceShower::clustered( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return a PYTHIA event record in which the emission 
of the particle with index <code>iEmt</code> in the input 
<code>event</code> (also changing the particles with index <code>iRad</code> 
and <code>iRec</code>) is undone. <code>name</code> 
is a string identifier for the splitting. Such reclustered events are 
crucial in setting up consistent parton shower histories. 
   
 
<a name="anchor43"></a>
<p/><strong>virtual map &lt;string,double&gt; SpaceShower::getStateVariables( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return a map of variables related to the splitting 
that produced the particles <code>event[iRad]</code>, <code>event[iRec]</code> 
and <code>event[iEmt]</code>. Four entries are compulsory, 
<br/>- the key <code>t</code> gives the evolution variable associated 
with the splitting, 
<br/>- the key <code>tRS</code> gives the evolution variable at 
which the shower should be restarted after a branching, 
<br/>- the key <code>scaleAS</code> gives the argument of 
<i>&alpha;<sub>s</sub></i> used for the branching, and 
<br/>- the key <code>scalePDF</code> gives the argument of PDF factors 
used for the branching. 
<br/>All other entries are optional and may be auxiliary variables related 
to the splitting (e.g. an energy sharing variable, an azimuthal angle, 
kinematical invariants etc.). 
<br/>The identifier <code>name</code> can be used for additional flexibility, 
e.g. if multiple kernels with identical post-branching states exist. 
   
 
<a name="anchor44"></a>
<p/><strong>virtual vector&lt;string&gt; SpaceShower::getSplittingName( const Event& event, int iRad, int iEmt, int iRec) &nbsp;</strong> <br/>
This function should return a string identifier of the splitting producing 
the particles with indices <code>iRad</code>, <code>iEmt</code> and 
<code>iRec</code> in the input <code>event</code>. 
If e.g. <code>iRad</code> is the index of an intial state quark 
and <code>iEmt</code> the index of a final state gluon, we could imagine 
returning the string "isr:Q2QG" (which is the name of such a branching 
in PYTHIA's UserHooks facilities). The return value is a vector to allow 
for multiple string identifiers, e.g. if multiple splittings lead to identical 
final states. 
   
 
<a name="anchor45"></a>
<p/><strong>virtual double SpaceShower::getSplittingProb( const Event& event, int iRad, int iEmt, int iRec, string name) &nbsp;</strong> <br/>
This function should return the probability of an emission 
of the particle with index <code>iEmt</code> from the particle 
with index <code>iRad</code> and <code>iRec</code>. All indices are 
relative to the input <code>event</code>. The identifier 
<code>name</code> can be used for additional flexibility. 
   
 
<a name="anchor46"></a>
<p/><strong>virtual bool SpaceShower::allowedSplitting( const Event& event, int iRad, int iEmt) &nbsp;</strong> <br/>
This function is not used in the Pythia core code, and can thus be omitted. 
The purpose of this function is to allow ME+PS plugins for Pythia to ask the 
spacelike shower if a particular combination of radiator and emission could 
have resulted from a splitting. The function arguments are the 
<code>event</code> which contains the radiator and emission, and the position 
of the radiator (<code>iRad</code>) and emission (<code>iEmt</code>) in said 
event. The function should return true if the combination was possible, and 
false otherwise. 
   
 
<a name="anchor47"></a>
<p/><strong>virtual vector&lt;int&gt; SpaceShower::getRecoilers( const Event& event, int iRad, int iEmt, string name) &nbsp;</strong> <br/>
This function is not used in the Pythia core code, and can thus be omitted. 
The purpose of this function is to allow ME+PS plugins for Pythia to ask the 
spacelike shower for all possible recoilers for a particular combination of 
radiator and emission. The function arguments are the 
<code>event</code> which contains the radiator and emission, and the position 
of the radiator (<code>iRad</code>) and emission (<code>iEmt</code>) in said 
event. Since this information might not be enough to uniquely determine the 
list of recoilers, a further string argument that encodes e.g. the splitting 
name is supplied. The function should return a vector of the positions of all 
possible allowed recoilers in the input event. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
