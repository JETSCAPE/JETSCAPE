<html>
<head>
<title>Diffraction</title>
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

<form method='post' action='Diffraction.php'>
 
<h2>Diffraction</h2> 
 
<h3>Introduction</h3> 
 
Diffraction is not well understood, and several alternative approaches 
have been proposed. Here we follow a fairly conventional Pomeron-based 
one, in the Ingelman-Schlein spirit [<a href="Bibliography.php#refIng85" target="page">Ing85</a>], 
but integrated to make full use of the standard PYTHIA machinery 
for multiparton interactions, parton showers and hadronization 
[<a href="Bibliography.php#refNav10" target="page">Nav10</a>,<a href="Bibliography.php#refCor10a" target="page">Cor10a</a>]. This is the approach pioneered in the PomPyt 
program by Ingelman and collaborators [<a href="Bibliography.php#refIng97" target="page">Ing97</a>]. 
 
<p/> 
For ease of use (and of modelling), the Pomeron-specific parts of the 
generation of inclusive ("soft") diffreactive events are subdivided into 
three sets of parameters that are rather independent of each other: 
<br/>(i) the total, elastic and diffractive cross sections are 
parametrized as functions of the CM energy, or can be set by the user 
to the desired values, see the 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>Total Cross Sections</a> page; 
<br/>(ii) once it has been decided to have a diffractive process, 
a Pomeron flux parametrization is used to pick the mass of the 
diffractive system(s) and the <i>t</i> of the exchanged Pomeron, 
see below; 
<br/>(iii) a diffractive system of a given mass is classified either 
as low-mass unresolved, which gives a simple low-<i>pT</i> string 
topology, or as high-mass resolved, for which the full machinery of 
multiparton interactions and parton showers are applied, making use of 
<?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>Pomeron PDFs</a>. 
<br/>The parameters related to multiparton interactions, parton showers 
and hadronization are kept the same as for normal nondiffractive events, 
with only one exception. This may be questioned, especially for the 
multiparton interactions, but we do not believe that there are currently 
enough good diffractive data that would allow detailed separate tunes. 
 
<p/> 
The above subdivision may not represent the way "physics comes about". 
For instance, the total diffractive cross section can be viewed as a 
convolution of a Pomeron flux with a Pomeron-proton total cross section. 
Since neither of the two is known from first principles there will be 
a significant amount of ambiguity in the flux factor. The picture is 
further complicated by the fact that the possibility of simultaneous 
further multiparton interactions ("cut Pomerons") will screen the rate of 
diffractive systems. In the end, our set of parameters refers to the 
effective description that emerges out of these effects, rather than 
to the underlying "bare" parameters. 
 
<p/> 
In the event record the diffractive system in the case of an excited 
proton is denoted <code>p_diffr</code>, code 9902210, whereas 
a central diffractive system is denoted <code>rho_diffr</code>, 
code 9900110. Apart from representing the correct charge and baryon 
numbers, no deeper meaning should be attributed to the names. 
 
<p/> 
PYTHIA also includes a possibility to select hard diffraction. This 
machinery relies on the same sets of parameters as described above, 
for the Pomeron flux and PDFs. The main difference between the hard 
and the soft diffractive frameworks is that the user can select any 
PYTHIA hard process in the former case, e.g. diffractive <i>Z</i>'s 
or <i>W</i>'s, whereas only QCD processes are generated in the latter. 
These QCD processes are generated inclusively, which means that mostly 
they occur in the low-<i>pT</i> region, even if a tail stretches to 
higher <i>pT</i> scales, thus overlapping with hard diffraction. 
Both hard and soft diffractive processes also include the normal PYTHIA 
machinery, such as MPIs and showers, but for the former the MPI 
framework can additionally be used to determine whether a hard process 
survives as a diffractive event or not. The different diffractive types 
- low mass soft, high mass soft and hard diffraction - are described 
in more detail below. 
 
<h3>Pomeron flux</h3> 
 
As already mentioned above, the total diffractive cross section is fixed 
by a default energy-dependent parametrization or by the user, see the 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>Total Cross Sections</a> page. 
Therefore we do not attribute any significance to the absolute 
normalization of the Pomeron flux. The choice of Pomeron flux model 
still will decide on the mass spectrum of diffractive states and the 
<i>t</i> spectrum of the Pomeron exchange. 
 
<br/><br/><table><tr><td><strong>Diffraction:PomFlux  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 7</code>)</td></tr></table>
Parametrization of the Pomeron flux <ei>f_Pom/p( x_Pom, t)</ei>. 
<br/>
<input type="radio" name="1" value="1" checked="checked"><strong>1 </strong>: Schuler and Sj&ouml;strand <ref>Sch94</ref>: based on a  critical Pomeron, giving a mass spectrum roughly like <ei>dm^2/m^2</ei>;  a mass-dependent exponential <ei>t</ei> slope that reduces the rate  of low-mass states; partly compensated by a very-low-mass (resonance region)  enhancement. Is currently the only one that contains a separate  <ei>t</ei> spectrum for double diffraction (along with MBR) and  separate parameters for pion beams.<br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: Bruni and Ingelman <ref>Bru93</ref>: also a critical  Pomeron giving close to <ei>dm^2/m^2</ei>,  with a <ei>t</ei> distribution  the sum of two exponentials. The original model only covers single  diffraction, but is here expanded by analogy to double and central  diffraction.<br/>
<input type="radio" name="1" value="3"><strong>3 </strong>: a conventional Pomeron description, in the RapGap  manual <ref>Jun95</ref> attributed to Berger et al. and Streng  <ref>Ber87a</ref>, but there (and here) with values updated to a  supercritical Pomeron with <ei>epsilon &gt; 0</ei> (see below),  which gives a stronger peaking towards low-mass diffractive states,  and with a mass-dependent (the <ei>alpha'</ei> below) exponential  <ei>t</ei> slope. The original model only covers single diffraction,  but is here expanded by analogy to double and central diffraction.  <br/>
<input type="radio" name="1" value="4"><strong>4 </strong>: a conventional Pomeron description, attributed to  Donnachie and Landshoff <ref>Don84</ref>, again with supercritical Pomeron,  with the same two parameters as option 3 above, but this time with a  power-law <ei>t</ei> distribution. The original model only covers single  diffraction, but is here expanded by analogy to double and central  diffraction.<br/>
<input type="radio" name="1" value="5"><strong>5 </strong>:  the MBR (Minimum Bias Rockefeller) simulation of  (anti)proton-proton interactions <ref>Cie12</ref>. The event  generation follows a renormalized-Regge-theory model, successfully tested  using CDF data. The simulation includes single and double diffraction,  as well as the central diffractive (double-Pomeron exchange) process (106).  Only <ei>p p</ei>, <ei>pbar p</ei> and <ei>p pbar</ei> beam combinations  are allowed for this option. Several parameters of this model are listed  below. <br/>
<input type="radio" name="1" value="6"><strong>6 </strong>:  The H1 Fit A parametrisation of the Pomeron flux  <ref>H1P06,H1P06a</ref>. The flux factors are motivated by Regge theory,  assuming a Regge trajectory as in options 3 and 4. The flux has been  normalised to 1 at <ei>x_Pomeron = 0.003</ei> and slope parameter and  Pomeron intercept has been fitted to H1 data.<br/>
<input type="radio" name="1" value="7"><strong>7 </strong>:  The H1 Fit B parametrisation of the Pomeron flux  <ref>H1P06,H1P06a</ref>.<br/>
 
<p/> 
In options 3, 4, 6, and 7 above, the Pomeron Regge trajectory is 
parametrized as 
<br/><i> 
alpha(t) = 1 + epsilon + alpha' t 
</i><br/> 
The <i>epsilon</i> and <i>alpha'</i> parameters can be set 
separately in options 3 and 4, and additionally <i>alpha'</i> 
is set in option 1: 
 
<br/><br/><table><tr><td><strong>Diffraction:PomFluxEpsilon </td><td></td><td> <input type="text" name="2" value="0.085" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.085</strong></code>; <code>minimum = 0.02</code>; <code>maximum = 0.15</code>)</td></tr></table>
The Pomeron trajectory intercept <i>epsilon</i> above for the 3 and 4 
flux options. For technical reasons <i>epsilon &gt; 0</i> is necessary 
in the current implementation. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:PomFluxAlphaPrime </td><td></td><td> <input type="text" name="3" value="0.25" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.25</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 0.4</code>)</td></tr></table>
The Pomeron trajectory slope <i>alpha'</i> above for the 1, 3 and 4 
flux options. 
   
 
Values are fixed in options 6 and 7. 
 
<p/> 
When option 5 is selected, the following parameters of the MBR model 
[<a href="Bibliography.php#refCie12" target="page">Cie12</a>] are used: 
 
<br/><br/><table><tr><td><strong>Diffraction:MBRepsilon </td><td></td><td> <input type="text" name="4" value="0.104" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.104</strong></code>; <code>minimum = 0.02</code>; <code>maximum = 0.15</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRalpha </td><td></td><td> <input type="text" name="5" value="0.25" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.25</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 0.4</code>)</td></tr></table>
the parameters of the Pomeron trajectory. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:MBRbeta0 </td><td></td><td> <input type="text" name="6" value="6.566" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>6.566</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRsigma0 </td><td></td><td> <input type="text" name="7" value="2.82" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.82</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
the Pomeron-proton coupling, and the total Pomeron-proton cross section. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:MBRm2Min </td><td></td><td> <input type="text" name="8" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 3.0</code>)</td></tr></table>
the lowest value of the mass squared of the dissociated system. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminSDflux </td><td></td><td> <input type="text" name="9" value="2.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminDDflux </td><td></td><td> <input type="text" name="10" value="2.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminCDflux </td><td></td><td> <input type="text" name="11" value="2.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
the minimum width of the rapidity gap used in the calculation of 
<i>Ngap(s)</i> (flux renormalization). 
   
 
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminSD </td><td></td><td> <input type="text" name="12" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminDD </td><td></td><td> <input type="text" name="13" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminCD </td><td></td><td> <input type="text" name="14" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
the minimum width of the rapidity gap used in the calculation of cross 
sections, i.e. the parameter <i>dy_S</i>, which suppresses the cross 
section at low <i>dy</i> (non-diffractive region). 
   
 
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminSigSD </td><td></td><td> <input type="text" name="15" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminSigDD </td><td></td><td> <input type="text" name="16" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5.0</code>)</td></tr></table>
<br/><br/><table><tr><td><strong>Diffraction:MBRdyminSigCD </td><td></td><td> <input type="text" name="17" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5.0</code>)</td></tr></table>
the parameter <i>sigma_S</i>, used for the cross section suppression at 
low <i>dy</i> (non-diffractive region). 
   
 
<h3>Separation of soft diffraction into low and high masses</h3> 
 
Preferably one would want to have a perturbative picture of the 
dynamics of Pomeron-proton collisions, like multiparton interactions 
provide for proton-proton ones. However, while PYTHIA by default 
will only allow collisions with a CM energy above 10 GeV, the 
mass spectrum of diffractive systems will stretch to down to 
the order of 1.2 GeV. It would not be feasible to attempt a 
perturbative description there. Therefore we do offer a simpler 
low-mass description, with only longitudinally stretched strings, 
with a gradual switch-over to the perturbative picture for higher 
masses. The probability for the latter picture is parametrized as 
<br/><i> 
P_pert = P_max ( 1 - exp( (m_diffr - m_min) / m_width ) ) 
</i><br/> 
which vanishes for the diffractive system mass 
<i>m_diffr &lt; m_min</i>, and is <i>1 - 1/e = 0.632</i> for 
<i>m_diffr = m_min + m_width</i>, assuming <i>P_max = 1</i>. 
 
<br/><br/><table><tr><td><strong>Diffraction:mMinPert </td><td></td><td> <input type="text" name="18" value="10." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.</strong></code>; <code>minimum = 5.</code>)</td></tr></table>
The abovementioned threshold mass <i>m_min</i> for phasing in a 
perturbative treatment. If you put this parameter to be bigger than 
the CM energy then there will be no perturbative description at all, 
but only the older low-<i>pt</i> description. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:mWidthPert </td><td></td><td> <input type="text" name="19" value="10." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned threshold width <i>m_width.</i> 
   
 
<br/><br/><table><tr><td><strong>Diffraction:probMaxPert </td><td></td><td> <input type="text" name="20" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
The abovementioned maximum probability <i>P_max.</i>. Would 
normally be assumed to be unity, but a somewhat lower value could 
be used to represent a small nonperturbative component also at 
high diffractive masses. 
   
 
<h3>Low-mass soft diffraction</h3> 
 
When an incoming hadron beam is diffractively excited, it is modeled 
as if either a valence quark or a gluon is kicked out from the hadron. 
In the former case this produces a simple string to the leftover 
remnant, in the latter it gives a hairpin arrangement where a string 
is stretched from one quark in the remnant, via the gluon, back to the 
rest of the remnant. The latter ought to dominate at higher mass of 
the diffractive system. Therefore an approximate behaviour like 
<br/><i> 
P_q / P_g = N / m^p 
</i><br/> 
is assumed. 
 
<br/><br/><table><tr><td><strong>Diffraction:pickQuarkNorm </td><td></td><td> <input type="text" name="21" value="5.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.0</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned normalization <i>N</i> for the relative quark 
rate in diffractive systems. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:pickQuarkPower </td><td></td><td> <input type="text" name="22" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
The abovementioned mass-dependence power <i>p</i> for the relative 
quark rate in diffractive systems. 
   
 
<p/> 
When a gluon is kicked out from the hadron, the longitudinal momentum 
sharing between the the two remnant partons is determined by the 
same parameters as above. It is plausible that the primordial 
<i>kT</i> may be lower than in perturbative processes, however: 
 
<br/><br/><table><tr><td><strong>Diffraction:primKTwidth </td><td></td><td> <input type="text" name="23" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of Gaussian distributions in <i>p_x</i> and <i>p_y</i> 
separately that is assigned as a primordial <i>kT</i> to the two 
beam remnants when a gluon is kicked out of a diffractive system. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:largeMassSuppress </td><td></td><td> <input type="text" name="24" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The choice of longitudinal and transverse structure of a diffractive 
beam remnant for a kicked-out gluon implies a remnant mass 
<i>m_rem</i> distribution (i.e. quark plus diquark invariant mass 
for a baryon beam) that knows no bounds. A suppression like 
<i>(1 - m_rem^2 / m_diff^2)^p</i> is therefore introduced, where 
<i>p</i> is the <code>diffLargeMassSuppress</code> parameter. 
   
 
<h3>High-mass soft diffraction</h3> 
 
The perturbative description need to use parton densities of the 
Pomeron. The options are described in the page on 
<?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>PDF Selection</a>. The standard 
perturbative multiparton interactions framework then provides 
cross sections for parton-parton interactions. In order to 
turn these cross section into probabilities one also needs an 
ansatz for the Pomeron-proton total cross section. In the literature 
one often finds low numbers for this, of the order of 2 mb. 
These, if taken at face value, would give way too much activity 
per event. There are ways to tame this, e.g. by a larger <i>pT0</i> 
than in the normal pp framework. Actually, there are many reasons 
to use a completely different set of parameters for MPI in 
diffraction than in pp collisions, especially with respect to the 
impact-parameter picture, see below. A lower number in some frameworks 
could alternatively be regarded as a consequence of screening, with 
a larger "bare" number. 
 
<p/> 
For now, however, an attempt at the most general solution would 
carry too far, and instead we patch up the problem by using a 
larger Pomeron-proton total cross section, such that average 
activity makes more sense. This should be viewed as the main 
tunable parameter in the description of high-mass diffraction. 
It is to be fitted to diffractive event-shape data such as the average 
charged multiplicity. It would be very closely tied to the choice of 
Pomeron PDF; we remind that some of these add up to less than unit 
momentum sum in the Pomeron, a choice that also affect the value 
one ends up with. Furthermore, like with hadronic cross sections, 
it is quite plausible that the Pomeron-proton cross section increases 
with energy, so we have allowed for a power-like dependence on the 
diffractive mass. 
 
<br/><br/><table><tr><td><strong>Diffraction:sigmaRefPomP </td><td></td><td> <input type="text" name="25" value="10." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.</strong></code>; <code>minimum = 2.</code>; <code>maximum = 40.</code>)</td></tr></table>
The assumed Pomeron-proton effective cross section, as used for 
multiparton interactions in diffractive systems. If this cross section 
is made to depend on the mass of the diffractive system then the above 
value refers to the cross section at the reference scale, and 
<br/><i> 
sigma_PomP(m) = sigma_PomP(m_ref) * (m / m_ref)^p 
</i><br/> 
where <i>m</i> is the mass of the diffractive system, <i>m_ref</i> 
is the reference mass scale <code>Diffraction:mRefPomP</code> below and 
<i>p</i> is the mass-dependence power <code>Diffraction:mPowPomP</code>. 
Note that a larger cross section value gives less MPI activity per event. 
There is no point in making the cross section too big, however, since 
then <i>pT0</i> will be adjusted downwards to ensure that the 
integrated perturbative cross section stays above this assumed total 
cross section. (The requirement of at least one perturbative interaction 
per event.) 
   
 
<br/><br/><table><tr><td><strong>Diffraction:mRefPomP </td><td></td><td> <input type="text" name="26" value="100.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.0</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The <i>mRef</i> reference mass scale introduced above. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:mPowPomP </td><td></td><td> <input type="text" name="27" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.5</code>)</td></tr></table>
The <i>p</i> mass rescaling pace introduced above. 
   
 
<p/> 
Also note that, even for a fixed CM energy of events, the diffractive 
subsystem will range from the abovementioned threshold mass 
<i>m_min</i> to the full CM energy, with a variation of parameters 
such as <i>pT0</i> along this mass range. Therefore multiparton 
interactions are initialized for a few different diffractive masses, 
currently five, and all relevant parameters are interpolated between 
them to obtain the behaviour at a specific diffractive mass. 
Furthermore, <i>A B &rarr; X B</i> and <i>A B &rarr; A X</i> are 
initialized separately, to allow for different beams or PDF's on the 
two sides. These two aspects mean that initialization of MPI is 
appreciably slower when perturbative high-mass diffraction is allowed. 
 
<p/> 
Diffraction tends to be peripheral, i.e. occur at intermediate impact 
parameter for the two protons. That aspect is implicit in the selection 
of diffractive cross section. For the simulation of the Pomeron-proton 
subcollision it is the impact-parameter distribution of that particular 
subsystem that should rather be modeled. That is, it also involves 
the transverse coordinate space of a Pomeron wavefunction. The outcome 
of the convolution therefore could be a different shape than for 
nondiffractive events. For simplicity we allow the same kind of 
options as for nondiffractive events, except that the 
<code>bProfile = 4</code> option for now is not implemented. 
 
<br/><br/><table><tr><td><strong>Diffraction:bProfile  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Choice of impact parameter profile for the incoming hadron beams. 
<br/>
<input type="radio" name="28" value="0"><strong>0 </strong>: no impact parameter dependence at all.<br/>
<input type="radio" name="28" value="1" checked="checked"><strong>1 </strong>: a simple Gaussian matter distribution;  no free parameters.<br/>
<input type="radio" name="28" value="2"><strong>2 </strong>: a double Gaussian matter distribution,  with the two free parameters <ei>coreRadius</ei> and  <ei>coreFraction</ei>.<br/>
<input type="radio" name="28" value="3"><strong>3 </strong>: an overlap function, i.e. the convolution of  the matter distributions of the two incoming hadrons, of the form  <ei>exp(- b^expPow)</ei>, where <ei>expPow</ei> is a free  parameter.<br/>
 
<br/><br/><table><tr><td><strong>Diffraction:coreRadius </td><td></td><td> <input type="text" name="29" value="0.4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.4</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 1.</code>)</td></tr></table>
When assuming a double Gaussian matter profile, <i>bProfile = 2</i>, 
the inner core is assumed to have a radius that is a factor 
<i>coreRadius</i> smaller than the rest. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:coreFraction </td><td></td><td> <input type="text" name="30" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
When assuming a double Gaussian matter profile, <i>bProfile = 2</i>, 
the inner core is assumed to have a fraction <i>coreFraction</i> 
of the matter content of the hadron. 
   
 
<br/><br/><table><tr><td><strong>Diffraction:expPow </td><td></td><td> <input type="text" name="31" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.4</code>; <code>maximum = 10.</code>)</td></tr></table>
When <i>bProfile = 3</i> it gives the power of the assumed overlap 
shape <i>exp(- b^expPow)</i>. Default corresponds to a simple 
exponential drop, which is not too dissimilar from the overlap 
obtained with the standard double Gaussian parameters. For 
<i>expPow = 2</i> we reduce to the simple Gaussian, <i>bProfile = 1</i>, 
and for <i>expPow &rarr; infinity</i> to no impact parameter dependence 
at all, <i>bProfile = 0</i>. For small <i>expPow</i> the program 
becomes slow and unstable, so the min limit must be respected. 
   
 
<h3>Hard diffraction</h3> 
 
When PYTHIA is requested to generate a hard process, by default it is 
assumed that the full perturbative cross section is associated with 
nondiffractive topologies. With the options in this section, PYTHIA 
includes a possibility for creating a perturbative process diffractively, 
however. This framework is denoted hard diffraction to distiguish it 
from soft diffraction, but recall that the latter contains a tail of 
high-<i>pT</i> processes that could alternatively be obtained as 
hard diffraction. The idea behind hard diffraction is similar to soft 
diffraction, as they are both based on the Pomeron model. The proton is 
thus modelled as having a Pomeron component, described by the Pomeron 
fluxes above, and the partonic content of the Pomeron is described by 
the Pomeron PDFs, also described above. From these components we can 
evaluate the probability for the chosen hard process to have been coming 
from a diffractively exited system, based on the ratio of the Pomeron flux 
convoluted with Pomeron PDF to the inclusive proton PDF. 
 
<p/> 
If the hard process is likely to have been created inside a diffractively 
excited system, then we also evaluate the momentum fraction carried by the 
Pomeron, <i>x_Pomeron</i>, and the momentum transfer, <i>t</i>, in 
the process. This information can also be extracted in the main programs, 
see eg. example <code>main61.cc</code>. 
 
<p/> 
Further, we distiguish between two alternative scenarios for the 
classification of hard diffraction. The first is based solely on the 
Pomeron flux and PDF, as described above. In the second an additional 
requirement is imposed, wherein the MPI machinery is not allowed to 
generate any extra MPIs at all, since presumably these would destroy 
the rapidity gap and thereby the diffractive nature. We refer to the 
former as MPI-unchecked and the latter as MPI-checked hard diffraction. 
The MPI-checked option is likely to be the more realistic one, but the 
MPI-unchecked one offers a convenient baseline for the study of MPI 
effects, which still are poorly understood. 
 
<p/> 
For the selected hard processes, the user can choose whether to generate 
the inclusive sample of both diffractive and nondiffractive topologies 
or diffractive only, and in each case with the diffractive ones 
distinguished either MPI-unchecked or MPI-checked. 
 
<br/><br/><strong>Diffraction:doHard</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allows for the possibility to include hard diffraction tests in a run. 
   
 
<p/> 
There is also the possibility to select only a specific subset of events 
in hard diffraction. 
 
<br/><br/><table><tr><td><strong>Diffraction:sampleType  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 4</code>)</td></tr></table>
Type of process the user wants to generate. Depends strongly on how an event 
is classified as diffractive. 
<br/>
<input type="radio" name="33" value="1"><strong>1 </strong>: Generate an inclusive sample of both diffractive and  nondiffractive hard processes, MPI-unchecked.  <br/>
<input type="radio" name="33" value="2" checked="checked"><strong>2 </strong>: Generate an inclusive sample of both diffractive and  nondiffractive hard processes, MPI-checked.  <br/>
<input type="radio" name="33" value="3"><strong>3 </strong>: Generate an exclusive diffractive sample, MPI-unchecked.  <br/>
<input type="radio" name="33" value="4"><strong>4 </strong>: Generate an exclusive diffractive sample, MPI-checked.  <br/>
 
<p/> 
The Pomeron PDFs have not been scaled to unit momentum sum by the 
H1 Collaboration, but instead they let the PDF normalization float 
after the flux had been normalized to unity at <i>x_Pom=0.03</i>. 
This means that the H1 Pomeron has a momentum sum that is about a half. 
It could be brought to unit momentum sum by picking the parameter 
<code>PDF:PomRescale</code> to be around 2. In order not to change the 
convolution of the flux and the PDFs, the flux then needs to be rescaled 
by the inverse. This introduces a new rescaling parameter: 
 
<br/><br/><table><tr><td><strong>Diffraction:PomFluxRescale </td><td></td><td> <input type="text" name="34" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 2.0</code>)</td></tr></table>
Rescale the Pomeron flux by this uniform factor. It should be 
<code>1 / PDF:PomRescale</code> to preserve the convolution of Pomeron 
flux and PDFs, but for greater flxibility the two can be set separately. 
   
 
<p/> 
When using the MBR flux, the model requires a renormalization of 
the Pomeron flux. This suppresses the flux with approximately a factor 
of ten, thus making it incompatible with the MPI suppression of the 
hard diffraction framework. We have thus implemented an option to 
renormalize the flux. If you wish to use the renormalized flux of the 
MBR model, you must generate the MPI-unchecked samples, otherwise 
diffractive events will be suppressed twice. 
 
<br/><br/><strong>Diffraction:useMBRrenormalization</strong>  <input type="radio" name="35" value="on"><strong>On</strong>
<input type="radio" name="35" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Use the renormalized MBR flux. 
   
 
<p/> 
The transverse matter profile of the Pomeron, relative to that of the 
proton, is not known. Generally a Pomeron is supposed to be a smaller 
object in a localized part of the proton, but one should keep an open 
mind. Therefore below you find three extreme scenarios, which can be 
compared to gauge the impact of this uncertainty. 
 
<br/><br/><table><tr><td><strong>Diffraction:bSelHard  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 3</code>)</td></tr></table>
Selection of impact parameter <ei>b</ei> and the related enhancement 
factor for the Pomeron-proton subsystem when the MPI check is carried 
out. This affects the underlying-event activity in hard diffractive 
events. 
<br/>
<input type="radio" name="36" value="1" checked="checked"><strong>1 </strong>: Use the same <ei>b</ei> as already assigned for the  proton-proton collision. This implicitly assumes that a Pomeron is  as big as a proton and centered in the same place. Since small  <ei>b</ei> values already have been suppressed, few events should  have high enhancement factors.  <br/>
<input type="radio" name="36" value="2"><strong>2 </strong>: Use the square root of the <ei>b</ei> as already  assigned for the proton-proton collision, thereby making the  enhancement factor fluctuate less between events. If the Pomeron  is very tiny then what matters is where it strikes the other proton,  not the details of its shape. Thus the variation with <ei>b</ei> is  of one proton, not two, and so the square root of the normal variation,  loosely speaking. Tecnhically this is difficult to implement, but  the current simple recipe provides the main effect of reducing the  variation, bringing all <ei>b</ei> values closer to the average.  <br/>
<input type="radio" name="36" value="3"><strong>3 </strong>: Pick a completely new <ei>b</ei>. This allows a broad  spread from central to peripheral values, and thereby also a more  varying MPI activity inside the diffractive system than the other two  options. This offers an extreme picture, even if not the most likely  one.  <br/>
 
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

if($_POST["1"] != "1")
{
$data = "Diffraction:PomFlux = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.085")
{
$data = "Diffraction:PomFluxEpsilon = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.25")
{
$data = "Diffraction:PomFluxAlphaPrime = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.104")
{
$data = "Diffraction:MBRepsilon = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.25")
{
$data = "Diffraction:MBRalpha = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "6.566")
{
$data = "Diffraction:MBRbeta0 = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "2.82")
{
$data = "Diffraction:MBRsigma0 = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "1.5")
{
$data = "Diffraction:MBRm2Min = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "2.3")
{
$data = "Diffraction:MBRdyminSDflux = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "2.3")
{
$data = "Diffraction:MBRdyminDDflux = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "2.3")
{
$data = "Diffraction:MBRdyminCDflux = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "2.0")
{
$data = "Diffraction:MBRdyminSD = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "2.0")
{
$data = "Diffraction:MBRdyminDD = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "2.0")
{
$data = "Diffraction:MBRdyminCD = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.5")
{
$data = "Diffraction:MBRdyminSigSD = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.5")
{
$data = "Diffraction:MBRdyminSigDD = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.5")
{
$data = "Diffraction:MBRdyminSigCD = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "10.")
{
$data = "Diffraction:mMinPert = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "10.")
{
$data = "Diffraction:mWidthPert = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "1.")
{
$data = "Diffraction:probMaxPert = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "5.0")
{
$data = "Diffraction:pickQuarkNorm = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "1.0")
{
$data = "Diffraction:pickQuarkPower = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.5")
{
$data = "Diffraction:primKTwidth = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "4.")
{
$data = "Diffraction:largeMassSuppress = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "10.")
{
$data = "Diffraction:sigmaRefPomP = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "100.0")
{
$data = "Diffraction:mRefPomP = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "0.0")
{
$data = "Diffraction:mPowPomP = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "1")
{
$data = "Diffraction:bProfile = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "0.4")
{
$data = "Diffraction:coreRadius = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "0.5")
{
$data = "Diffraction:coreFraction = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "1.")
{
$data = "Diffraction:expPow = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "Diffraction:doHard = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "2")
{
$data = "Diffraction:sampleType = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "1.0")
{
$data = "Diffraction:PomFluxRescale = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "off")
{
$data = "Diffraction:useMBRrenormalization = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "1")
{
$data = "Diffraction:bSelHard = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
