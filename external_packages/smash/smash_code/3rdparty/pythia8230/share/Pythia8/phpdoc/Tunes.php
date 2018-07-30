<html>
<head>
<title>Tunes</title>
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

<form method='post' action='Tunes.php'>
 
<h2>Tunes</h2> 
 
Since some physics aspects cannot be derived from first principles, 
this program contains many parameters that represent a true 
uncertainty in our understanding of nature. Particularly afflicted 
are the areas of hadronization and multiparton interactions, which both 
involve nonperturbative QCD physics. 
 
<p/> 
Technically, PYTHIA  parameters can be varied independently of each 
other, but the physical requirement of a sensible description of a set 
of data leads to correlations and anticorrelations between the 
parameters. Hence the need to produce tunes, not of one parameter at 
a time, but simultaneously for a group of them. A well-known (separate) 
such example is parton densities, where combined tunes to a wide range 
of data have been produced, that can then be obtained prepackaged. 
 
<p/> 
Given the many PYTHIA parameters to be tuned, it is convenient to 
divide the task into subtasks. Firstly, if we assume jet universality, 
hadronization and final-state parton showers should be tuned to 
<i>e^+e^-</i> annihilation data, notably from LEP1, since this 
offers the cleanest environment. Secondly, with such parameters fixed, 
hadron collider data should be studied to pin down multiparton interactions 
and other further aspects, such as initial-state radiation. Ideally this 
would be done separately for diffractive and non-diffractive events, 
although it is not possible to have a clean separation. (Thirdly 
would come anything else, such as physics with photon beams, which 
involve further parameters, but that is beyond the current scope.) 
 
<p/> 
The first step was taken, with a tune to LEP1 data by Hendrik Hoeth, 
using the Rivet + Professor framework. Starting with version 8.125 it 
defined the default values for hadronization parameters and timelike 
showers. 
 
<p/> 
The situation is more complicated for hadronic interactions in general 
and multiparton interactions in particular, where PYTHIA 8 is more 
different from PYTHIA 6, and therefore more work is needed. Specifically, 
it is <i>not</i> possible to "port" a PYTHIA 6 tune to PYTHIA 8. 
 
<p/> 
A first simple tune, appropriately called "Tune 1", became default 
starting with version 8.127. It was noted, in particular by Hendrik 
Hoeth, that this tune had a tension between parameters needed to 
describe minimum-bias and underlying-event activity. Therefore some 
further physics features were introduced in the code itself 
[<a href="Bibliography.php#refCor10a" target="page">Cor10a</a>], which were made default as of 8.140. This version 
also included two new tunes, 2C and 2M, based on the CTEQ 6L1 and the 
MRST LO** PDF sets, respectively. These have been made by hand, as a 
prequel to complete Professor-style tunings. 
 
<p/> 
The very first data to come out of the LHC showed a higher rapidity 
plateau than predicted for current PYTHIA 6 tunes, also for the lower 
energies. This may suggest some tension in the data. Two alternatives, 
3C and 3M, were produced by a few brute-force changes of 2C and 2M. 
These were introduced in 8.140, but discontinued in 8.145 in favour of 
the newer 4C tune, that is based on a more serious study of some early 
LHC data, see [<a href="Bibliography.php#refCor10a" target="page">Cor10a</a>]. Following the comparative studies in 
[<a href="Bibliography.php#refBuc11" target="page">Buc11</a>], which independently confirmed a reasonable agreement 
with LHC data, tune 4C was made the default from 8.150. A variant is 
tune 4Cx, where the Gaussian matter profile has an <i>x</i>-dependent 
width [<a href="Bibliography.php#refCor11" target="page">Cor11</a>]. 
 
<p/> 
Tune 4C was the basis for many subsequent LHC tunes. Several ATLAS tunes 
have been included here, obtained with different PDFs and with different 
emphasis on minimum-bias and underlying-event data [<a href="Bibliography.php#refATL12" target="page">ATL12</a>]. 
These typically require LHAPDF to be linked, but this can be avoided 
in cases where the same PDF set is implemented internally. Also two CMS 
underlying-event tunes are implemented [<a href="Bibliography.php#refCMS14" target="page">CMS14</a>]. The ATLAS AZ tune 
[<a href="Bibliography.php#refATL14" target="page">ATL14</a>] instead puts emphasis on the <i>pT</i> spectrum of 
the <i>Z^/gamma^*0</i> boson. 
 
<p/> 
The Monash 2013 tune [<a href="Bibliography.php#refSka14" target="page">Ska14</a>] is based on a larger set of LHC 
distributions. It starts out from a more careful comparison with and 
tuning to LEP data, and so involves several parameter changes. The PDF 
used is the NNPDF2.3 QCD+QED LO one with <i>alpha_s(M_Z) = 0.130</i>, 
which includes more recent data than used in the previous default, and 
opens up for processes with incoming photons to the hard process. 
It is the default starting from 8.200. 
 
<p/> 
Recent tunes by the LHC collaborations are based on the Monash 2013 one. 
This includes the CMS tune MonashStar, or formally CUETP8M1-NNPDF2.3LO 
(currently unpublished). More significantly, ATLAS has produced and 
published a whole family for underlying-event tunes based on a major 
effort, including simultaneous fits of ten parameters [<a href="Bibliography.php#refATL14a" target="page">ATL14a</a>]. 
This includes four central tunes, with four different PDF sets, and 
ten variations around the NNPDF2.3 QCD+QED LO central tune, grouped in 
five pairs of variations up and down. The publication [<a href="Bibliography.php#refATL14a" target="page">ATL14a</a>] 
should be consulted for further details, like with what additional 
settings various processes have been generated, which should be respected 
to reap full benefit of the tunes. 
 
<p/> 
Central diffraction is a recent addition to the "soft QCD" process palette, 
and is thus not yet included in tunes; indeed its cross section is actively 
zeroed. You can switch it back on <i>after</i> you have selected your tune, 
with <code>SigmaTotal:zeroAXB = off</code>. But note that, since the 
total cross section is assumed unchanged, the nondiffractive cross section 
is reduced and thus also the MPI machinery affected, even if effects 
should not be big (for a small central diffractive cross section). 
 
<p/> 
Note that comparisons with data also require that other aspects agree, 
such as that decay chains are stopped at an agreed-on level. For instance, 
in the ATLAS tunes all particles with a lifetime above 10 mm 
are considered stable, <code>ParticleDecays:limitTau0 = on</code>, 
<code>ParticleDecays:tau0Max = 10</code>. We have chosen not to 
include this as part of the tune settings itself, since the tune as 
such could still be used with any other choice of stable and 
unstable particles. 
 
<p/> 
Further comparisons have been posted on the 
<a href="http://mcplots.cern.ch/">MCPLOTS</a> pages. 
They have been produced with help of the 
<a href="http://projects.hepforge.org/rivet/">Rivet</a> package 
[<a href="Bibliography.php#refBuc10" target="page">Buc10</a>]. 
 
<p/> 
To simplify comparisons for the user, we propose to collect some of 
the tunes here, in a prepackaged form. Of course, in all cases it is 
a matter of setting values for parameters already defined elsewhere, 
so the tunes offer no new functionality, only a more convenient setup. 
 
<p/> 
You should be aware that the evolution of the program will not guarantee 
complete backwards compatibility between versions. Most obviously this 
concerns bug fixes. But also for some other major changes, like the 
introduction of the new diffractive machinery, the default behaviour 
of old tunes has been changed retroactively. (Which should be fine for 
diffraction, since previous tunes were not based on data strongly 
influenced by diffraction.) 
 
<p/> 
The setup of the tunes is special, in that the choice of a tune forces 
the change of several different flags, modes and parameters. Furthermore 
a design principle has been that it should be possible to start out 
from a tune and then change a few of its settings. This gives power 
and flexibility at the expense of requiring a more careful ordering 
of commands. We therefore here sketch the order in which operations 
are carried out. 
<ol> 
<li>The constructor of a <code>Pythia</code> instance will read in 
all settings, and initialize them with their default values. 
</li> 
<li>At the end of this operation, the <code>Tune:ee</code> and 
<code>Tune:pp</code> modes (see further below) are checked. If either 
of them are positive the methods <code>Settings::initTuneEE(...)</code> 
and <code>Settings::initTunePP(...)</code>, respectively, are called 
to overwrite the whole collection of settings in the relevant tune. 
Zero (or negative) means that nothing will be done. 
Since most <i>pp/ppbar</i> tunes have been made in the context 
of an <i>e^+e^-</i> one, <code>initTunePP(...)</code> usually 
calls <code>initTuneEE(...)</code> to provide this synchronization. 
</li> 
<li>After the <code>Pythia</code> constructor all the relevant values 
for the default tune(s) have thus been set up. 
</li> 
<li>You as a user can now start to overwrite the values at will, 
using <code>Pythia::readFile(...)</code> to read a configuration file, 
or a list of <code>Pythia::readString(...)</code> commands, 
or the lower-level <code>Settings</code> methods. All changes 
are made in the order in which the commands are encountered during 
the execution. A given variable can be changed multiple times, 
but it is the latest change that sets the current value. 
</li> 
<li>The two <code>Tune:ee</code> and <code>Tune:pp</code> modes can also 
be changed in exactly the same way as described for all other settings 
above. Unique for them, however, is that when one of them is encountered 
it also initiates a call to the <code>initTuneEE(...)</code> or 
<code>initTunePP(...)</code> method, respectively. In such cases all 
settings affected by the <i>e^+e^-</i> or <i>pp/ppbar</i> tune 
are first reset to the default values (the <code>-1</code> options) 
and thereafter the relevant tune is set up. 
Recall that <code>initTunePP(...)</code> in its turn is allowed to call 
<code>initTuneEE(...)</code>. 
</li> 
<li>It is possible to mix commands of type 4 and 5 in any order; it 
is always the last change that counts. That is, any changes you have 
made to variables of a tune <i>before</i> a <code>Tune:ee</code> or 
<code>Tune:pp</code> command are overwritten by it, while variables 
you set <i>after</i> will overwrite the tune values. Further, 
the <code>Tune:pp</code> command usually implies an <i>e^+e^-</i> 
tune as well. Therefore <code>Tune:ee</code> would rarely be used for 
LHC applications. As a rule, instead, you want to begin with the 
<code>Tune:pp</code> choice, and thereafter modify only a small part 
of its settings. 
</li> 
<li>Needless to say, the flexibility can lead to unwanted setups if 
you do not exercise some discipline. It is therefore recommended that 
you always check the listing obtained with 
<code>Pythia::settings.listChanged()</code> to confirm that the 
final set of changes is the intended one. 
</li> 
</ol> 
 
<br/><br/><table><tr><td><strong>Tune:ee  </td><td>  &nbsp;&nbsp;(<code>default = <strong>7</strong></code>; <code>minimum = -1</code>; <code>maximum = 7</code>)</td></tr></table>
Choice of tune to <ei>e^+e^-</ei> data, mainly for the hadronization 
and timelike-showering aspects of PYTHIA. You should study the 
<code>Settings::initTuneEE(...)</code> method to find exactly which 
are the settings for the respective tune. 
<br/>
<input type="radio" name="1" value="-1"><strong>-1 </strong>: reset all values that are affected by any of the  <ei>e^+e^-</ei> tunes to the default values. This option can be used  on its own, but is also automatically used as a first step for either  of the positive tune values below, to undo the effect of previous tune  settings.  <br/>
<input type="radio" name="1" value="0"><strong>0 </strong>: no values are overwritten during the initial setup,  step 2 above. Note that changing to <code>0</code> in the user code  has no effect; if you want to restore the individual settings you  should instead use <code>-1</code>.  <br/>
<input type="radio" name="1" value="1"><strong>1 </strong>: the original PYTHIA 8 parameter set, based on some  very old flavour studies (with JETSET around 1990) and a simple tune  <ei>of alpha_strong</ei> to three-jet shapes to the new  <ei>pT</ei>-ordered shower. These were the default values before  version 8.125.  <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: a tune by Marc Montull to the LEP 1 particle  composition, as published in the RPP (August 2007). No related (re)tune  to event shapes has been performed, however.  <br/>
<input type="radio" name="1" value="3"><strong>3 </strong>: a tune to a wide selection of LEP1 data by Hendrik  Hoeth within the Rivet + Professor framework, both to hadronization and  timelike-shower parameters (June 2009). These were the default values  starting from version 8.125.  <br/>
<input type="radio" name="1" value="4"><strong>4 </strong>: a tune to LEP data by Peter Skands, by hand, both  to hadronization and timelike-shower parameters (September 2013).  Note the use of the CMW convention for the shower <ei>alpha_s</ei>  scale.  <br/>
<input type="radio" name="1" value="5"><strong>5 </strong>: first tune to LEP data by Nadine Fischer  (September 2013), based on the default flavour-composition  parameters. Input is event shapes (ALEPH and DELPHI),  identified particle spectra (ALEPH), multiplicities (PDG),  and B hadron fragmentation functions (ALEPH).  <br/>
<input type="radio" name="1" value="6"><strong>6 </strong>: second tune to LEP data by Nadine Fischer  (September 2013). Similar to the first one, but event shapes  are weighted up significantly, and multiplicites not included.  <br/>
<input type="radio" name="1" value="7" checked="checked"><strong>7 </strong>: the Monash 2013 tune by Peter Skands at al.  <ref>Ska14</ref>, to both <ei>e^+e^-</ei> and <ei>pp/pbarp</ei> data.  <br/>
 
<br/><br/><table><tr><td><strong>Tune:preferLHAPDF  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Tunes made by experimental collaborations typically use the LHAPDF 
package to obtain their PDF values, and so PYTHIA must be built 
accordingly. See the <aloc href="PDFSelection">PDF 
documentation</aloc> for more information. For PDFs implemented 
natively in PYTHIA it is possible to use the respective tunes, without 
having to use LHAPDF, if you set <code>Tune:preferLHAPDF = 0</code> 
<i>before</i> the <code>Tune:pp</code> choice. 
<br/>
<input type="radio" name="2" value="0"><strong>0 </strong>:  Use the internal PYTHIA PDFs.  <br/>
<input type="radio" name="2" value="1" checked="checked"><strong>1 </strong>:  Use LHAPDF5 PDFs.  <br/>
<input type="radio" name="2" value="2"><strong>2 </strong>:  Use LHAPDF6 PDFs.  <br/>
 
<br/><br/><table><tr><td><strong>Tune:pp  </td><td>  &nbsp;&nbsp;(<code>default = <strong>14</strong></code>; <code>minimum = -1</code>; <code>maximum = 34</code>)</td></tr></table>
Choice of tune to <ei>pp/ppbar</ei> data, mainly for the 
initial-state-radiation, multiparton-interactions and  beam-remnants 
aspects of PYTHIA. You should study the 
<code>Settings::initTunePP(...)</code> method to find exactly which 
are the settings for the respective tune. Note that all early tunes, 
including those done by the LHC collaborations based on tune 4C, 
imply the settings of <code>Tune:ee = 3</code>, while the Monash 2013 
tune and the further tunes based on it imply <code>Tune:ee = 7</code>. 
This is set automatically, and has to be overridden afterwards if not 
the wanted behaviour. 
<br/>
<input type="radio" name="3" value="-1"><strong>-1 </strong>: reset all values that are affected by any of the  <ei>pp/ppbar</ei> tunes to the default values. This option can be used  on its own, but is also automatically used as a first step for either  of the positive tune values below, to undo the effect of previous tune  settings.  <br/>
<input type="radio" name="3" value="0"><strong>0 </strong>: no values are overwritten during the initial setup,  step 2 above. Note that changing to <code>0</code> in the user code  has no effect; if you want to restore the individual settings you  should instead use <code>-1</code>.  <br/>
<input type="radio" name="3" value="1"><strong>1 </strong>: default used up to version 8.126, based on  some early and primitive comparisons with data.  <br/>
<input type="radio" name="3" value="2"><strong>2 </strong>: "Tune 1", default in 8.127 - 8.139, based on some  data comparisons by Peter Skands. Largely but not wholly overlaps  with the default option 0.  <br/>
<input type="radio" name="3" value="3"><strong>3 </strong>: "Tune 2C", introduced with 8.140 <ref>Cor10a</ref>.  It uses the CTEQ 6L1 PDF, and is intended to give good agreement with  much of the published CDF data.  <br/>
<input type="radio" name="3" value="4"><strong>4 </strong>: "Tune 2M", introduced with 8.140 <ref>Cor10a</ref>.  It is uses the MRST LO** PDF, which has a momentum sum somewhat above  unity, which is compensated by a smaller <ei>alpha_s</ei> than in the  previous tune. Again it is intended to give good agreement with much of  the published CDF data.  <br/>
<input type="radio" name="3" value="5"><strong>5 </strong>: "Tune 4C", newer tune, introduced with 8.145  <ref>Cor10a</ref>. Starts out from tune 2C, but with a reduced cross  section for diffraction, plus modified multiparton interactions parameters  to give a higher and more rapidly increasing charged pseudorapidity  plateau, for better agreement with some early key LHC numbers.  See also the comparative study in <ref>Buc11</ref>.  The starting point for many later tunes.  <br/>
<input type="radio" name="3" value="6"><strong>6 </strong>: "Tune 4Cx", based on tune 4C, but using the x-dependent  matter profile, <code>MultipartonInteractions:bProfile = 4</code> and an  increased <code>MultipartonInteractions:pT0Ref</code> <ref>Cor11</ref>.  <br/>
<input type="radio" name="3" value="7"><strong>7 </strong>: "ATLAS MB Tune A2-CTEQ6L1", a minimum-bias tune based  on tune 4Cx, but without rapidity-ordered spacelike emissions  <ref>ATL12</ref>. Uses CTEQ 6L1, by default from LHAPDF.  <br/>
<input type="radio" name="3" value="8"><strong>8 </strong>: "ATLAS MB Tune A2-MSTW2008LO", as above,  but uses MSTW 2008 LO, by default from LHAPDF.  <br/>
<input type="radio" name="3" value="9"><strong>9 </strong>: "ATLAS UE Tune AU2-CTEQ6L1", an underlying-event tune  based on tune 4Cx, but without rapidity-ordered spacelike emissions  <ref>ATL12</ref>. Uses CTEQ 6L1, by default from LHAPDF.  <br/>
<input type="radio" name="3" value="10"><strong>10 </strong>: "ATLAS UE Tune AU2-MSTW2008LO", as above,  but uses MSTW 2008 LO, by default from LHAPDF.  <br/>
<input type="radio" name="3" value="11"><strong>11 </strong>: "ATLAS UE Tune AU2-CT10", as above,  but uses CT 10, which is not currently implemented in PYTHIA,  so you must link LHAPDF.  <br/>
<input type="radio" name="3" value="12"><strong>12 </strong>: "ATLAS UE Tune AU2-MRST2007LO*", as above,  but uses MRST 2007 LO*, by default from LHAPDF.  <br/>
<input type="radio" name="3" value="13"><strong>13 </strong>: "ATLAS UE Tune AU2-MRST2007LO**", as above,  but uses MRST 2007 LO**, by default from LHAPDF.  <br/>
<input type="radio" name="3" value="14" checked="checked"><strong>14 </strong>: the Monash 2013 tune by Peter Skands at al.  <ref>Ska14</ref>, to both <ei>e^+e^-</ei> and <ei>pp/pbarp</ei> data.  The starting point for many later tunes.  <br/>
<input type="radio" name="3" value="15"><strong>15 </strong>: "CMS UE Tune CUETP8S1-CTEQ6L1", an underlying-event  tune based on tune 4C <ref>CMS14</ref>. Uses CTEQ 6L1, by default  from LHAPDF.  <br/>
<input type="radio" name="3" value="16"><strong>16 </strong>: "CMS UE Tune CUETP8S1-HERAPDF1.5LO", an underlying-event  tune based on tune 4C <ref>CMS14</ref>. Uses HERAPDF1.5LO, which is not  currently implemented in PYTHIA, so you must link LHAPDF.  <br/>
<input type="radio" name="3" value="17"><strong>17 </strong>: "ATLAS Tune AZ", is tuned to the <ei>pT</ei> spectrum  of the <ei>Z^/gamma^*0</ei> boson in a set of rapidity bins <ref>ATL14</ref>.  <br/>
<input type="radio" name="3" value="18"><strong>18 </strong>: "CMS Tune MonashStar", alias CUETP8M1-NNPDF2.3LO,  an underlying-event tune based on the Monash 2013 tune.  <br/>
<input type="radio" name="3" value="19"><strong>19 </strong>: "ATLAS A14 central tune with CTEQL1", a full-scale  tune to most ATLAS jet and underlying-event observables <ref>ATL14a</ref>,  starting out from the Monash 2013 tune. The following tunes 20 - 32  belong to the same group.  <br/>
<input type="radio" name="3" value="20"><strong>20 </strong>: "ATLAS A14 central tune with MSTW2008LO",  see above tune 19.  <br/>
<input type="radio" name="3" value="21"><strong>21 </strong>: "ATLAS A14 central tune with NNPDF2.3LO",  see above tune 19. Defines the center of the 23 - 32 variations,  so would be a good choice if you only want to study one tune from  the A14 family.  <br/>
<input type="radio" name="3" value="22"><strong>22 </strong>: "ATLAS A14 central tune with HERAPDF1.5LO",  see above tune 19. Uses HERAPDF1.5LO, which is not currently implemented  in PYTHIA, so you must link LHAPDF.  <br/>
<input type="radio" name="3" value="23"><strong>23 </strong>: "ATLAS A14 variation 1+" of tune 21.  <br/>
<input type="radio" name="3" value="24"><strong>24 </strong>: "ATLAS A14 variation 1-" of tune 21.  <br/>
<input type="radio" name="3" value="25"><strong>25 </strong>: "ATLAS A14 variation 2+" of tune 21.  <br/>
<input type="radio" name="3" value="26"><strong>26 </strong>: "ATLAS A14 variation 2-" of tune 21.  <br/>
<input type="radio" name="3" value="27"><strong>27 </strong>: "ATLAS A14 variation 3a+" of tune 21.  <br/>
<input type="radio" name="3" value="28"><strong>28 </strong>: "ATLAS A14 variation 3a-" of tune 21.  <br/>
<input type="radio" name="3" value="29"><strong>29 </strong>: "ATLAS A14 variation 3b+" of tune 21.  <br/>
<input type="radio" name="3" value="30"><strong>30 </strong>: "ATLAS A14 variation 3b-" of tune 21.  <br/>
<input type="radio" name="3" value="31"><strong>31 </strong>: "ATLAS A14 variation 3c+" of tune 21.  <br/>
<input type="radio" name="3" value="32"><strong>32 </strong>: "ATLAS A14 variation 3c-" of tune 21.  <br/>
<input type="radio" name="3" value="33"><strong>33 </strong>: tune that includes close-packing of strings and  hadron rescattering, Gaussian model for pT and flavour selection  <ref>Fis16</ref>. Based on Monash tune, mainly tuned to pT spectra.  <br/>
<input type="radio" name="3" value="34"><strong>34 </strong>: tune that includes close-packing of strings and  hadron rescattering, thermodynamical model for pT and flavour selection  <ref>Fis16</ref>. Based on Monash tune, mainly tuned to pT spectra.  <br/>
 
 
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

if($_POST["1"] != "7")
{
$data = "Tune:ee = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1")
{
$data = "Tune:preferLHAPDF = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "14")
{
$data = "Tune:pp = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
