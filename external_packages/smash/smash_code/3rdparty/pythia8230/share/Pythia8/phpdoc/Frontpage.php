<html>
<head>
<title>Front</title>
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

<form method='post' action='Frontpage.php'>
 
<h1>PYTHIA 8</h1> 
 
<h2>Welcome to PYTHIA - The Lund Monte Carlo!</h2> 
 
<p/> 
PYTHIA 8 is the successor to PYTHIA 6, rewritten from scratch in C++. 
At the release of the first public version, 8.100, it was untried and 
also still missed a few relevant features. This has changed over the 
years since then, and the release of 8.200 marks the end of a transition 
period. PYTHIA 8.2 has a significantly extended physics scope, notably 
for LHC physics, relative to what PYTHIA 6.4 could offer. There are only 
a few areas left where 6.4 still has a role to fill. 
 
<h2>Documentation</h2> 
 
On these webpages you will find the up-to-date manual for PYTHIA 8.2. 
Use the left-hand index to navigate this documentation of program 
elements, especially of all possible program settings. All parameters 
are provided with sensible default values, however, so you need only 
change those of relevance to your particular study, such as choice of 
beams, processes and phase space cuts. The pages also contain a fairly 
extensive survey of all methods available to the user, e.g. to study 
the produced events. What is lacking on these webpages is an overview, 
on the one hand, and an in-depth physics description, on the other. 
 
<p/> 
The overview can be found in the attached PDF file 
<br/><a href="../pdfdoc/pythia8200.pdf" target="page"> 
<b>An Introduction to PYTHIA 8.2</b></a> 
<br/>T. Sj&ouml;strand et al, Comput. Phys.Commun. 191 (2015) 159 
[arXiv:1410.3012 [hep-ph]]. 
<br/>You are strongly recommended to read this summary when you 
start out to learn how to use PYTHIA 8.2. 
 
<p/> 
For the physics description we refer to the complete 
<br/><b>PYTHIA 6.4 Physics and Manual</b> 
<br/>T. Sj&ouml;strand, S. Mrenna and P. Skands, JHEP05 (2006) 026, 
<br/>which in detail describes the physics (largely) implemented also in 
PYTHIA 8, and also provides a more extensive bibliography than found 
here. When you use PYTHIA 8.2, you should therefore cite both. 
 
<p/> 
Furthermore, a separate 
<br/><a href="../pdfdoc/worksheet8200.pdf" target="page"> 
<b>PYTHIA 8.2 Worksheet</b></a>, 
<br/>also an attached PDF file, offers a practical introduction to 
using the generator. It has been developed for and used at a few 
summer schools, with minor variations, but is also suited for 
self-study. 
 
<h2>Authors</h2> 
 
<p/> 
<b>Christian Bierlich</b><br/> 
Department of Astronomy and Theoretical Physics, Lund University, 
S&ouml;lvegatan 14A, SE-223 62 Lund, Sweden<br/> 
e-mail: christian.bierlich@thep.lu.se 
 
<p/> 
<b>Nishita Desai</b><br/> 
Laboratoire Charles Coulomb (L2C) & Laboratoire Univers et Particules 
de Montpellier (LUPM), CNRS-Universit&eacute; de Montpellier, 
34090 Montpellier, France<br/> 
e-mail: nishita.desai@umontpellier.fr 
 
<p/> 
<b>Nadine Fischer</b><br/> 
School of Physics, Monash University, PO Box 27, 3800 Melbourne, 
Australia<br/> 
e-mail: nadine.fischer@monash.edu 
 
<p/> 
<b>Ilkka Helenius</b><br/> 
Institute for Theoretical Physics, Tuebingen University, 
Auf der Morgenstelle 14, D-72076 Tuebingen, Germany<br/> 
e-mail: ilkka.helenius@uni-tuebingen.de 
 
<p/> 
<b>Philip Ilten</b><br/> 
School of Physics and Astronomy, University of Birmingham, 
Birmingham, B152 2TT, UK<br/> 
e-mail: philten@cern.ch 
 
<p/> 
<b>Leif L&ouml;nnblad</b><br/> 
Department of Astronomy and Theoretical Physics, Lund University, 
S&ouml;lvegatan 14A, SE-223 62 Lund, Sweden<br/> 
e-mail: leif.lonnblad@thep.lu.se 
 
<p/> 
<b>Stephen Mrenna</b><br/> 
Computing Division, Simulations Group, 
Fermi National Accelerator Laboratory, 
MS 234, Batavia, IL 60510, USA<br/> 
e-mail: mrenna@fnal.gov 
 
<p/> 
<b>Stefan Prestel</b><br/> 
Theoretical Physics Department, Fermi National Accelerator Laboratory, 
MS 106, Batavia, IL 60510, USA<br/> 
e-mail: sprestel@fnal.gov 
 
<p/> 
<b>Christine O. Rasmussen</b><br/> 
Department of Astronomy and Theoretical Physics, Lund University, 
S&ouml;lvegatan 14A, SE-223 62 Lund, Sweden<br/> 
e-mail: christine.rasmussen@thep.lu.se 
 
<p/> 
<b>Torbj&ouml;rn Sj&ouml;strand</b><br/> 
Department of Astronomy and Theoretical Physics, Lund University, 
S&ouml;lvegatan 14A, SE-223 62 Lund, Sweden<br/> 
e-mail: torbjorn@thep.lu.se 
 
<p/> 
<b>Peter Skands</b><br/> 
School of Physics, Monash University, PO Box 27, 3800 Melbourne, 
Australia<br/> 
e-mail: peter.skands@monash.edu 
 
<h2>Former authors</h2> 
 
<p/><b>Stefan Ask</b>, e-mail: ask.stefan@gmail.com 
 
<p/><b>Jesper Roy Christiansen</b>, e-mail: Jesper.Roy.Christiansen@thep.lu.se 
 
<p/><b>Richard Corke</b>, e-mail: r.corke@errno.net 
 
<h2>Further contributions</h2> 
 
Makefiles, configure scripts and HepMC interface by <b>Mikhail Kirsanov</b>. 
<br/>Conversion of XML files to PHP ones by <b>Ben Lloyd</b>. 
<br/>Simple Makefile for Win32/NMAKE by <b>Bertrand Bellenot</b>. 
<br/>Extended Higgs sector partly implemented by <b>Marc Montull</b>. 
<br/>Parts of charm and bottom decay tables courtesy <b>DELPHI</b> and 
<b>LHCb</b> collaborations. 
<br/>Tunes and comparisons with data, based on Rivet and Professor, 
by <b>Hendrik Hoeth</b>. 
<br/>Text and code on the use of ROOT in conjunction with PYTHIA 
by <b>Rene Brun</b>, <b>Andreas Morsch</b> and <b>Axel Naumann</b>. 
<br/>Code and data for MRST/MSTW PDFs by <b>Robert Thorne</b> and 
<b>Graeme Watt</b>. 
<br/>Code and data for the CTEQ/CT PDFs by <b>Joey Huston</b> 
and colleagues. 
<br/>Help with implementing new proton PDFs by <b>Tomas Kasemets</b>. 
<br/>Code and data for Pomeron PDFs by <b>H1</b> collaboration and 
especially <b>Paul Newman</b>. 
<br/>Help with implementing new Pomeron fluxes and PDFs by 
<b>Sparsh Navin</b>. 
<br/>The new Hidden Valley code developed together with <b>Lisa Carloni</b>. 
<br/>Code for a Kaluza-Klein electroweak gauge boson provided by 
<b>Noam Hod</b> and <b>Mark Sutton</b>. 
<br/>Code for equivalent photon flux around an unresolved proton by 
<b>Oystein Alvestad</b>. 
<br/>The MBR diffractive model and central diffraction by 
<b>Robert Ciesielski</b>. 
<br/>2012 branching ratios for most light hadrons, and the tau lepton, 
by <b>Anil Pratap Singh</b>. 
<br/>The pythia8-config script has been contributed by 
<b>Andy Buckley</b>, along with many other helpful suggestions. 
<br/>Code and data for several of the NNPDF2.3 QCD+QED sets provided by 
<b>Juan Rojo</b> and <b>Stefano Carrazza</b>. 
<br/>The fjcore code from FastJet provided by <b>Matteo Cacciari</b>, 
<b>Gavin Salam</b> and <b>Gregory Soyez</b>. 
<br/>The initial-final dipole approach has been developed and 
implemented by <b>Baptiste Cabouat</b>. 
<br/>The MixMax random number generator has been contributed by 
<b>Konstantin Savvidy</b> and <b>George Savvidy</b>. 
 
<br/><b>Note</b>: in several cases modifications have been made to 
the original code, in order to integrate it with PYTHIA. In these cases 
the blame for any mistakes has to rest with the regular authors. 
 
<h2>Licence</h2> 
 
PYTHIA 8 is licensed under the 
<a href="COPYING" target="page"><b>GNU General Public Licence 
version 2</b></a>. 
<br/>Please respect the 
<a href="GUIDELINES" target="page"><b>MCnet Guidelines</b></a> 
for Event Generator Authors and Users. 
 
<p/> 
The program and the documentation is 
Copyright &copy; 2017 Torbj&ouml;rn Sj&ouml;strand 
 
 
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
