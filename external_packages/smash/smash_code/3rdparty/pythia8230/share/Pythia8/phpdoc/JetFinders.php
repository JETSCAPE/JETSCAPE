<html>
<head>
<title>Jet Finders</title>
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

<form method='post' action='JetFinders.php'>
 
<h2>Jet Finders</h2> 
 
<code>Pythia</code> comes with three <?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>built-in 
jet finders</a>, <code>ClusterJet</code> for <i>e^+e^-</i> events 
and <code>SlowJet</code> and <code>CellJet</code>for hadron collider ones. 
Especially the latter is not so well matched to the standards of its field, 
however. (But it is closely related to the anti-<i>kT</i> algorithm, 
so is also not completely disconnected [<a href="Bibliography.php#refCac08" target="page">Cac08</a>].) 
 
<p/> 
<code>SlowJet</code> can do jet finding according to the current-day 
<i>kT</i>, Cambridge/Aachen and anti-<i>kT</i> algorithms. 
It can be run in two modes. The original one is a native implementation 
which, as the name indicates, is rather slow. However, with the release 
of the <code>fjcore</code> code from <code>FastJet</code> 
[<a href="Bibliography.php#refCac06" target="page">Cac06</a>, <a href="Bibliography.php#refCac12" target="page">Cac12</a>], the default mode has become to use the 
<code>fjcore</code> methods. This is transparent to the user. 
 
<h3>FastJet</h3> 
 
<code>SlowJet</code> does not exhaust all the posssibilities of the 
<code>fjcore</code> code, so users are welcome to extend on the existing 
functionality by a direct usage of the <code>fjcore</code> methods. 
 
<p/> 
Missing from <code>fjcore</code> is a number of aspects, such as 
jet areas functionality, background estimation, access to other algorithms 
via plugins, interface to CGAL and tools such as filters and taggers. 
Therefore, for more sophisticated jet studies the complete 
<code>FastJet</code> package needs to be linked. This is foreseen in the 
configure file in the <code>examples</code> subdirectory, and the 
<code>main71.cc</code> and <code>main72.cc</code> programs contain 
examples how it can be used with <code>Pythia</code> events. (Even if 
these examples do not go beyond the functionality that 
<code>SlowJet</code> can offer.) 
 
<p/> 
The latter program makes use of the 
<code>include/Pythia8Plugins/FastJet3.h</code> 
header file, contributed by Gavin Salam. This allows simple input 
of a <code>Pythia</code> particle into a <code>FastJet</code> one, 
either retaining only the four-momentum or the full particle information. 
Thereby more sophisticated selectors become possible at the 
<code>FastJet</code> level. This code could be duplicated, with trivial 
modifications, to augment the <code>fjcore</code> package functionality 
in an identical manner, should the need arise. 
 
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
