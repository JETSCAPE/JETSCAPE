<html>
<head>
<title>Four-Vectors</title>
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

<form method='post' action='FourVectors.php'>
 
<h2>Four-Vectors</h2> 
 
The <code>Vec4</code> class gives a simple implementation of four-vectors. 
The member function names are based on the assumption that these 
represent four-momentum vectors. Thus one can get or set 
<i>p_x, p_y, p_z</i> and <i>e</i>, but not <i>x, y, z</i> 
or <i>t</i>. This is only a matter of naming, however; a 
<code>Vec4</code> can equally well be used to store a space-time 
four-vector. 
 
<p/> 
The <code>Particle</code> object contains a <code>Vec4 p</code> that 
stores the particle four-momentum, and another <code>Vec4 vProd</code> 
for the production vertex. For the latter the input/output method 
names are adapted to the space-time character rather than the normal 
energy-momentum one. Thus a user would not normally access the 
<code>Vec4</code> classes directly, but only via the methods of the 
<code>Particle</code> class, 
see <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>Particle Properties</a>. 
 
<p/> 
Nevertheless you are free to use the PYTHIA four-vectors, e.g. as 
part of some simple analysis code based directly on the PYTHIA output, 
say to define the four-vector sum of a set of particles. But note that 
this class was never set up to allow complete generality, only  to 
provide the operations that are of use inside PYTHIA. There is no 
separate class for three-vectors, since such can easily be represented 
by four-vectors where the fourth component is not used. 
 
<p/> 
Four-vectors have the expected functionality: they can be created, 
copied, added, multiplied, rotated, boosted, and manipulated in other 
ways. Operator overloading is implemented where reasonable. Properties 
can be read out, not only the components themselves but also for derived 
quantities such as absolute momentum and direction angles. 
 
<h3>Constructors and basic operators</h3> 
 
A few methods are available to create or copy a four-vector: 
 
<a name="anchor1"></a>
<p/><strong>Vec4::Vec4(double x = 0., double y = 0., double z = 0., double t = 0.) &nbsp;</strong> <br/>
creates a four-vector, by default with all components set to 0. 
   
 
<a name="anchor2"></a>
<p/><strong>Vec4::Vec4(const Vec4& v) &nbsp;</strong> <br/>
creates a four-vector copy of the input four-vector. 
   
 
<a name="anchor3"></a>
<p/><strong>Vec4& Vec4::operator=(const Vec4& v) &nbsp;</strong> <br/>
copies the input four-vector. 
   
 
<a name="anchor4"></a>
<p/><strong>Vec4& Vec4::operator=(double value) &nbsp;</strong> <br/>
gives a  four-vector with all components set to <i>value</i>. 
   
 
<h3>Member methods for input</h3> 
 
The values stored in a four-vector can be modified in a few different 
ways: 
 
<a name="anchor5"></a>
<p/><strong>void Vec4::reset() &nbsp;</strong> <br/>
sets all components to 0. 
   
 
<a name="anchor6"></a>
<p/><strong>void Vec4::p(double pxIn, double pyIn, double pzIn, double eIn) &nbsp;</strong> <br/>
sets all components to their input values. 
   
 
<a name="anchor7"></a>
<p/><strong>void Vec4::p(Vec4 pIn) &nbsp;</strong> <br/>
sets all components equal to those of the input four-vector. 
   
 
<a name="anchor8"></a>
<p/><strong>void Vec4::px(double pxIn) &nbsp;</strong> <br/>
   
<a name="anchor9"></a>
<strong>void Vec4::py(double pyIn) &nbsp;</strong> <br/>
   
<a name="anchor10"></a>
<strong>void Vec4::pz(double pzIn) &nbsp;</strong> <br/>
   
<a name="anchor11"></a>
<strong>void Vec4::e(double eIn) &nbsp;</strong> <br/>
sets the respective component to the input value. 
   
 
<h3>Member methods for output</h3> 
 
A number of methods provides output of basic or derived quantities: 
 
<a name="anchor12"></a>
<p/><strong>double Vec4::px() &nbsp;</strong> <br/>
   
<a name="anchor13"></a>
<strong>double Vec4::py() &nbsp;</strong> <br/>
   
<a name="anchor14"></a>
<strong>double Vec4::pz() &nbsp;</strong> <br/>
   
<a name="anchor15"></a>
<strong>double Vec4::e() &nbsp;</strong> <br/>
gets the respective component. 
   
 
<a name="anchor16"></a>
<p/><strong>double& operator[](int i) &nbsp;</strong> <br/>
returns component by index, where 1 gives <i>p_x</i>, 2 gives <i>p_y</i>, 
3 gives <i>p_z</i>, and anything else gives <i>e</i>. 
   
 
<a name="anchor17"></a>
<p/><strong>double Vec4::mCalc() &nbsp;</strong> <br/>
   
<a name="anchor18"></a>
<strong>double Vec4::m2Calc() &nbsp;</strong> <br/>
the (squared) mass, calculated from the four-vectors. 
If <i>m^2 &lt; 0</i> the mass is given with a minus sign, 
<i>-sqrt(-m^2)</i>.  Note the possible loss of precision 
in the calculation of <i>E^2 - p^2</i>; for particles the 
correct mass is stored separately to avoid such problems. 
   
 
<a name="anchor19"></a>
<p/><strong>double Vec4::pT() &nbsp;</strong> <br/>
   
<a name="anchor20"></a>
<strong>double Vec4::pT2() &nbsp;</strong> <br/>
the (squared) transverse momentum. 
   
 
<a name="anchor21"></a>
<p/><strong>double Vec4::pAbs() &nbsp;</strong> <br/>
   
<a name="anchor22"></a>
<strong>double Vec4::pAbs2() &nbsp;</strong> <br/>
the (squared) absolute momentum. 
   
 
<a name="anchor23"></a>
<p/><strong>double Vec4::eT() &nbsp;</strong> <br/>
   
<a name="anchor24"></a>
<strong>double Vec4::eT2() &nbsp;</strong> <br/>
the (squared) transverse energy, 
<i>eT = e * sin(theta) = e * pT / pAbs</i>. 
   
 
<a name="anchor25"></a>
<p/><strong>double Vec4::theta() &nbsp;</strong> <br/>
the polar angle, in the range 0 through 
<i>pi</i>. 
   
 
<a name="anchor26"></a>
<p/><strong>double Vec4::phi() &nbsp;</strong> <br/>
the azimuthal angle, in the range <i>-pi</i> through <i>pi</i>. 
   
 
<a name="anchor27"></a>
<p/><strong>double Vec4::thetaXZ() &nbsp;</strong> <br/>
the angle in the <i>xz</i> plane, in the range <i>-pi</i> through 
<i>pi</i>, with 0 along the <i>+z</i> axis. 
   
 
<a name="anchor28"></a>
<p/><strong>double Vec4::pPos() &nbsp;</strong> <br/>
   
<a name="anchor29"></a>
<strong>double Vec4::pNeg() &nbsp;</strong> <br/>
the combinations <i>E+-p_z</i>.   
 
<a name="anchor30"></a>
<p/><strong>double Vec4::rap() &nbsp;</strong> <br/>
   
<a name="anchor31"></a>
<strong>double Vec4::eta() &nbsp;</strong> <br/>
true rapidity <i>y</i> and pseudorapidity <i>eta</i>. 
   
 
<h3>Friend methods for output</h3> 
 
There are also some <code>friend</code> methods that take one, two 
or three four-vectors as argument. Several of them only use the 
three-vector part of the four-vector. 
 
<a name="anchor32"></a>
<p/><strong>friend ostream& operator&lt;&lt;(ostream&, const Vec4& v) &nbsp;</strong> <br/>
writes out the values of the four components of a <code>Vec4</code> and, 
within brackets, a fifth component being the invariant length of the 
four-vector, as provided by <code>mCalc()</code> above, and it all 
ended with a newline. 
   
 
<a name="anchor33"></a>
<p/><strong>friend double m(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
   
<a name="anchor34"></a>
<strong>friend double m2(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
the (squared) invariant mass. 
   
 
<a name="anchor35"></a>
<p/><strong>friend double dot3(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
the three-product. 
   
 
<a name="anchor36"></a>
<p/><strong>friend double cross3(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
the cross-product. 
   
 
<a name="anchor37"></a>
<p/><strong>friend double cross4(const Vec4& v1, const Vec4& v2, const Vec4& v3) &nbsp;</strong> <br/>
the cross-product of three four-vectors: 
<i>v_i = epsilon_{iabc} v1_a v2_b v3_c</i>. 
   
 
<a name="anchor38"></a>
<p/><strong>friend double theta(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
   
<a name="anchor39"></a>
<strong>friend double costheta(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
the (cosine) of the opening angle between the vectors, 
in the range 0 through <i>pi</i>. 
   
 
<a name="anchor40"></a>
<p/><strong>friend double phi(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
   
<a name="anchor41"></a>
<strong>friend double cosphi(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
the (cosine) of the azimuthal angle between the vectors around the 
<i>z</i> axis, in the range 0 through <i>pi</i>. 
   
 
<a name="anchor42"></a>
<p/><strong>friend double phi(const Vec4& v1, const Vec4& v2, const Vec4& v3) &nbsp;</strong> <br/>
   
<a name="anchor43"></a>
<strong>friend double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& v3) &nbsp;</strong> <br/>
the (cosine) of the azimuthal angle between the first two vectors 
around the direction of the third, in the range 0 through <i>pi</i>. 
   
 
<a name="anchor44"></a>
<p/><strong>friend double RRapPhi(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
   
<a name="anchor45"></a>
<strong>friend double REtaPhi(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
the <i>R</i> distance measure, in <i>(y, phi)</i> or 
<i>(eta, phi)</i> cylindrical coordinates, i.e. 
<i>R^2 = (y_1 - y_2)^2 + (phi_1 - phi_2)^2</i> and equivalent. 
   
 
<a name="anchor46"></a>
<p/><strong>friend bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New) &nbsp;</strong> <br/>
transfer four-momentum between the two four-vectors so that they get 
the masses <code>m1New</code> and <code>m2New</code>, respectively. 
Note that <code>p1Move</code> and <code>p2Move</code> act both as 
input and output arguments. The method will return false if the invariant 
mass of the four-vectors is too small to accommodate the new masses, 
and then the four-vectors are not changed. 
   
 
<a name="anchor47"></a>
<p/><strong>friend pair&lt;Vec4,Vec4&gt; getTwoPerpendicular(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
create a pair of four-vectors that are perpendicular to both input vectors 
and to each other, and have the squared norm <i>-1</i>. 
   
 
<h3>Operations with four-vectors</h3> 
 
Of course one should be able to add, subtract and scale four-vectors, 
and more: 
 
<a name="anchor48"></a>
<p/><strong>Vec4 Vec4::operator-() &nbsp;</strong> <br/>
return a vector with flipped sign for all components, while leaving 
the original vector unchanged. 
   
 
<a name="anchor49"></a>
<p/><strong>Vec4& Vec4::operator+=(const Vec4& v) &nbsp;</strong> <br/>
add a four-vector to an existing one. 
   
 
<a name="anchor50"></a>
<p/><strong>Vec4& Vec4::operator-=(const Vec4& v) &nbsp;</strong> <br/>
subtract a four-vector from an existing one. 
   
 
<a name="anchor51"></a>
<p/><strong>Vec4& Vec4::operator*=(double f) &nbsp;</strong> <br/>
multiply all four-vector components by a real number. 
   
 
<a name="anchor52"></a>
<p/><strong>Vec4& Vec4::operator/=(double f) &nbsp;</strong> <br/>
divide all four-vector components by a real number. 
   
 
<a name="anchor53"></a>
<p/><strong>friend Vec4 operator+(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
add two four-vectors. 
   
 
<a name="anchor54"></a>
<p/><strong>friend Vec4 operator-(const Vec4& v1, const Vec4& v2) &nbsp;</strong> <br/>
subtract two four-vectors. 
   
 
<a name="anchor55"></a>
<p/><strong>friend Vec4 operator*(double f, const Vec4& v) &nbsp;</strong> <br/>
   
<a name="anchor56"></a>
<strong>friend Vec4 operator*(const Vec4& v, double f) &nbsp;</strong> <br/>
multiply a four-vector by a real number. 
   
 
<a name="anchor57"></a>
<p/><strong>friend Vec4 operator/(const Vec4& v, double f) &nbsp;</strong> <br/>
divide a four-vector by a real number. 
   
 
<a name="anchor58"></a>
<p/><strong>friend double operator*(const Vec4& v1, const Vec4 v2) &nbsp;</strong> <br/>
four-vector product. 
   
 
<p/> 
There are also a few related operations that are normal member methods: 
 
<a name="anchor59"></a>
<p/><strong>void Vec4::rescale3(double f) &nbsp;</strong> <br/>
multiply the three-vector components by <i>f</i>, but keep the 
fourth component unchanged. 
   
 
<a name="anchor60"></a>
<p/><strong>void Vec4::rescale4(double f) &nbsp;</strong> <br/>
multiply all four-vector components by <i>f</i>. 
   
 
<a name="anchor61"></a>
<p/><strong>void Vec4::flip3() &nbsp;</strong> <br/>
flip the sign of the three-vector components, but keep the 
fourth component unchanged. 
   
 
<a name="anchor62"></a>
<p/><strong>void Vec4::flip4() &nbsp;</strong> <br/>
flip the sign of all four-vector components. 
   
 
<h3>Rotations and boosts</h3> 
 
A common task is to rotate or boost four-vectors. In case only one 
four-vector is affected the operation may be performed directly on it. 
However, in case many particles are affected, the helper class 
<code>RotBstMatrix</code> can be used to speed up operations. 
 
<a name="anchor63"></a>
<p/><strong>void Vec4::rot(double theta, double phi) &nbsp;</strong> <br/>
rotate the three-momentum with the polar angle <i>theta</i> 
and the azimuthal angle <i>phi</i>. 
   
 
<a name="anchor64"></a>
<p/><strong>void Vec4::rotaxis(double phi, double nx, double ny, double nz) &nbsp;</strong> <br/>
rotate the three-momentum with the azimuthal angle <i>phi</i> 
around the direction defined by the <i>(n_x, n_y, n_z)</i> 
three-vector. 
   
 
<a name="anchor65"></a>
<p/><strong>void Vec4::rotaxis(double phi, Vec4& n) &nbsp;</strong> <br/>
rotate the three-momentum with the azimuthal angle <i>phi</i> 
around the direction defined by the three-vector part of <i>n</i>. 
   
 
<a name="anchor66"></a>
<p/><strong>void Vec4::bst(double betaX, double betaY, double betaZ) &nbsp;</strong> <br/>
boost the four-momentum by <i>beta = (beta_x, beta_y, beta_z)</i>. 
   
 
<a name="anchor67"></a>
<p/><strong>void Vec4::bst(double betaX, double betaY, double betaZ, double gamma) &nbsp;</strong> <br/>
boost the four-momentum by <i>beta = (beta_x, beta_y, beta_z)</i>, 
where the <i>gamma = 1/sqrt(1 - beta^2)</i> is also input to allow 
better precision when <i>beta</i> is close to unity. 
   
 
<a name="anchor68"></a>
<p/><strong>void Vec4::bst(const Vec4& p) &nbsp;</strong> <br/>
boost the four-momentum by <i>beta = (p_x/E, p_y/E, p_z/E)</i>. 
   
 
<a name="anchor69"></a>
<p/><strong>void Vec4::bst(const Vec4& p, double m) &nbsp;</strong> <br/>
boost the four-momentum by <i>beta = (p_x/E, p_y/E, p_z/E)</i>, 
where the <i>gamma = E/m</i> is also calculated from input to allow 
better precision when <i>beta</i> is close to unity. 
   
 
<a name="anchor70"></a>
<p/><strong>void Vec4::bstback(const Vec4& p) &nbsp;</strong> <br/>
boost the four-momentum by <i>beta = (-p_x/E, -p_y/E, -p_z/E)</i>. 
   
 
<a name="anchor71"></a>
<p/><strong>void Vec4::bstback(const Vec4& p, double m) &nbsp;</strong> <br/>
boost the four-momentum by <i>beta = (-p_x/E, -p_y/E, -p_z/E)</i>, 
where the <i>gamma = E/m</i> is also calculated from input to allow 
better precision when <i>beta</i> is close to unity. 
   
 
<a name="anchor72"></a>
<p/><strong>void Vec4::rotbst(const RotBstMatrix& M) &nbsp;</strong> <br/>
perform a combined rotation and boost; see below for a description 
of the <code>RotBstMatrix</code>. 
   
 
<p/> 
For a longer sequence of rotations and boosts, and where several 
<code>Vec4</code> are to be rotated and boosted in the same way, 
a more efficient approach is to define a <code>RotBstMatrix</code>, 
which forms a separate auxiliary class. You can build up this 
4-by-4 matrix by successive calls to the methods of the class, 
such that the matrix encodes the full sequence of operations. 
The order in which you do these calls must agree with the imagined 
order in which the rotations/boosts should be applied to a 
four-momentum, since in general the operations do not commute. 
<br/>(Mathematically you would e.g. define <i>M = M_3 M_2 M_1</i> 
in that <i>M p = M_3( M_2( M_1 p) ) )</i>. That is, operations 
on the four-vector <i>p</i> are carried out in the order first 
<i>M_1</i>, then <i>M_2</i> and finally <i>M_3</i>. Thus 
<i>M_1, M_2, M_3</i> is also the order in which you should input 
rotations and boosts to <i>M</i>.) 
 
<a name="anchor73"></a>
<p/><strong>RotBstMatrix::RotBstMatrix() &nbsp;</strong> <br/>
creates a diagonal unit matrix, i.e. one that leaves a four-vector 
unchanged. 
   
 
<a name="anchor74"></a>
<p/><strong>RotBstMatrix::RotBstMatrix(const RotBstMatrix& Min) &nbsp;</strong> <br/>
creates a copy of the input matrix. 
   
 
<a name="anchor75"></a>
<p/><strong>RotBstMatrix& RotBstMatrix::operator=(const RotBstMatrix4& Min) &nbsp;</strong> <br/>
copies the input matrix. 
   
 
<a name="anchor76"></a>
<p/><strong>void RotBstMatrix::rot(double theta = 0., double phi = 0.) &nbsp;</strong> <br/>
rotate by this polar and azimuthal angle. 
   
 
<a name="anchor77"></a>
<p/><strong>void RotBstMatrix::rot(const Vec4& p) &nbsp;</strong> <br/>
rotate so that a vector originally along the <i>+z</i> axis becomes 
parallel with <i>p</i>. More specifically, rotate by <i>-phi</i>, 
<i>theta</i> and <i>phi</i>, with angles defined by <i>p</i>. 
   
 
<a name="anchor78"></a>
<p/><strong>void RotBstMatrix::bst(double betaX = 0., double betaY = 0., double betaZ = 0.) &nbsp;</strong> <br/>
boost by this <i>beta</i> vector. 
   
 
<a name="anchor79"></a>
<p/><strong>void RotBstMatrix::bst(const Vec4&) &nbsp;</strong> <br/>
   
<a name="anchor80"></a>
<strong>void RotBstMatrix::bstback(const Vec4&) &nbsp;</strong> <br/>
boost with a <i>beta = p/E</i> or <i>beta = -p/E</i>, respectively. 
   
 
<a name="anchor81"></a>
<p/><strong>void RotBstMatrix::bst(const Vec4& p1, const Vec4& p2) &nbsp;</strong> <br/>
boost so that <i>p_1</i> is transformed to <i>p_2</i>. It is assumed 
that the two vectors obey <i>p_1^2 = p_2^2</i>. 
   
 
<a name="anchor82"></a>
<p/><strong>void RotBstMatrix::toCMframe(const Vec4& p1, const Vec4& p2) &nbsp;</strong> <br/>
boost and rotate to the rest frame of <i>p_1</i> and <i>p_2</i>, 
with <i>p_1</i> along the <i>+z</i> axis. 
   
 
<a name="anchor83"></a>
<p/><strong>void RotBstMatrix::fromCMframe(const Vec4& p1, const Vec4& p2) &nbsp;</strong> <br/>
rotate and boost from the rest frame of <i>p_1</i> and <i>p_2</i>, 
with <i>p_1</i> along the <i>+z</i> axis, to the actual frame of 
<i>p_1</i> and <i>p_2</i>, i.e. the inverse of the above. 
   
 
<a name="anchor84"></a>
<p/><strong>void RotBstMatrix::rotbst(const RotBstMatrix& Min); &nbsp;</strong> <br/>
combine the current matrix with another one. 
   
 
<a name="anchor85"></a>
<p/><strong>void RotBstMatrix::invert() &nbsp;</strong> <br/>
invert the matrix, which corresponds to an opposite sequence and sign 
of rotations and boosts. 
   
 
<a name="anchor86"></a>
<p/><strong>void RotBstMatrix::reset() &nbsp;</strong> <br/>
reset to no rotation/boost; i.e. the default at creation. 
   
 
<a name="anchor87"></a>
<p/><strong>double RotBstMatrix::value(int i, int j) &nbsp;</strong> <br/>
return the value of the (i,j) = [i][j] matrix element. 
   
 
<a name="anchor88"></a>
<p/><strong>double RotBstMatrix::deviation() &nbsp;</strong> <br/>
crude estimate how much a matrix deviates from the unit matrix: 
the sum of the absolute values of all non-diagonal matrix elements 
plus the sum of the absolute deviation of the diagonal matrix 
elements from unity. 
   
 
<a name="anchor89"></a>
<p/><strong>friend ostream& operator&lt;&lt;(ostream&, const RotBstMatrix& M) &nbsp;</strong> <br/>
writes out the values of the sixteen components of a 
<code>RotBstMatrix</code>, on four consecutive lines and 
ended with a newline. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
