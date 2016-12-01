/*******************************************************************************

This is the class that holds the configuration information and parameters to run
Jetscape. This class will be built up as different modules require different 
parameters. Members set to public since modules will be reading these parameters

*******************************************************************************/

class Jet {
 public :
  Jet(float p_in[4]);
  void set_jet_p(float p_in[4]);
  ~Jet(){};
  float[] get_jet_p() { return jet_p; ;
  float get_jet_pt() { return jet_pt; };
  float get_jet_eta() { return jet_eta; ;
  float get_jet_phi() { return jet_phi; ;
 private :
  float jet_p[4]; // momenta of jet
  float jet_pt;   // pt of jet
  float jet_eta;  // eta of jet
  float jet_phi;  // phi of jet
  // may be useful to keep track of e loss model applied to a particular jet
  string e_loss_model = "";
  void setPtEtaPhi();
}

Jet::Jet(float p_in[4]) {
  jet_p[0]=p_in[0];
  jet_p[1]=p_in[1];
  jet_p[2]=p_in[2];
  jet_p[3]=p_in[3];
  setPtEtaPhi();
}

void setPtEtaPhi() {
  // to be implemented, set pt eta phi from jet_p four momentum vector
}
