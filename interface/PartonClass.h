/*******************************************************************************

This is the data structure class for Partons. This class should evolve as new
needs arsise.

*******************************************************************************/

class Parton {
 public :
  Parton(int label, int id, int stat, float p_in[4], float x[4]);
  ~Parton(){};
  void set_label(int label) { plabel = label; };
  void set_id(int id) { pid = id; };
  void set_stat(int stat) { pstat = stat; };
  void set_p(float p_in[4]) {};
  void set_x(float x[4]) {};
  void set_pl() {};
  void set_mean_form_time() {};


 private :
  int pid;                  // particle id
  int pstat;                // status of particle
  int plabel;               // line number in event record
  float mass;               // particle mass
  float e;                  // particle energy
  float pj[4];              // particle momenta in jet frame
  float pp, pm, t;          // light cone momenta
  float pt, eta, phi;       // pt, pseudorapidity, azimuthal angle
  float p[4];               // particle momenta in CM frame
  float xp[4];              // particle position in CM frame
  float pl;                 // modulus of 3-momentum
  float mean_form_time;     // mean formation time
  float form_time;          // event by event formation time

}

Parton::Parton(int label, int id, int stat, float p_in[4], float x[4]) {
  p[0]=p_in[0];
  p[1]=p_in[1];
  p[2]=p_in[2];
  p[3]=p_in[3];
  setPtEtaPhi();

}

void setPtEtaPhi() {
  // to be implemented, set pt eta phi from p four momentum vector
}
