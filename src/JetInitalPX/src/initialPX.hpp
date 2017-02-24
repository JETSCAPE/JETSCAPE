
class initialPX {

      public:

          void next();


          double px() { return px0_; }
          double py() { return py0_; }
          double pz() { return pz0_; }
          double e() { return e0_; }
          double x() { return x0_; }
          double y() { return y0_; }
          double z() { return z0_; }
          double t() { return t0_; }
          double m() { return m0_; }
          int id() { return id0_; }

          void set_px(double ppx) { px0_ = ppx; }
          void set_py(double ppy) { py0_ = ppy; }
          void set_pz(double ppz) { pz0_ = ppz; }
          void set_e(double pe) { e0_ = pe; }
          void set_x(double prx) { x0_ = prx; }
          void set_y(double pry) { y0_ = pry; }
          void set_z(double prz) { z0_ = prz; }
          void set_t(double pt) { t0_ = pt; }
          void set_m(double pmass) { m0_ = pmass; }
          void set_id(int pid) { id0_ = pid; }

      private:

          int id0_;
          double px0_;
          double py0_;
          double pz0_;
          double e0_;
          double x0_;
          double y0_;
          double z0_;
          double t0_;
          double m0_;

};

