#include<tuple>
#include<vector>

typedef std::tuple<float, float, float> float3;

typedef struct {
    float energy_density;
    float temperature;
    float entropy_density;
    float qgp_fraction;
    float vx, vy, vz;
    // do we need pi^{mu nu}, Pi, net_baryon, net_charge?
    // for thermal photon or other kinds of studies
} BulkElement;

class EvolutionHistory{
  public:
    float tmin, nt, dt;
    float xmin, nx, dx;
    float ymin, ny, dy;
    float zmin, nz, dz;
    // default: set using_tz_for_tau_eta=true
    bool  using_tz_for_tau_eta;
    // the bulk information
    std::vector<BulkElement> data;

    EvolutionHistory();
};


class FluidDynamics{
  public:
    // How to store this data? In memory or hard disk?
    // 3D hydro may eat out the memory, 
    // for large dataset, std::deque is better than std::vector.
    EvolutionHistory bulk_info;

    /*Keep this interface open in the beginning.*/
    //FreezeOutHyperSf hyper_sf;
    
    /* currently we have no standard for passing configurations */
    // pure virtual function; to be implemented by users
    // should make it easy to save evolution history to bulk_info
    virtual void evolution(const EnergyMomentumTensor & tmn,
                           float xmax, float ymax, float hmax,
                           float tau0, float dtau, float dx,
                           float dy, float dh, float etaos,
                           float dtau_out, float dx_out, float dy_out,
                           float dh_out) const=0;

    // the following functions should be implemented in Jetscape
    float get_energy_density(float time, float x, float y, float z);
    float get_entropy_density(float time, float x, float y, float z);
    float get_temperature(float time, float x, float y, float z);
    float get_qgp_fraction(float time, float x, float y, float z);
    // float3 return std::make_tuple(vx, vy, vz)
    float3 get_3fluid_velocity(float time, float x, float y, float z);
    // float4 return std::make_tuple(ut, ux, uy, uz)
    float4 get_4fluid_velocity(float time, float x, float y, float z);
    //float get_net_baryon_density(float time, float x, float y, float z);
    //float get_net_charge_density(float time, float x, float y, float z);
};
