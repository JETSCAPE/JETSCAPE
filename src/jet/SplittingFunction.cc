#include "SplittingFunction.h"
#include "Matter.h"

namespace SplittingFunction{

//need to put in the matter_obj.length information
double P_z_gg_int(double cg, double cg1, double loc_e, double cg3,double l_fac, double E2, Matter& matter_obj){

  double t3, t4, t5, t10, t11, t12, t15, t9, qL, tau, res, limit_factor, lz, uz, m_fac;

  //if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
  //JSINFO<<"P_z_gg_int cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" cg3 "<<cg3<<" l_fac "<<l_fac<<" E2 "<<E2;
  t3 = std::log((1.0 - cg1));
  t4 = std::log(cg1);
  t5 = std::pow(cg1, 2);
  t10 = std::log((1.0 - cg));
  t11 = std::log(cg);
  t12 = std::pow(cg, 2);
  t15 = -(2.0 * cg1) - t3 + t4 - (2.0 / 3.0) * t5 * cg1 + t5 + (2.0 * cg) +
        t10 - t11 + (2.0 / 3.0) * t12 * cg - t12;

  res = t15;

  limit_factor = 2.0 * std::sqrt(2.0) * cg3 / E2 / 0.1;

  if (limit_factor < 0.0) {
    cerr << " error in z limit factor for medium calculation = " << limit_factor
         << endl;
    throw std::runtime_error(" error in z limit factor for medium calculation");
  }

  tau = l_fac;

  if ((matter_obj.length - loc_e) < tau)
    tau = (matter_obj.length - loc_e);

  //             if ((matter_obj.length - loc_e) < tau) tau = 0;

  if (loc_e > matter_obj.length)
    tau = 0.0;

  m_fac = 1.0;

  //            if ((qhat*tau<1.0)&&(in_vac==false)) m_fac = 1.0/qhat/tau;

  // SC
  //qL = m_fac*qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  t9 = 2.0 * cg1 - 1.0 / (-1.0 + cg1) - 1.0 / cg1 - 2.0 * cg +
       1.0 / (-1.0 + cg) + 1.0 / cg;

  if (t9 < 0.0) {
    cerr << "ERROR: medium contribution negative in P_z_gg_int : t9 = " << t9
         << endl;

    cerr << " cg, cg1 = " << cg << "  " << cg1 << endl;
    throw std::runtime_error(
        "ERROR: medium contribution negative in P_z_gg_int");
  }

  res = t15 + 2.0 * t9 * qL / cg3;

  return (res);
}

//need to put in matter_obj.length information
double P_z_qq_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2,Matter& matter_obj){
  double t_q1, t_q3, t_q4, t_q6, t_q8, t_q9, t_q12, q_q1, q_q4, q_q6, q_q9,
      q_q11, qL, tau, res;

  //JSINFO<<"P_z_qq_int cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" cg3 "<<cg3<<" l_fac "<<l_fac<<" E2 "<<E2;

  if ((cg < cg1 / (2.0 * E2 * E2 / cg1 + 1.0)))
    cg = cg1 / (2.0 * E2 * E2 / cg1 + 1.0);

  t_q1 = std::pow(cg1, 2);
  t_q3 = 1.0 - cg1;
  t_q4 = t_q3 * t_q3;
  t_q6 = std::pow(cg, 2);
  t_q8 = 1.0 - cg;
  t_q9 = t_q8 * t_q8;
  t_q12 = t_q1 * cg1 / 6.0 - t_q4 * t_q3 / 6.0 - t_q6 * cg / 6.0 +
          t_q9 * t_q8 / 6.0;

  tau = l_fac;

  if ((matter_obj.length - loc_e) < tau)
    tau = (matter_obj.length - loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  q_q1 = std::log(cg1);
  q_q4 = std::log(1.0 - cg1);
  q_q6 = std::log(cg);
  q_q9 = std::log(1.0 - cg);
  q_q11 = -cg1 + q_q1 / 2.0 - q_q4 / 2.0 + cg - q_q6 / 2.0 + q_q9 / 2.0;

  if (q_q11 < 0.0) {
    cerr << "ERROR: medium contribution negative in P_z_gg_int : q_q11 = "
         << q_q11 << endl;
    throw std::runtime_error(
        "ERROR: medium contribution negative in P_z_gg_int");
  }

  res = t_q12 * Tf / Ca + 2.0 * qL * q_q11 / cg3 * (Tf * Cf / Ca / Ca);

  return (res);
}

double P_z_qq_int_w_M_vac_only(double M, double cg, double cg1,double loc_e, double cg3, double l_fac,double E2,Matter& matter_obj) {
  double t_q1, t_q3, t_q4, t_q6, t_q8, t_q9, t_q12, q_q1, q_q4, q_q6, q_q9,q_q11, qL, tau, res;

  //JSINFO<<"P_z_qq_int_w_M_vac_only cg "<<cg<<" M "<<M<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" cg3 "<<cg3<<" l_fac "<<l_fac<<" E2 "<<E2;

  if ((cg < cg1 / (2.0 * E2 * E2 / cg1 + 1.0)))
    cg = cg1 / (2.0 * E2 * E2 / cg1 + 1.0);

  t_q1 = std::pow(cg1, 2);
  t_q3 = 1.0 - cg1;
  t_q4 = t_q3 * t_q3;
  t_q6 = std::pow(cg, 2);
  t_q8 = 1.0 - cg;
  t_q9 = t_q8 * t_q8;
  t_q12 = t_q1 * cg1 / 6.0 - t_q4 * t_q3 / 6.0 - t_q6 * cg / 6.0 +
          t_q9 * t_q8 / 6.0;

  tau = l_fac;

  if ((matter_obj.length - loc_e) < tau)
    tau = (matter_obj.length - loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  q_q1 = std::log(cg1);
  q_q4 = std::log(1.0 - cg1);
  q_q6 = std::log(cg);
  q_q9 = std::log(1.0 - cg);
  q_q11 = -cg1 + q_q1 / 2.0 - q_q4 / 2.0 + cg - q_q6 / 2.0 + q_q9 / 2.0;

  if (q_q11 < 0.0) {
    cerr << "ERROR: medium contribution negative in P_z_gg_int_w_M : q_q11 = "
         << q_q11 << endl;
    cout << " z_low = " << cg << " z_hi = " << cg1 << endl;
    throw std::runtime_error(
        "ERROR: medium contribution negative in P_z_gg_int");
  }

  res = t_q12 * Tf / Ca + 2.0 * qL * q_q11 / cg3 * (Tf * Cf / Ca / Ca);

  return (res);
}

double P_z_qp_int(double cg, double cg1, double loc_e, double cg3,double l_fac, double E2,Matter& matter_obj) {

  //JSINFO<<"P_z_qp_int cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" cg3 "<<cg3<<" l_fac "<<l_fac<<" E2 "<<E2;

  double t2, t5, t7, t10, t12, q2, q6, q10, tau, qL, res;

  if ((cg < cg1 / (2.0 * E2 * E2 / cg1 + 1.0)))
    cg = cg1 / (2.0 * E2 * E2 / cg1 + 1.0);

  t2 = std::pow(cg1, 2);
  t5 = std::log(1.0 - cg1);
  t7 = std::pow(cg, 2);
  t10 = std::log(1.0 - cg);
  t12 = -cg1 - t2 / 2.0 - 2.0 * t5 + cg + t7 / 2.0 + 2.0 * t10;

  //    return(t12);

  q10 = 0.0;
  tau = l_fac;

  if ((matter_obj.length - loc_e) < tau)
    tau = (matter_obj.length - loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  res = t12 + 2.0 * qL * q10 / cg3;

  return (res);
}

double P_z_qg_int(double cg, double cg1, double loc_e, double cg3,double l_fac, double E2,Matter& matter_obj) {

  double t2, t5, t7, t10, t12, q2, q6, q10, tau, qL, res;

  //JSINFO<<"P_z_qg_int cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" cg3 "<<cg3<<" l_fac "<<l_fac<<" E2 "<<E2;

  if ((cg < cg1 / (2.0 * E2 * E2 / cg1 + 1.0)))
    cg = cg1 / (2.0 * E2 * E2 / cg1 + 1.0);

  t2 = std::pow(cg1, 2);
  t5 = std::log(1.0 - cg1);
  t7 = std::pow(cg, 2);
  t10 = std::log(1.0 - cg);
  t12 = -cg1 - t2 / 2.0 - 2.0 * t5 + cg + t7 / 2.0 + 2.0 * t10;

  //	return(t12);

  q2 = std::log(cg1);
  q6 = std::log(cg);
  q10 = q2 - 2.0 / (cg1 - 1.0) - q6 + 2.0 / (cg - 1.0);

  tau = l_fac;

  if ((matter_obj.length - loc_e) < tau)
    tau = (matter_obj.length - loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  res = t12 + 2.0 * qL * q10 / cg3;

  return (res);
}

double P_z_qg_int_w_M(double M, double cg, double cg1, double loc_e,double cg3, double l_fac, double E2,Matter& matter_obj) {

  double t2, t5, t7, t10, t12, tau, qL, res;
  //JSINFO<<"P_z_qg_int_w_M cg "<<cg<<" M "<<M<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" cg3 "<<cg3<<" l_fac "<<l_fac<<" E2 "<<E2;

  if (std::abs(cg - cg1) < rounding_error)
    return (cg);
  //if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );

  t2 = std::pow(cg1, 2);
  t5 = std::log(1.0 - cg1);
  t7 = std::pow(cg, 2);
  t10 = std::log(1.0 - cg);
  t12 = -cg1 - t2 / 2.0 - 2.0 * t5 + cg + t7 / 2.0 + 2.0 * t10;

  //	return(t12);

  double q1 = M * M;
  double q2 = 1.0 / cg3;
  double q3 = q2 * q1;
  double q5 = 4.0 * q3 + 1.0;
  double q6 = 1.0 - cg1;
  double q7 = std::log(q6);
  double q9 = q1 * q1;
  double q10 = cg3 * cg3;
  double q12 = 1.0 / q10 * q9;
  double q14 = q12 - 2.0 * q3 + 1.0 / 2.0;
  double q15 = std::log(cg1);
  double q17 = q3 + 1.0;
  double q26 = std::pow(cg1, 2.0);
  double q30 = 1.0 - cg;
  double q31 = std::log(q30);
  double q33 = std::log(cg);
  double q43 = std::pow(cg, 2.0);
  double q47 = q7 * q5 + q15 * q14 + cg1 * q17 / 2.0 + 2.0 / q6 +
               3.0 / 2.0 * q2 / cg1 * q1 - 1.0 / q26 * q12 / 2.0 - q31 * q5 -
               q33 * q14 - cg * q17 / 2.0 - 2.0 / q30 -
               3.0 / 2.0 * q2 / cg * q1 + 1.0 / q43 * q12 / 2.0;

  tau = l_fac;

  if ((matter_obj.length - loc_e) < tau)
    tau = (matter_obj.length - loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * 2.0 * tau;
  }

  double e1 = M * M;
  double e3 = 1.0 / cg3 * e1;
  double e4 = 1.0 - cg1;
  double e5 = std::log(e4);
  double e10 = 1.0 - cg;
  double e11 = std::log(e10);
  double e17 = 2.0 * (e5 + 1.0 / e4 + cg1 / 2.0) * e3 -
               2.0 * (e11 + 1.0 / e10 + cg / 2.0) * e3;

  double eL;

  if (tau < rounding_error) {
    eL = 0.0;
  } else {
    matter_obj.ehat = 0.0; //fncAvrEhat(loc_e,tau);
    eL = matter_obj.ehat * 4.0;
  }

  double f1 = M * M;
  double f2 = 1.0 / cg3;
  double f3 = f2 * f1;
  double f4 = 13.0 * f3;
  double f6 = f1 * f1;
  double f7 = f6 * (f4 + 15.0);
  double f8 = 1.0 - cg1;
  double f9 = std::log(f8);
  double f10 = cg3 * cg3;
  double f11 = 1.0 / f10;
  double f15 = std::log(cg1);
  double f23 = f1 * (39.0 / 4.0 * f11 * f6 + 15.0 / 2.0 * f3 + 1.0);
  double f28 = f6 * (f4 + 15.0 / 2.0);
  double f29 = f8 * f8;
  double f37 = 1.0 / f10 / cg3 * f6 * f1;
  double f42 = 1.0 - cg;
  double f43 = std::log(f42);
  double f47 = std::log(cg);
  double f54 = f42 * f42;
  double f63 = f11 * f9 * f7 / 4.0 + f2 * f1 * f15 / 2.0 + 1.0 / f8 * f2 * f23 -
               1.0 / f29 * f11 * f28 / 2.0 + 13.0 / 6.0 / f29 / f8 * f37 -
               f11 * f43 * f7 / 4.0 - f2 * f1 * f47 / 2.0 -
               1.0 / f42 * f2 * f23 + 1.0 / f54 * f11 * f28 / 2.0 -
               13.0 / 6.0 / f54 / f42 * f37;

  double e2L;

  if (tau < rounding_error) {
    e2L = 0.0;
  } else {
    matter_obj.e2hat = matter_obj.qhat / 2.0; //fncAvrE2hat(loc_e,tau);
    e2L = matter_obj.e2hat * 8.0 / (tau * cg3);
  }

  res = t12 + qL * q47 / cg3;
  //+ eL*e17/cg3 + e2L*f63/cg3;
  // Uncomment only if you have an eL larger than 2 times e2L for charm, and derive expression for bottom.
  // MC simulation is not valid for all choices of e-hat and e2-hat.

  return (res);
}

}