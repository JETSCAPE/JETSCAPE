#include "SudakovFunction.h"
#include "Matter.h"

namespace SudakovFunction{

double sudakov_Pgg(double g0, double g1, double loc_c, double E, Matter& matter_obj) {
  double sud, g;
  int blurb;
  //JSINFO<<"sudakov_Pgg g0 "<<g0<<" g1 "<<g1<<" loc_c "<<loc_c<<" E "<<E;
  sud = 1.0;

  if (g1 < 2.0 * g0) {
    cerr << " warning: the lower limit of the sudakov > 1/2 upper limit, "
            "returning 1 "
         << endl;
    cerr << " in sudakov_P glue glue, g0, g1 = " << g0 << "  " << g1 << endl;
    throw std::runtime_error(" warning: the lower limit of the sudakov > 1/2 "
                             "upper limit, returning 1");

    return (sud);
  }
  g = 2.0 * g0;

  if (g1 > g) {

    sud = exp(-1.0 * (Ca / 2.0 / pi) * sud_val_GG(g0, g, g1, loc_c, E,matter_obj));
  }
  return (sud);
}

double sud_val_GG(double h0, double h1, double h2, double loc_d,double E1, Matter& matter_obj) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  //JSINFO<<"sud_val_GG h0 "<<h0<<" h1 "<<h1<<" h2 "<<h2<<" loc_d "<<loc_d<<" E1 "<<E1;
  val = 0.0;

  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  t_form = 2.0 * E1 / h;

  val = matter_obj.alpha_s(h) * sud_z_GG(h0, h, loc_d, t_form, E1,matter_obj);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  t_form = 2.0 * E1 / hL;

  intg_L = matter_obj.alpha_s(hL) * sud_z_GG(h0, hL, loc_d, t_form, E1,matter_obj) * (h - h1);

  hR = (h + h2) / 2.0;

  t_form = 2.0 * E1 / hR;

  intg_R = matter_obj.alpha_s(hR) * sud_z_GG(h0, hR, loc_d, t_form, E1,matter_obj) * (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;

  if ((diff > approx) || (span > error)) {
    intg = sud_val_GG(h0, h1, h, loc_d, E1,matter_obj) + sud_val_GG(h0, h, h2, loc_d, E1,matter_obj);
  }

  //	cout << " returning with intg = " << intg << endl;

  return (intg);
}

double sud_z_GG(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj) {

  double t2, t3, t7, t11, t12, t15, t21, t25, q2, q3, q8, q12, qL, tau, res,
      z_min, limit_factor, lz, uz, mz, m_fac;
  double t_q1, t_q3, t_q4, t_q6, t_q8, t_q9, t_q12, q_q1, q_q4, q_q6, q_q9,
      q_q11;
  //JSINFO<<"sud_z_GG cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" l_fac "<<l_fac<<" E2 "<<E2;
  z_min = std::sqrt(2) * E_minimum / E2;

  if (cg1 < 2.0 * cg) {

    //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
    return (0.0);
  };

  t2 = std::pow(cg1, 2);
  t3 = t2 * cg1;
  t7 = std::log(cg);
  t11 = std::abs(cg - cg1);
  t12 = std::log(t11);
  t15 = std::pow(cg, 2);
  t21 = t2 * t2;
  t25 = -(5.0 * t3 - 12.0 * cg * t2 + 6.0 * t7 * t3 - 6.0 * t12 * t3 -
          4.0 * t15 * cg + 6.0 * t15 * cg1) /
        t21 / 3.0;

  res = t25;

  limit_factor = 2.0 * std::sqrt(2.0) * cg1 / E2 / 0.1;

  if (limit_factor < 0.0) {
    cerr << " error in z limit factor for medium calculation in sud-z-gg = "
         << limit_factor << endl;
    throw std::runtime_error(
        "error in z limit factor for medium calculation in sud-z-gg");
  }

  q2 = 1.0 / cg1;
  q3 = cg * q2;
  q8 = 1.0 - q3;
  q12 = (2.0 - 4.0 * q3 + 2.0 / cg * cg1 - 2.0 / q8) * q2;

  if (q12 < 0.0) {
    cerr << "ERROR: medium contribution negative in sud_z_GG: q12 = " << q12
         << endl;
    cerr << "cg, cg1 = " << cg << "  " << cg1 << endl;
    cerr << " t25 = " << t25 << endl;
    throw std::runtime_error("ERROR: medium contribution negative in sud_z_GG");
  }

  tau = l_fac;

  if ((matter_obj.length- loc_e) < tau)
    tau = (matter_obj.length- loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  m_fac = 1.0;

  // SC
  //qL = m_fac*matter_obj.qhat*0.6*tau*profile(loc_e+tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  res = t25 + 2.0 * qL * q12 / cg1;
  //        }
  //        else{
  //            cout << " z trap for medium enabled in sud-val-z " << endl ;
  //        }
  //    }

  return (res);
}

double sudakov_Pqq(double q0, double q1, double loc_c, double E, Matter& matter_obj) {
  double sud, q;
  //JSINFO<<"sudakov_Pqq q0 "<<q0<<" q1 "<<q1<<" loc_c "<<loc_c<<" E "<<E;
  sud = 1.0;

  if (q1 < 2.0 * q0)
  //	if (g1<g0)
  {
    JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, "
              "returning 1 ";
    JSWARN << " in sudakov_Pquark quark, q0, q1 = " << q0 << "  " << q1;
    return (sud);
  }
  q = 2.0 * q0;

  //	g = g0 ;

  sud = exp(-1.0 * (Tf / 2.0 / pi) * sud_val_QQ(q0, q, q1, loc_c, E,matter_obj));

  return (sud);
}

double sudakov_Pqq_w_M_vac_only(double M, double q0, double q1,double loc_c, double E, Matter& matter_obj) {
  double sud, q;
  //JSINFO<<"sudakov_Pqq_w_M_vac_only q0 "<<q0<<" q1 "<<q1<<" loc_c "<<loc_c<<" E "<<E;
  sud = 1.0;

  if (q1 < 2.0 * (q0 + M * M)) {
    JSWARN << " warning: the upper limit of the sudakov q1<2.0*(q0+M*M), "
              "returning 1 ";
    JSWARN << " in sudakov_Pquark quark, q0, q1 = " << q0 << "  " << q1;
    return (sud);
  } else {
    q = 2.0 * (q0 + M * M);
    sud = exp(-1.0 * (Tf / 2.0 / pi) *
              sud_val_QQ_w_M_vac_only(M, q0, q, q1, loc_c, E,matter_obj));
    return (sud);
  }
}

double sud_val_QQ(double h0, double h1, double h2, double loc_d,double E1, Matter& matter_obj) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  //JSINFO<<"sud_val_QQ h0 "<<h0<<" h1 "<<h1<<" h2 "<<h2<<" loc_d "<<loc_d<<" E1 "<<E1;
  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  t_form = 2.0 * E1 / h;

  val = matter_obj.alpha_s(h) * sud_z_QQ(h0, h, loc_d, t_form, E1,matter_obj);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  t_form = 2.0 * E1 / hL;

  intg_L = matter_obj.alpha_s(hL) * sud_z_QQ(h0, hL, loc_d, t_form, E1,matter_obj) * (h - h1);

  hR = (h + h2) / 2.0;

  t_form = 2.0 * E1 / hR;

  intg_R = matter_obj.alpha_s(hR) * sud_z_QQ(h0, hR, loc_d, t_form, E1,matter_obj) * (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  //JSINFO << BOLDYELLOW << " h2,  h1, diff = " << h2 << " , " << h1 << " , " << diff  ;
  //JSINFO << BOLDYELLOW << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R ;

  if ((diff > approx) || (span > error)) {
    intg = sud_val_QQ(h0, h1, h, loc_d, E1,matter_obj) + sud_val_QQ(h0, h, h2, loc_d, E1,matter_obj);
  }

  //	cout << " returning with intg = " << intg << endl;

  return (intg);
}

double sud_val_QQ_w_M_vac_only(double M, double h0, double h1,double h2, double loc_d, double E1, Matter& matter_obj) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  //JSINFO<<"sud_val_QQ_w_M_vac_only h0 "<<h0<<" h1 "<<h1<<" h2 "<<h2<<" loc_d "<<loc_d<<" E1 "<<E1;
  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  t_form = 2.0 * E1 / h;

  val = matter_obj.alpha_s(h) * sud_z_QQ_w_M_vac_only(M, h0, h, loc_d, t_form, E1,matter_obj);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  t_form = 2.0 * E1 / hL;

  intg_L = matter_obj.alpha_s(hL) * sud_z_QQ_w_M_vac_only(M, h0, hL, loc_d, t_form, E1,matter_obj) *
           (h - h1);

  hR = (h + h2) / 2.0;

  t_form = 2.0 * E1 / hR;

  intg_R = matter_obj.alpha_s(hR) * sud_z_QQ_w_M_vac_only(M, h0, hR, loc_d, t_form, E1,matter_obj) *
           (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;

  if ((diff > approx) || (span > error)) {
    intg = sud_val_QQ_w_M_vac_only(M, h0, h1, h, loc_d, E1,matter_obj) +
           sud_val_QQ_w_M_vac_only(M, h0, h, h2, loc_d, E1,matter_obj);
  }

  //	cout << " returning with intg = " << intg << endl;

  return (intg);
}

double sud_z_QQ(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj) {

  //JSINFO<<"sud_z_QQ cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" l_fac "<<l_fac<<" E2 "<<E2;
  double t2, t4, t5, t7, t9, t14, q2, q3, q5, q6, q8, q15, qL, tau, res, z_min;

  z_min = std::sqrt(2) * E_minimum / E2;

  //    if (cg<cg1*z_min) cg = cg1*z_min;

  //    if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );

  if (cg1 < 2.0 * cg) {

    //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
    return (0.0);
  };

  t2 = 1.0 / cg1;
  t4 = 1.0 - cg * t2;
  t5 = t4 * t4;
  t7 = std::pow(cg, 2.0);
  t9 = std::pow(cg1, 2.0);
  t14 = ((t5 * t4) - t7 * cg / t9 / cg1) * t2 / 3.0;

  //	return(t25);

  q2 = 1.0 / cg1;
  q3 = (cg * q2);
  q5 = 1.0 - q3;
  q6 = std::log(std::abs(q5));
  q8 = std::log(q3);
  /*	q10 = std::log(q3);
	q12 = std::log(q5);*/
  q15 = (-1.0 + (2.0 * q3) + q6 - q8) * q2;

  if (q15 < 0.0) {
    cerr << "ERROR: medium contribution negative in sud_z_QQ: q15 = " << q15
         << endl;
    cerr << "cg, cg1 = " << cg << "  " << cg1 << endl;
    cerr << " t14 = " << t14 << endl;
    throw std::runtime_error("ERROR: medium contribution negative in sud_z_QQ");
  }

  tau = l_fac;

  if ((matter_obj.length- loc_e) < tau)
    tau = (matter_obj.length- loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = matter_obj.qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  res = t14 + 2.0 * qL * q15 / cg1;

  return (res);
}

double sud_z_QQ_w_M_vac_only(double M, double cg, double cg1,double loc_e, double l_fac, double E2, Matter& matter_obj) {

  //JSINFO<<"sud_z_QQ_w_M_vac_only cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" l_fac "<<l_fac<<" E2 "<<E2;
  double q2, q3, q5, q6, q8, q15, qL, tau, res, z_min;

  z_min = std::sqrt(2) * E_minimum / E2;

  //    if (cg<cg1*z_min) cg = cg1*z_min;

  //    if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );

  if (cg1 < 2.0 * (cg + M * M)) {

    //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
    return (0.0);
  };

  double t1 = M * M;
  double t2 = t1 + cg;
  double t3 = 1.0 / cg1;
  double t5 = -t2 * t3 + 1.0;
  double t6 = t5 * t5;
  double t8 = t2 * t2;
  double t10 = pow(cg1, 2.0);
  double t15 = 2.0 / 3.0 * (t6 * t5 - t8 * t2 / t10 / cg1) * t3;

  q2 = 1.0 / cg1;
  q3 = (cg * q2);
  q5 = 1.0 - q3;
  q6 = std::log(std::abs(q5));
  q8 = std::log(q3);
  /*	q10 = std::log(q3);
	q12 = std::log(q5);*/
  q15 = (-1.0 + (2.0 * q3) + q6 - q8) * q2;

  if (q15 < 0.0) {
    cerr << "ERROR: medium contribution negative in sud_z_QQ: q15 = " << q15
         << endl;
    cerr << "cg, cg1 = " << cg << "  " << cg1 << endl;
    cerr << " t15 = " << t15 << endl;
    throw std::runtime_error("ERROR: medium contribution negative in sud_z_QQ");
  }

  tau = l_fac;

  if ((matter_obj.length- loc_e) < tau)
    tau = (matter_obj.length- loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = matter_obj.qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    qL = matter_obj.qhat * tau;
  }

  res = t15 + 2.0 * qL * q15 / cg1;

  return (res);
}

double sudakov_Pqp(double g0, double g1, double loc_c, double E, Matter& matter_obj) {
  double sud, g;
  int blurb;
  //JSINFO<<"sudakov_Pqp g0 "<<g0<<" g1 "<<g1<<" loc_c "<<loc_c<<" E "<<E;

  sud = 1.0;

  if (g1 < 2.0 * g0) {
    JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, "
              "returning 1 ";
    JSWARN << " in sudakov_Pquark Photon, g0, g1 = " << g0 << "  " << g1;
    return (sud);
  }
  g = 2.0 * g0;

  double logsud = sud_val_QP(g0, g, g1, loc_c, E,matter_obj);

  sud = exp((-1.0 / 2.0 / pi) * logsud);

  return (sud);
}

double sud_val_QP(double h0, double h1, double h2, double loc_d,double E1, Matter& matter_obj) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  int blurb;
  //JSINFO<<"sud_val_QP h0 "<<h0<<" h1 "<<h1<<" h2 "<<h2<<" loc_d "<<loc_d<<" E1 "<<E1;

  double alphaEM = 1.0 / 137.0;

  val = 0.0;

  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  t_form = 2.0 * E1 / h;

  val = alphaEM * sud_z_QP(h0, h, loc_d, t_form, E1,matter_obj);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  t_form = 2.0 * E1 / hL;

  intg_L = alphaEM * sud_z_QP(h0, hL, loc_d, t_form, E1,matter_obj) * (h - h1);

  hR = (h + h2) / 2.0;

  t_form = 2.0 * E1 / hR;

  intg_R = alphaEM * sud_z_QP(h0, hR, loc_d, t_form, E1,matter_obj) * (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  if ((diff > approx) || (span > error)) {
    intg = sud_val_QP(h0, h1, h, loc_d, E1,matter_obj) + sud_val_QP(h0, h, h2, loc_d, E1,matter_obj);
  }

  return (intg);
}

double sud_z_QP(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj) {

  //JSINFO<<"sud_z_QP cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" l_fac "<<l_fac<<" E2 "<<E2;
  double t2, t6, t10, t11, t17, q2, q3, q4, q5, q6, q10, q14, qL, tau, res,
      z_min;
  int blurb;

  z_min = std::sqrt(2) * E_minimum / E2;

  if (cg1 < 2.0 * cg) {
    return (0.0);
  };

  t2 = std::pow(cg1, 2);
  t6 = std::log(cg);
  t10 = std::abs(cg - cg1);
  t11 = std::log(t10);
  t17 = -1.0 / t2 * (3.0 * cg1 - 6.0 * cg + 4.0 * t6 * cg1 - 4.0 * t11 * cg1) /
        2.0;

  //    return(t17);

  q14 = 0.0;

  tau = l_fac;

  if ((matter_obj.length- loc_e) < tau)
    tau = (matter_obj.length- loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = matter_obj.qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau) * Cf /
           Ca; //for photon production, only the quark scatters
    qL = matter_obj.qhat * tau;
  }

  //JSINFO << BOLDRED << " matter_obj.qhat L = " << qL << " location = " << loc_e << " tau = " << tau << " matter_obj.length= " << length;

  res = t17 + 2.0 * qL * q14 / cg1;

  //   cout << " t0 , t , res = " << cg << "  "  << cg1 << "   " << res << endl ;

  if (q14 < 0.0) {
    cerr << "ERROR: medium contribution negative in sud_z_QG : q14 = " << q14
         << endl;
    throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG");
  }

  return (res);
}

double sudakov_Pqg(double g0, double g1, double loc_c, double E, Matter& matter_obj) {
  double sud, g;
  int blurb;
  //JSINFO<<"sudakov_Pqg g0 "<<g0<<" g1 "<<g1<<" loc_c "<<loc_c<<" E "<<E;
  sud = 1.0;

  if (g1 < 2.0 * g0) {
    JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, "
              "returning 1 ";
    JSWARN << " in sudakov_Pquark gluon, g0, g1 = " << g0 << "  " << g1;
    return (sud);
  }
  g = 2.0 * g0;

  sud = exp(-1.0 * (Cf / 2.0 / pi) * sud_val_QG(g0, g, g1, loc_c, E,matter_obj));

  return (sud);
}

double sudakov_Pqg_w_M(double M, double g0, double g1, double loc_c,double E, Matter& matter_obj) {
  double sud, g;
  int blurb;
  //JSINFO<<"sudakov_Pqg_w_M g0 "<<g0<<" g1 "<<g1<<" loc_c "<<loc_c<<" E "<<E;
  sud = 1.0;

  if (g1 < g0 * (1.0 + std::sqrt(1.0 + 2.0 * M * M / g0))) {
    JSWARN << " warning: Not enough separation between upper and lower limits "
              "of Sudakov to have resolvable radiation ";
    JSWARN << " in sudakov_Pquark gluon, g0*( 1.0 + std::sqrt( 1.0 + "
              "2.0*M*M/g0 ) ) = "
           << g0 * (1.0 + std::sqrt(1.0 + 2.0 * M * M / g0)) << " g1 =  " << g1;
    JSWARN << " M = " << M;

    return (sud);
  }
  g = g0 * (1.0 + std::sqrt(1.0 + 2.0 * M * M / g0));

  sud = exp(-1.0 * (Cf / 2.0 / pi) * sud_val_QG_w_M(M, g0, g, g1, loc_c, E,matter_obj));

  return (sud);
}

double sud_val_QG(double h0, double h1, double h2, double loc_d, double E1, Matter& matter_obj) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  int blurb;
  //JSINFO<<"sud_val_QG h0 "<<h0<<" h1 "<<h1<<" h2 "<<h2<<" loc_d "<<loc_d<<" E1 "<<E1;
  val = 0.0;

  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  t_form = 2.0 * E1 / h;

  val = matter_obj.alpha_s(h) * sud_z_QG(h0, h, loc_d, t_form, E1,matter_obj);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  t_form = 2.0 * E1 / hL;

  intg_L = matter_obj.alpha_s(hL) * sud_z_QG(h0, hL, loc_d, t_form, E1,matter_obj) * (h - h1);

  hR = (h + h2) / 2.0;

  t_form = 2.0 * E1 / hR;

  intg_R = matter_obj.alpha_s(hR) * sud_z_QG(h0, hR, loc_d, t_form, E1,matter_obj) * (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  if ((diff > approx) || (span > error)) {
    intg = sud_val_QG(h0, h1, h, loc_d, E1,matter_obj) + sud_val_QG(h0, h, h2, loc_d, E1,matter_obj);
  }

  return (intg);
}

double sud_val_QG_w_M(double M, double h0, double h1, double h2,double loc_d, double E1, Matter& matter_obj) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  int blurb;
  //JSINFO<<"sud_val_QG_w_M h0 "<<h0<<" h1 "<<h1<<" h2 "<<h2<<" loc_d "<<loc_d<<" E1 "<<E1;  

  val = 0.0;

  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  t_form = 2.0 * E1 / h;

  val = matter_obj.alpha_s(h) * sud_z_QG_w_M(M, h0, h, loc_d, t_form, E1,matter_obj);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  t_form = 2.0 * E1 / hL;

  intg_L = matter_obj.alpha_s(hL) * sud_z_QG_w_M(M, h0, hL, loc_d, t_form, E1,matter_obj) * (h - h1);

  hR = (h + h2) / 2.0;

  t_form = 2.0 * E1 / hR;

  intg_R = matter_obj.alpha_s(hR) * sud_z_QG_w_M(M, h0, hR, loc_d, t_form, E1,matter_obj) * (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;

  if ((diff > approx) || (span > error)) {
    intg = sud_val_QG_w_M(M, h0, h1, h, loc_d, E1,matter_obj) +
           sud_val_QG_w_M(M, h0, h, h2, loc_d, E1,matter_obj);
  }

  //    cout << " returning with intg = " << intg << endl;

  return (intg);
}

double sud_z_QG(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj) {
  //JSINFO<<"sud_z_QG cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" l_fac "<<l_fac<<" E2 "<<E2;  

  double t2, t6, t10, t11, t17, q2, q3, q4, q5, q6, q10, q14, qL, tau, res,
      z_min;
  int blurb;

  z_min = std::sqrt(2) * E_minimum / E2;

  if (cg1 < 2.0 * cg) {
    return (0.0);
  };

  t2 = std::pow(cg1, 2);
  t6 = std::log(cg);
  t10 = std::abs(cg - cg1);
  t11 = std::log(t10);
  t17 = -1.0 / t2 * (3.0 * cg1 - 6.0 * cg + 4.0 * t6 * cg1 - 4.0 * t11 * cg1) /
        2.0;

  //	return(t17);

  q2 = 1.0 / cg1;
  q3 = cg * q2;
  q4 = q3 - 1.0;
  q5 = std::abs(q4);
  q6 = std::log(q5);
  q10 = std::log(q3);
  q14 = (q6 + 2.0 / cg * cg1 - q10 + 2.0 / q4) * q2;

  tau = l_fac;

  if ((matter_obj.length- loc_e) < tau)
    tau = (matter_obj.length- loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = matter_obj.qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    if (matter_obj.qhat * sqrt(2) > 0.6) {
      // JSINFO << BOLDYELLOW << " matter_obj.length= " << matter_obj.length<< " loc = " << loc_e << " tau = " << tau ;
      //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
      // JSINFO << BOLDYELLOW << " mean matter_obj.qhat for sudakov in GeV^2/fm = " << matter_obj.qhat*5*sqrt(2) ;
    }
    qL = matter_obj.qhat * tau;
  }

  //JSINFO << BOLDRED << " matter_obj.qhat L = " << qL << " location = " << loc_e << " tau = " << tau << " matter_obj.length= " << length;

  res = t17 + 2.0 * qL * q14 / cg1;

  //   cout << " t0 , t , res = " << cg << "  "  << cg1 << "   " << res << endl ;

  if (q14 < 0.0) {
    cerr << "ERROR: medium contribution negative in sud_z_QG : q14 = " << q14
         << endl;
    throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG");
  }

  return (res);
}

double sud_z_QG_w_M(double M, double cg, double cg1, double loc_e,double l_fac, double E2, Matter& matter_obj){
  //JSINFO<<"sud_z_QG_w_M cg "<<cg<<" cg1 "<<cg1<<" loc_e "<<loc_e<<" l_fac "<<l_fac<<" E2 "<<E2;  
  double qL, tau, res, z_min;
  int blurb;

  z_min = std::sqrt(2) * E_minimum / E2;

  if (cg1 < 2.0 * cg + M * M / (1.0 + M * M / cg1)) {

    //JSINFO << MAGENTA << " returning with cg, cg1 = " << cg << "   " << cg1
    //       << "    " << E_minimum << "  " << E2;
    return (M * M);
  };

  double t1 = 1.0 / cg1;
  double t2 = t1 * cg;
  double t4 = std::pow(1.0 - t2, 2.0);
  double t7 = std::log(t2);
  double t9 = M * M;
  double t10 = t1 * t9;
  double t13 = 1.0 / (t10 + 1.0) * t10;
  double t15 = std::pow(t2 + t13, 2.0);
  double t18 = std::log(1.0 - t2 - t13);
  double t21 = t1 * (-t4 / 2.0 - 1.0 + 2.0 * t2 - 2.0 * t7 + t15 / 2.0 + t13 +
                     2.0 * t18);

  double q1 = M * M;
  double q2 = 1.0 / cg1;
  double q3 = q2 * q1;
  double q5 = 4.0 * q3 + 1.0;
  double q6 = q2 * cg;
  double q7 = std::log(q6);
  double q9 = q1 * q1;
  double q10 = std::pow(cg1, 2.0);
  double q12 = 1.0 / q10 * q9;
  double q14 = q12 - 2.0 * q3 + 1.0 / 2.0;
  double q15 = 1.0 - q6;
  double q16 = std::log(q15);
  double q18 = q3 + 1.0;
  double q28 = q15 * q15;
  double q33 = 1.0 / q18 * q3;
  double q34 = 1.0 - q6 - q33;
  double q35 = std::log(q34);
  double q37 = q6 + q33;
  double q38 = std::log(q37);
  double q48 = q37 * q37;
  double q52 = q7 * q5 + q16 * q14 + q15 * q18 / 2.0 + 2.0 / cg * cg1 +
               3.0 / 2.0 * q2 / q15 * q1 - 1.0 / q28 * q12 / 2.0 - q35 * q5 -
               q38 * q14 - q37 * q18 / 2.0 - 2.0 / q34 -
               3.0 / 2.0 * q2 / q37 * q1 + 1.0 / q48 * q12 / 2.0;
  double q53 = q2 * q52;

  tau = l_fac;

  if ((matter_obj.length- loc_e) < tau)
    tau = (matter_obj.length- loc_e);

  if (loc_e > matter_obj.length)
    tau = 0.0;

  // SC
  //qL = matter_obj.qhat*0.6*tau*profile(loc_e + tau) ;
  if (tau < rounding_error) {
    qL = 0.0;
  } else {
    matter_obj.qhat = matter_obj.fncAvrQhat(loc_e, tau);
    // if (matter_obj.qhat*sqrt(2)>0.6)
    // {
    //   JSINFO << MAGENTA << " Big q-hat warning ! ";
    //   JSINFO << BOLDYELLOW << " matter_obj.length= " << matter_obj.length<< " loc = " << loc_e << " tau = " << tau ;
    //   JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
    //   JSINFO << BOLDYELLOW << " mean matter_obj.qhat for sudakov in GeV^2/fm = " << matter_obj.qhat*5*sqrt(2) ;
    // }
    qL = matter_obj.qhat * 2.0 * tau;
  }

  double e1 = M * M;
  double e2 = 1.0 / cg1;
  double e3 = e2 * e1;
  double e4 = e2 * cg;
  double e5 = std::log(e4);
  double e8 = std::log(1.0 - e4);
  double e13 = 1.0 / (e3 + 1.0) * e3;
  double e15 = std::log(1.0 - e4 - e13);
  double e18 = std::log(e4 + e13);
  double e22 = e2 * (-(2.0 * e5 - e8 + 1.0 - e4) * e3 +
                     (2.0 * e15 - e18 + e4 + e13) * e3);

  double eL;

  if (tau < rounding_error) {
    eL = 0.0;
  } else {
    matter_obj.ehat = 0.0; //fncAvrEhat(loc_e,tau);
    if (matter_obj.ehat * sqrt(2) > 0.6) {
      // JSINFO << BOLDYELLOW << " matter_obj.length= " << matter_obj.length<< " loc = " << loc_e << " tau = " << tau ;
      //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
      // JSINFO << BOLDYELLOW << " mean matter_obj.qhat for sudakov in GeV^2/fm = " << matter_obj.qhat*5*sqrt(2) ;
    }
    eL = matter_obj.ehat * 4.0;
  }

  double f1 = M * M;
  double f2 = 1.0 / cg1;
  double f3 = f2 * f1;
  double f4 = f2 * cg;
  double f5 = std::log(f4);
  double f8 = f1 * f1;
  double f9 = std::pow(cg1, 2.0);
  double f11 = 1.0 / f9 * f8;
  double f14 = 13.0 / 4.0 * f11 - 15.0 / 4.0 * f3 + 1.0 / 2.0;
  double f15 = 1.0 - f4;
  double f16 = std::log(f15);
  double f24 = f15 * f15;
  double f32 = 1.0 / (f3 + 1.0) * f3;
  double f33 = 1.0 - f4 - f32;
  double f34 = std::log(f33);
  double f37 = f4 + f32;
  double f38 = std::log(f37);
  double f45 = f37 * f37;
  double f52 = f2 * ((15.0 / 2.0 * f5 * f3 + f16 * f14 + 1.0 / cg * cg1 +
                      15.0 / 4.0 * f2 / f15 * f1 - 13.0 / 8.0 / f24 * f11) *
                         f3 -
                     (15.0 / 2.0 * f34 * f3 + f38 * f14 + 1.0 / f33 +
                      15.0 / 4.0 * f2 / f37 * f1 - 13.0 / 8.0 / f45 * f11) *
                         f3);
  double e2L;

  if (tau < rounding_error) {
    e2L = 0.0;
  } else {
    matter_obj.e2hat = matter_obj.qhat / 2.0; //fncAvrE2hat(loc_e,tau);
    if (matter_obj.e2hat * sqrt(2) > 0.6) {
      // JSINFO << BOLDYELLOW << " matter_obj.length= " << matter_obj.length<< " loc = " << loc_e << " tau = " << tau ;
      //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
      // JSINFO << BOLDYELLOW << " mean matter_obj.qhat for sudakov in GeV^2/fm = " << matter_obj.qhat*5*sqrt(2) ;
    }
    e2L = matter_obj.e2hat * 8.0 / (tau * cg1);
  }

  //    JSINFO << BOLDRED << " matter_obj.qhat L = " << qL << " location = " << loc_e << " tau = " << tau << " matter_obj.length= " << length;

  res = t21 + qL * q53 / cg1;

  //+ eL*e22/cg1 +e2L*f52/cg1;
  // Uncomment only if you have an eL larger than 2 times e2L for charm, and derive expression for bottom.
  // MC simulation is not valid for all choices of e-hat and e2-hat.

  if (res < 0.0) {
    cerr << "ERROR: medium contribution negative in sud_z_QG : res = " << res
         << endl;

    throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG");
  }

  return (res);
}

}