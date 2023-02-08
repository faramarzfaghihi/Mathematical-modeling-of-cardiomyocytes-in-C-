

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <array>
#include "matplotlibcpp.h"
using namespace std;
#include <vector>
namespace plt = matplotlibcpp;


//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

double A0_bck = 1.0278; // dimensionless (in A0_bc)
double Ampl_gain = 1; // dimensionless (in Ampl_gai)
double Cext=1; // picoF (in Cex)
double E_l=1; // millivolt (in E_)
double I=1; // dimensionless (in )
double R_seal=1; // ohm (in R_sea)
double Scale_bck=1; // dimensionless (in Scale_bc)
double Scaling=1; // dimensionless (in Scalin)
double k_bck = 0.098599999999999993; // one_over_millivolt (in k_bc)
double leak_comp_perc=1; // dimensionless (in leak_comp_per)
double Vol_leak = 0.00036000000000000002; // per_millisecond (in Ileak_Iup_Ixfer)
double Vol_rel = 0.30599999999999999; // per_millisecond (in Irel)
double Vmax_up = 0.0063749999999999996; // millimolar_per_millisecond (in Ileak_Iup_Ixfer)
double Cm = 0.115;      

double conc_clamp=1; // dimensionless (in Na)
double Ca_o=2; // millimolar (in Environment)
double K_o = 5.4000000000000004; // millimolar (in Environment)
double Na_o=140; // millimolar (in Environment)
double g_CaL = 2.0000000000000002e-5; // litre_per_farad_millisecond (in ICaL)
double g_Na=11; // microS_per_nanoF (in INa)
double perc_reduced_inact_for_IpNa=0; // dimensionless (in INa)
double shift_INa_inact=0; // millivolt (in INa)
double g_K1_0 = 0.68210000000000004; // microS_per_nanoF (in IK1)
double g_Kr_0 = 0.024; // microS_per_nanoF (in IKr)
double g_Ks = 0.039199999999999999; // microS_per_nanoF (in IKs)
double K_NaCa = 200; // nanoA_per_nanoF (in INaCa)

double stim_amplitude = -12; // nanoA_per_nanoF (in cell)
double stim_duration = 3; // millisecond (in cell)
double stim_start = 100; // millisecond (in cell)
double stim_period = 7000; // millisecond (in cell)
double g_to = 0.20000000000000001; // microS_per_nanoF (in Ito)  
    
double V_sr = 0.0010939999999999999; // nanolitre (in Ca)
double V_ss = 5.4679999999999998e-5; // nanolitre (in Ca)
double Buf_c = 0.20000000000000001; // millimolar (in Ca_buffer)
double Buf_sr = 10; // millimolar (in Ca_buffer)
double Buf_ss = 0.40000000000000002; // millimolar (in Ca_buffer)
double K_buf_c = 0.001; // millimolar (in Ca_buffer)
double K_buf_sr = 0.29999999999999999; // millimolar (in Ca_buffer)
double K_buf_ss = 0.00025000000000000001; // millimolar (in Ca_buffer)
double F = 96485.341499999995; // coulomb_per_mole (in Environment)
double R = 8314.4719999999998; // millijoule_per_mole_kelvin (in Environment)
double T = 310; // kelvin (in Environment)
double g_bca = 0.00047360000000000002; // microS_per_nanoF (in ICab)
double K_sat = 0.10000000000000001; // dimensionless (in INaCa)
double Km_Ca = 1.3799999999999999; // millimolar (in INaCa)
double Km_Nai = 87.5; // millimolar (in INaCa)
double alpha = 2.5; // dimensionless (in INaCa)
double gama = 0.34999999999999998; // dimensionless (in INaCa)
double K_mNa = 40; // millimolar (in INaK)
double K_mk = 1; // millimolar (in INaK)
double P_NaK = 1.2969999999999999; // nanoA_per_nanoF (in INaK)
double g_bna = 0.00029; // microS_per_nanoF (in INab)
double K_up = 0.00025000000000000001; // millimolar (in Ileak_Iup_Ixfer)

double Vol_xfer = 0.0038; // per_millisecond (in Ileak_Iup_Ixfer)
double K_pCa = 0.00050000000000000001; // millimolar (in IpCa)
double g_pCa = 0.061899999999999997;// nanoA_per_nanoF (in IpCa)
double g_pK = 0.0097300000000000008; // microS_per_nanoF (in IpK)
double EC = 1.5; // millimolar (in Irel)
double k1_prime = 0.14999999999999999; // per_millimolar2_per_millisecond (in Irel)
double k2_prime = 0.044999999999999998; // per_millimolar_per_millisecond (in Irel)
double k3 = 0.059999999999999998; // per_millisecond (in Irel)
double k4 = 0.0050000000000000001; // per_millisecond (in Irel)
double max_sr = 2.5; // dimensionless (in Irel)
double min_sr = 1; // dimensionless (in Irel)
double Vol_c = 0.016403999999999998; // nanolitre (in cell)
double stim_end = 100000;// millisecond (in cell)
double d_inf_shift = 5; // millivolt (in iCaL_d_gate)
double Mg_Buf = 0.0356; // millimolar (in iK1_rectification)
double SPM = 0.0014613; // millimolar (in iK1_rectification)
double fac = 1.0648; // dimensionless (in iK1_rectification)
double phi = 0.88380000000000003; // dimensionless (in iK1_rectification)
double T_Base = 310; // kelvin (in iKr_Markov)
double kBinding = 0.0050000000000000001; // per_millimolar_per_millisecond (in iKr_Markov_Sotalol_block)
double kDiss = 0.00125; // per_millisecond (in iKr_Markov_Sotalol_block)
double P_kna = 0.029999999999999999; // dimensionless (in reversal_potentials)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------
int main()
{

double dCr1;
double dCr2;
double dCr3;
double dIr5;
double dBCr1;
double dBCr2;
double dBCr3;
double dBIr5;
double dOr4;
double dBOr4;
double dXs;
double dCa_i;
double dNa_i;
double dK_i;

double dm;
double dh;
double djj;
double dr;
double ds;
double dV;

double dCa_ss;
double Ca_i_bufc; // dimensionless (in Ca_buffer)
double Ca_sr_bufsr; // dimensionless (in Ca_buffer)
double Ca_ss_bufss; // dimensionless (in Ca_buffer)
double i_leak; // millimolar_per_millisecond (in Ileak_Iup_Ixfer)
double i_up; // millimolar_per_millisecond (in Ileak_Iup_Ixfer)
double i_xfer; // millimolar_per_millisecond (in Ileak_Iup_Ixfer)
double i_p_Ca; // nanoA_per_nanoF (in IpCa)
double kcasr; // dimensionless (in Irel)
double k1; // per_millimolar2_per_millisecond (in Irel)
double O; // dimensionless (in Irel)
double i_rel; // millimolar_per_millisecond (in Irel)
double k2; // per_millimolar_per_millisecond (in Irel)
double i_NaCa; // nanoA_per_nanoF (in INaCa)
double i_NaK; // nanoA_per_nanoF (in INaK)
double i_Stim; // nanoA_per_nanoF (in cell)
double alpha_d; // dimensionless (in iCaL_d_gate)
double beta_d; // dimensionless (in iCaL_d_gate)
double d_inf; // dimensionless (in iCaL_d_gate)
double gamma_d; // millisecond (in iCaL_d_gate)
double tau_d; // millisecond (in iCaL_d_gate)
double f2_inf; // dimensionless (in iCaL_f2_gate)
double tau_f2; // millisecond (in iCaL_f2_gate)
double fCass_inf; // dimensionless (in iCaL_fCass_gate)
double tau_fCass; // millisecond (in iCaL_fCass_gate)
double i_CaL; // nanoA_per_nanoF (in ICaL)
double f_inf; // dimensionless (in iCaL_f_gate)
double tau_f; // millisecond (in iCaL_f_gate)
double alpha_xr1; // per_millisecond (in iKr_Markov)
double alpha_xr2; // per_millisecond (in iKr_Markov)
double alpha_xr3; // per_millisecond (in iKr_Markov)
double alpha_xr4; // per_millisecond (in iKr_Markov)
double beta_xr1; // per_millisecond (in iKr_Markov)
double beta_xr2; // per_millisecond (in iKr_Markov)
double beta_xr3; // per_millisecond (in iKr_Markov)
double beta_xr4; // per_millisecond (in iKr_Markov)
double Sotalol_mM; // millimolar (in iKr_Markov_Sotalol_block)
double OtoB; // per_millisecond (in iKr_Markov_Sotalol_block)
double BtoO; // per_millisecond (in iKr_Markov_Sotalol_block)
double alpha_xs; // dimensionless (in iKs_Xs_gate)
double beta_xs; // dimensionless (in iKs_Xs_gate)
double tau_xs; // millisecond (in iKs_Xs_gate)
double xs_inf; // dimensionless (in iKs_Xs_gate)
double alpha_h; // per_millisecond (in iNa_h_gate)
double beta_h; // per_millisecond (in iNa_h_gate)
double h_inf; // dimensionless (in iNa_h_gate)
double tau_h; // millisecond (in iNa_h_gate)
double alpha_j; // per_millisecond (in iNa_j_gate)
double beta_j; // per_millisecond (in iNa_j_gate)
double j_inf; // dimensionless (in iNa_j_gate)
double tau_j; // millisecond (in iNa_j_gate)
double alpha_m; // dimensionless (in iNa_m_gate)
double beta_m; // dimensionless (in iNa_m_gate)
double m_inf; // dimensionless (in iNa_m_gate)
double tau_m; // millisecond (in iNa_m_gate)
double r_inf; // dimensionless (in ito_r_gate)
double tau_r; // millisecond (in ito_r_gate)
double s_inf; // dimensionless (in ito_s_gate)
double tau_s; // millisecond (in ito_s_gate)
double E_Ca; // millivolt (in reversal_potentials)
double i_b_Ca; // nanoA_per_nanoF (in ICab)
double E_K; // millivolt (in reversal_potentials)
double i_Kr; // nanoA_per_nanoF (in IKr)
double i_p_K; // nanoA_per_nanoF (in IpK)
double i_to; // nanoA_per_nanoF (in Ito)
double KbMg; // millimolar (in iK1_rectification)
double Kd1SPM; // millimolar (in iK1_rectification)
double Kd2SPM; // millimolar (in iK1_rectification)
double KiMg; // millimolar (in iK1_rectification)
double rec2; // dimensionless (in iK1_rectification)
double temp; // dimensionless (in iK1_rectification)
double rec1; // dimensionless (in iK1_rectification)
double xK1_inf; // dimensionless (in iK1_rectification)
double i_K1; // nanoA_per_nanoF (in IK1)
double E_Na; // millivolt (in reversal_potentials)
double i_Na; // nanoA_per_nanoF (in INa)
double i_b_Na; // nanoA_per_nanoF (in INab)
double E_Ks; // millivolt (in reversal_potentials)
double i_Ks; // nanoA_per_nanoF (in IKs)
double i_tot; // nanoA_per_nanoF (in cell)
//------------------------------------------------------------------------------
std::vector<double> t(5000);
std::vector<double> V(t.size());

//------------------------------------------------------------------------------
double dCa_SR;
double dR_prime;
std::vector<double> dd(t.size());
std::vector<double> df2(t.size());
std::vector<double> df(t.size());
double dfCass;

std::vector<double> Ca_i(t.size()); 
std::vector<double> Cr1(t.size()); 
std::vector<double> Cr2(t.size()); 
std::vector<double> Cr3(t.size()); 
std::vector<double> Or4(t.size()); 
std::vector<double> Ir5(t.size()); 
std::vector<double> BCr1(t.size()); 
std::vector<double> BCr2(t.size()); 
std::vector<double> BCr3(t.size()); 
std::vector<double> BOr4(t.size()); 
std::vector<double> BIr5(t.size()); 
std::vector<double> Xs(t.size()); 
std::vector<double> h(t.size()); 
std::vector<double> r(t.size()); 
std::vector<double> m(t.size()); 
std::vector<double> s(t.size()); 
std::vector<double> jj(t.size()); 
std::vector<double> d(t.size()); 
std::vector<double> f(t.size()); 
std::vector<double> f2(t.size()); 
std::vector<double> fCass(t.size()); 
std::vector<double> Ca_SR(t.size()); 
std::vector<double> Ca_ss(t.size()); 
std::vector<double> R_prime(t.size()); 
std::vector<double> Na_i(t.size()); 
std::vector<double> K_i(t.size()); 


//------------------------------------------------------------------------------
    V[0] = -86.45; // 0 V (millivolt) (in cell)
    Ca_i[0] = 0.0001092; // 1 Ca_i (millimolar) (in Ca)
    Cr1[0] = 0.9786; // 2 Cr1 (dimensionless) (in iKr_Markov)
    Cr2[0] = 0.0031; // 3 Cr2 (dimensionless) (in iKr_Markov)
    Cr3[0] = 0.0029; // 4 Cr3 (dimensionless) (in iKr_Markov)
    Or4[0] = 0.014; // 5 Or4 (dimensionless) (in iKr_Markov)
    Ir5[0] = 0.0014; // 6 Ir5 (dimensionless) (in iKr_Markov)
    BCr1[0] = 0.0; // 7 BCr1 (dimensionless) (in iKr_Markov_Sotalol_block)
    BCr2[0] = 0.0; // 8 BCr2 (dimensionless) (in iKr_Markov_Sotalol_block)
    BCr3[0] = 0.0; // 9 BCr3 (dimensionless) (in iKr_Markov_Sotalol_block)
    BOr4[0] = 0.0; // 10 BOr4 (dimensionless) (in iKr_Markov_Sotalol_block)
    BIr5[0] = 0.0; // 11 BIr5 (dimensionless) (in iKr_Markov_Sotalol_block)
    Xs[0] = 0.00303; // 12 Xs (dimensionless) (in iKs_Xs_gate)
    s[0] = 1.0; // 13 s (dimensionless) (in ito_s_gate)
    r[0] = 2.11e-08; // 14 r (dimensionless) (in ito_r_gate)
    m[0] = 0.00132; // 15 m (dimensionless) (in iNa_m_gate)
    h[0] = 0.7768; // 16 h (dimensionless) (in iNa_h_gate)
    jj[0] = 0.7766; // 17 j (dimensionless) (in iNa_j_gate)
    d[0] = 5.06e-06; // 18 d (dimensionless) (in iCaL_d_gate)
    f[0] = 0.9999; // 19 f (dimensionless) (in iCaL_f_gate)
    f2[0] = 0.9995; // 20 f2 (dimensionless) (in iCaL_f2_gate)
    fCass[0] = 1.0; // 21 fCass (dimensionless) (in iCaL_fCass_gate)
    Ca_SR[0] = 2.7656; // 22 Ca_SR (millimolar) (in Ca)
    Ca_ss[0] = 0.0001893; // 23 Ca_ss (millimolar) (in Ca)
    R_prime[0] = 0.9864; // 24 R_prime (dimensionless) (in Irel)
    Na_i[0] = 7.940167; // 25 Na_i (millimolar) (in Na)
    K_i[0] = 141.0167; // 26 K_i (millimolar) (in K)
    //time[0]; // (time} (milliseconds)

//------------------------------------------------------------------------------

//-------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------


  double dt = 0.08;

  for (size_t j = 0; j<t.size(); j++)
   {  

     if ((j >= stim_start) && (j <= stim_end) && (stim_duration >= -stim_start - stim_period * floor((-stim_start + j) / stim_period) + j))
        i_Stim = 8*stim_amplitude;
      else 
        i_Stim = 0;
    //---------------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------------

    Ca_i_bufc = 1 / (1 + Buf_c * K_buf_c / pow((K_buf_c + Ca_i[j]), 2));

    Ca_sr_bufsr = 1 / (1 + Buf_sr * K_buf_sr / pow((K_buf_sr + Ca_SR[j]), 2));

    Ca_ss_bufss = 1 / (1 + Buf_ss * K_buf_ss / pow((K_buf_ss + Ca_ss[j]), 2));

    i_leak = (-Ca_i[j] + Ca_SR[j]) * Vol_leak;

    i_up = Vmax_up / (1 + pow(K_up, 2) / pow(Ca_i[j], 2));

    i_xfer = (-Ca_i[j] + Ca_ss[j]) * Vol_xfer;

    i_p_Ca = g_pCa * Ca_i[j] / (K_pCa + Ca_i[j]);

    kcasr = -(-min_sr + max_sr) / (1 + pow(EC, 2) / pow(Ca_SR[j], 2)) + max_sr;

    k1 = k1_prime / kcasr;

    O = pow(Ca_ss[j], 2) * k1 * R_prime[j] / (pow(Ca_ss[j], 2) * k1 + k3);

    i_rel = (-Ca_ss[j] + Ca_SR[j]) * O * Vol_rel;

    dCa_SR = (-i_leak - i_rel + i_up) * Ca_sr_bufsr;

    k2 = k2_prime * kcasr;

    dR_prime = (1 - R_prime[j]) * k4 - k2 * Ca_ss[j] * R_prime[j];

    i_NaCa = (pow(Na_i[j], 3) * Ca_o * exp(F * gama * V[j] / (R * T)) - pow(Na_o, 3) * alpha * Ca_i[j] * exp((-1 + gama) * F * V[j] / (R * T))) * K_NaCa / ((1 + K_sat * exp((-1 + gama) * F * V[j] / (R * T))) * (pow(Na_o, 3) + pow(Km_Nai, 3)) * (Ca_o + Km_Ca));
    
    i_NaK = K_o * P_NaK * Na_i[j] / ((K_o + K_mk) * (K_mNa + Na_i[j]) * (1 + 0.035299999999999998 * exp(-F * V[j] / (R * T)) + 0.1245 * exp(-0.10000000000000001 * F * V[j] / (R * T))));
    
    

    
    alpha_d = 0.25 + 1.3999999999999999 / (1 + exp(-2.6923076923076925 - 0.076923076923076927 * V[j]));
    
    beta_d = 1.3999999999999999 / (1 + exp(1 + 0.20000000000000001 * V[j]));

    d_inf = 1 / (1 + exp(0.13333333333333333 * d_inf_shift - 0.13333333333333333 * V[j]));

    gamma_d = 1 / (1 + exp(2.5 - 0.050000000000000003 * V[j]));

    tau_d = alpha_d * beta_d + gamma_d;

    dd[j] = (-d[j] + d_inf) / tau_d;

    f2_inf = 0.25 + 0.75 / (1 + exp(5 + 0.14285714285714285 * V[j]));

    
    tau_f2 = 40 / (1 + exp(3 + 0.10000000000000001 * V[j])) + 15.5 / (1 + exp(2.5 - 0.10000000000000001 * V[j])) + 281 * exp(-3.0375000000000001 * pow((1 + 0.037037037037037035 * V[j]), 2));
    df2[j] = (-f2[j] + f2_inf) / tau_f2;

    fCass_inf = 0.59999999999999998 + 0.40000000000000002 / (1 + 399.99999999999994 * pow(Ca_ss[j], 2));
    
    tau_fCass = 2 + 80 / (1 + 399.99999999999994 * pow(Ca_ss[j], 2));
    
    dfCass = (-fCass[j] + fCass_inf) / tau_fCass;

    if ((V[j] >= 14.999998664311967) && (V[j] <= 15.000001335688033))
        i_CaL = 374338.90822798351 * (-14.999998664311967 + V[j]) * (0.019297068299972742 * (-Ca_o + 0.25 * Ca_ss[j] * exp(9.9999999999858739e-8)) * g_CaL * d[j] * f[j] * f2[j] * fCass[j] / (-1 + exp(9.9999999999858739e-8)) + 0.019297068299972742 * (-Ca_o + 0.25 * Ca_ss[j] * exp(-9.9999999999858739e-8)) * g_CaL * d[j] * f[j] * f2[j] * fCass[j] / (-1 + exp(-9.9999999999858739e-8))) - 0.019297068299972742 * (-Ca_o + 0.25 * Ca_ss[j] * exp(-9.9999999999858739e-8)) * g_CaL * d[j] * f[j] * f2[j] * fCass[j] / (-1 + exp(-9.9999999999858739e-8));
    else 
        i_CaL = 14447.286958825251 * (-15 + V[j]) * (-Ca_o + 0.25 * Ca_ss[j] * exp(-1.1230167246823641 + 0.074867781645490947 * V[j])) * g_CaL * d[j] * f[j] * f2[j] * fCass[j] / (-1 + exp(-1.1230167246823641 + 0.074867781645490947 * V[j]));

    
    dCa_ss = (V_sr * i_rel / V_ss - i_xfer * Vol_c / V_ss - 0.5 * i_CaL * Cm / (V_ss * F)) * Ca_ss_bufss;
    
    f_inf = 1 / (1 + exp(2.8571428571428572 + 0.14285714285714285 * V[j]));
    
    tau_f = 5 + 50 / (1 + exp(1.3 - 0.10000000000000001 * V[j])) + 45 / (1 + exp(3 + 0.10000000000000001 * V[j])) + 275.625 * exp(-3.2400000000000002 * pow((1 + 0.037037037037037035 * V[j]), 2));
    
    df[j] = (-f[j] + f_inf) / tau_f;
    
    alpha_xr1 = T * exp(24.335000000000001 + (-25.914000000000001 + 0.0112 * V[j]) * T_Base / T) / T_Base;
    
    alpha_xr2 = T * exp(22.745999999999999 - 25.914000000000001 * T_Base / T) / T_Base;
    
    alpha_xr3 = T * exp(22.097999999999999 + (-25.914000000000001 + 0.036499999999999998 * V[j]) * T_Base / T) / T_Base;
    
    alpha_xr4 = 1.9631681698237122 * pow((1 / K_o), 0.40000000000000002) * T * exp(30.015999999999998 + (-30.888000000000002 + 0.0223 * V[j]) * T_Base / T) / T_Base;
    
    beta_xr1 = T * exp(13.688000000000001 + (-15.707000000000001 - 0.060299999999999999 * V[j]) * T_Base / T) / T_Base;
    
    dCr1 = Cr2[j] * beta_xr1 - Cr1[j] * alpha_xr1;
    
    beta_xr2 = T * exp(13.193 - 15.707000000000001 * T_Base / T) / T_Base;
    
    dCr2 = Cr1[j] * alpha_xr1 + Cr3[j] * beta_xr2 - (alpha_xr2 + beta_xr1) * Cr2[j];
    
    beta_xr3 = T * exp(7.3129999999999997 + (-15.707000000000001 - 0.039899999999999998 * V[j]) * T_Base / T) / T_Base;
    
    dCr3 = Cr2[j] * alpha_xr2 + Or4[j] * beta_xr3 - (alpha_xr3 + beta_xr2) * Cr3[j];
    
    beta_xr4 = T * exp(30.061 + (-33.243000000000002 - 0.031199999999999999 * V[j]) * T_Base / T) / T_Base;
    
    dIr5 = Or4[j] * alpha_xr4 - Ir5[j] * beta_xr4;
    
    dBCr1 = BCr2[j] * beta_xr1 - BCr1[j] * alpha_xr1;
    
    dBCr2 = BCr1[j] * alpha_xr1 + BCr3[j] * beta_xr2 - (alpha_xr2 + beta_xr1) * BCr2[j];
    
    dBCr3 = BOr4[j] * beta_xr3 + BCr2[j] * alpha_xr2 - (alpha_xr3 + beta_xr2) * BCr3[j];
    
    dBIr5 = BOr4[j] * alpha_xr4 - BIr5[j] * beta_xr4;
    
    Sotalol_mM = 0;
    
    OtoB = Or4[j] * Sotalol_mM * kBinding;
    
    BtoO = BOr4[j] * kDiss;
    
    dOr4 = -OtoB + Cr3[j] * alpha_xr3 + Ir5[j] * beta_xr4 - (alpha_xr4 + beta_xr3) * Or4[j] + BtoO;
    
    dBOr4 = -BtoO + BIr5[j] * beta_xr4 + BCr3[j] * alpha_xr3 - (alpha_xr4 + beta_xr3) * BOr4[j] + OtoB;
    
    alpha_xs = 1400 / sqrt(1 + exp(0.83333333333333337 - 0.16666666666666666 * V[j]));
    
    beta_xs = 1 / (1 + exp(-2.3333333333333335 + 0.066666666666666666 * V[j]));
    
    tau_xs = 80 + alpha_xs * beta_xs;
    
    xs_inf = 1 / (1 + exp(-0.35714285714285715 - 0.071428571428571425 * V[j]));
    
    dXs = (-Xs[j] + xs_inf) / tau_xs;
    
    if (V[j] < -40 + shift_INa_inact)
        alpha_h = 0.057000000000000002 * exp(-11.764705882352942 + 0.14705882352941177 * shift_INa_inact - 0.14705882352941177 * V[j]);
    else 
        alpha_h = 0;

    
    if (V[j] < -40 + shift_INa_inact)
        beta_h = 310000 * exp(0.34849999999999998 * V[j] - 0.34849999999999998 * shift_INa_inact) + 2.7000000000000002 * exp(0.079000000000000001 * V[j] - 0.079000000000000001 * shift_INa_inact);
    else 
        beta_h = 5.9230769230769234 / (1 + exp(-0.96036036036036043 + 0.0900900900900901 * shift_INa_inact - 0.0900900900900901 * V[j]));

    
    h_inf = 0.01 * perc_reduced_inact_for_IpNa + (1 - 0.01 * perc_reduced_inact_for_IpNa) / pow((1 + exp(9.6298788694481825 + 0.13458950201884254 * V[j] - 0.13458950201884254 * shift_INa_inact)), 2);
    
    tau_h = 1 / (alpha_h + beta_h);
    
    dh = (-h[j] + h_inf) / tau_h;
    
    if (V[j] < -40 + shift_INa_inact)
        alpha_j = (37.780000000000001 + V[j]) * (-25428 * exp(0.24440000000000001 * V[j] - 0.24440000000000001 * shift_INa_inact) - 6.9480000000000002e-6 * exp(0.043909999999999998 * shift_INa_inact - 0.043909999999999998 * V[j])) / (1 + exp(24.640530000000002 + 0.311 * V[j] - 0.311 * shift_INa_inact));
    else 
        alpha_j = 0;

    
    if (V[j] < -40 + shift_INa_inact)
        beta_j = 0.024240000000000001 * exp(0.01052 * shift_INa_inact - 0.01052 * V[j]) / (1 + exp(-5.5312920000000005 + 0.13780000000000001 * shift_INa_inact - 0.13780000000000001 * V[j]));
    else 
        beta_j = 0.59999999999999998 * exp(0.057000000000000002 * V[j] - 0.057000000000000002 * shift_INa_inact) / (1 + exp(-3.2000000000000002 + 0.10000000000000001 * shift_INa_inact - 0.10000000000000001 * V[j]));

    
    j_inf = 0.01 * perc_reduced_inact_for_IpNa + (1 - 0.01 * perc_reduced_inact_for_IpNa) / pow((1 + exp(9.6298788694481825 + 0.13458950201884254 * V[j] - 0.13458950201884254 * shift_INa_inact)), 2);
    
    tau_j = 1 / (alpha_j + beta_j);
    
    djj = (-jj[j] + j_inf) / tau_j;
    
    alpha_m = 1 / (1 + exp(-12 - 0.20000000000000001 * V[j]));
    
    beta_m = 0.10000000000000001 / (1 + exp(7 + 0.20000000000000001 * V[j])) + 0.10000000000000001 / (1 + exp(-0.25 + 0.0050000000000000001 * V[j]));
    
    m_inf = 1 / pow((1 + exp(-6.2967884828349945 - 0.11074197120708749 * V[j])), 2);
    
    tau_m = alpha_m * beta_m;
    
    dm = (-m[j] + m_inf) / tau_m;
    
    r_inf = 1 / (1 + exp(3.3333333333333335 - 0.16666666666666666 * V[j]));
    
    tau_r = 0.80000000000000004 + 9.5 * exp(-0.88888888888888884 * pow((1 + 0.025000000000000001 * V[j]), 2));
    
    dr = (-r[j] + r_inf) / tau_r;
    
    s_inf = 1 / (1 + exp(4 + 0.20000000000000001 * V[j]));
    
    tau_s = 3 + 5 / (1 + exp(-4 + 0.20000000000000001 * V[j])) + 85 * exp(-6.328125 * pow((1 + 0.022222222222222223 * V[j]), 2));
    
    ds = (-s[j] + s_inf) / tau_s;
    
    E_Ca = 0.5 * R * T * log(Ca_o / Ca_i[j]) / F;

    i_b_Ca = (-E_Ca + V[j]) * g_bca;

    
    dCa_i = ((-i_up + i_leak) * V_sr / Vol_c - 0.5 * (-2 * i_NaCa + i_b_Ca + i_p_Ca) * Cm / (F * Vol_c) + i_xfer) * Ca_i_bufc;
    
    E_K = R * T * log(K_o / K_i[j]) / F;
    
    i_Kr = 0.43033148291193518 * sqrt(K_o) * (-7.8571428571428568 + 0.028571428571428571 * T) * (-E_K + V[j]) * g_Kr_0 * Or4[j];
    
    i_p_K = (-E_K + V[j]) * g_pK / (1 + exp(4.1806020066889626 - 0.16722408026755853 * V[j]));
    
    i_to = (-E_K + V[j]) * g_to * s[j] * r[j];
    
    KbMg = 0.45000000000000001 * exp(-0.050000000000000003 * V[j] + 0.050000000000000003 * fac * E_K);
    
    Kd1SPM = 0.00069999999999999999 * exp(-0.20833333333333334 * V[j] - 1.6666666666666667 * Mg_Buf + 0.20833333333333334 * fac * E_K);
    
    Kd2SPM = 0.040000000000000001 * exp(-0.10989010989010989 * V[j] + 0.10989010989010989 * fac * E_K);
    
    KiMg = 2.7999999999999998 * exp(-0.0055555555555555558 * V[j] + 0.0055555555555555558 * fac * E_K);
    
    rec2 = 1 / (1 + SPM / Kd2SPM);
    
    temp = 1 + Mg_Buf / KbMg;
    
    rec1 = pow(temp, 2) / (pow(temp, 3) + SPM / Kd1SPM + Mg_Buf / KiMg);
    
    xK1_inf = (1 - phi) * rec2 + phi * rec1;
    
    i_K1 = 0.43033148291193518 * sqrt(K_o) * (-7.8571428571428568 + 0.028571428571428571 * T) * (-E_K + V[j]) * g_K1_0 * xK1_inf;
    
    E_Na = R * T * log(Na_o / Na_i[j]) / F;
    
    i_Na = pow(m[j], 3) * (-E_Na + V[j]) * g_Na * h[j] * jj[j];
    
    i_b_Na = (-E_Na + V[j]) * g_bna;
    
    dNa_i = -(3 * i_NaCa + 3 * i_NaK + i_Na + i_b_Na) * conc_clamp * Cm / (F * Vol_c);
    
    E_Ks = R * T * log((Na_o * P_kna + K_o) / (Na_i[j] * P_kna + K_i[j])) / F;
    
    i_Ks = pow(Xs[j], 2) * (-E_Ks + V[j]) * g_Ks;
    
    dK_i = -(-2 * i_NaK + i_K1 + i_Kr + i_Ks + i_p_K + i_to + i_Stim) * conc_clamp * Cm / (F * Vol_c);
    
    i_tot = i_CaL + i_b_Ca + i_K1 + i_Kr + i_Ks + i_Na + i_NaCa + i_NaK + i_b_Na + i_p_Ca + i_p_K + i_to + i_Stim;
    
    dV = -i_tot;


  //----Rush-Larsen --------
  Xs[j+1] = xs_inf + (Xs[j] - xs_inf)*exp(-dt/tau_xs);

  s[j+1] = s_inf + (s[j] - s_inf)*exp(-dt/tau_s);

  r[j+1] = r_inf + (r[j] - r_inf)*exp(-dt/tau_r);

  m[j+1] = m_inf + (m[j] - m_inf)*exp(-dt/tau_m);

  h[j+1] = h_inf + (h[j] - h_inf)*exp(-dt/tau_h);

  jj[j+1] = j_inf + (jj[j] - j_inf)*exp(-dt/tau_j);

  d[j+1] = d_inf + (d[j] - d_inf)*exp(-dt/tau_d);

  f[j+1] = f_inf + (f[j] - f_inf)*exp(-dt/tau_f);

  f2[j+1] = f2_inf + (f2[j] - f2_inf)*exp(-dt/tau_f2);


//------ dynamics ----------




  V[j+1] = V[j] + dt * dV;

  Ca_i[j+1] = Ca_i[j] + dt * dCa_i;

  Cr1[j+1] = Cr1[j] + dt * dCr1;

  Cr2[j+1] = Cr2[j] + dt * dCr2;

  Cr3[j+1] = Cr3[j] + dt * dCr3;

  Or4[j+1] = Or4[j] + dt * dOr4;

  Ir5[j+1] = Ir5[j] + dt * dIr5;

  BCr1[j+1] = BCr1[j] + dt * dBCr1;

  BCr2[j+1] = BCr2[j] + dt * dBCr2;

  BCr3[j+1] = BCr3[j] + dt * dBCr3;

  BOr4[j+1] = BOr4[j] + dt * dBOr4;

  BIr5[j+1] = BIr5[j] + dt * dBIr5;

  fCass[j+1] = fCass[j] + dt * dfCass;

  Ca_SR[j+1] = Ca_SR[j] + dt * dCa_SR;

  Ca_ss[j+1] = Ca_ss[j] + dt * dCa_ss;

  R_prime[j+1] = R_prime[j] + dt * dR_prime;

  Na_i[j+1] = Na_i[j] + dt * dNa_i;

  K_i[j+1] = K_i[j] + dt * dK_i;
  
  
   }

   plt::figure();
   plt::plot(V);
   plt::show();
}
