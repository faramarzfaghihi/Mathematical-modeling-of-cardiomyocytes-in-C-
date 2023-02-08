// C++ code for model : Ten Tusscher 2004
// Written by Faramarz Faghihi using Chaste_code from CellML model

//------------------------------- Headers -------------------------------------

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
// Computed variables
//------------------------------------------------------------------------------

double I_ext; // dimensionless (in I_ex)

double alpha_fCa; // dimensionless (in L_type_Ca_current_fCa_gate)
double beta_fCa; // dimensionless (in L_type_Ca_current_fCa_gate)
double gama_fCa; // dimensionless (in L_type_Ca_current_fCa_gate)
//----------
double fCa_inf; // dimensionless (in L_type_Ca_current_fCa_gate)
//----------
double d_fCa; // per_millisecond (in L_type_Ca_current_fCa_gate)

double Ca_i_bufc; // dimensionless (in calcium_dynamics)
double Ca_sr_bufsr; // dimensionless (in calcium_dynamics)
//-----------
double g_inf; // dimensionless (in calcium_dynamics)
double tau_g = 2;
//-----------
double i_leak; // millimolar_per_millisecond (in calcium_dynamics)
double i_rel; // millimolar_per_millisecond (in calcium_dynamics)
double i_up; // millimolar_per_millisecond (in calcium_dynamics)
double d_g; // per_millisecond (in calcium_dynamics)
double i_p_Ca; // picoA_per_picoF (in calcium_pump_current)
double g_leak; // one_over_ohm (in g_lea)
double i_CaL; // picoA_per_picoF (in L_type_Ca_current)

double alpha_d; // dimensionless (in L_type_Ca_current_d_gate)
double beta_d; // dimensionless (in L_type_Ca_current_d_gate)
//------------
double d_inf; // dimensionless (in L_type_Ca_current_d_gate)
double gamma_d; // millisecond (in L_type_Ca_current_d_gate)
double tau_d; // millisecond (in L_type_Ca_current_d_gate)
//----------
double f_inf; // dimensionless (in L_type_Ca_current_f_gate)
double tau_f; // millisecond (in L_type_Ca_current_f_gate)
//--------
double alpha_h; // per_millisecond (in fast_sodium_current_h_gate)
double beta_h; // per_millisecond (in fast_sodium_current_h_gate)
//-----------
double h_inf; // dimensionless (in fast_sodium_current_h_gate)
double tau_h; // millisecond (in fast_sodium_current_h_gate)
//----------
double alpha_j; // per_millisecond (in fast_sodium_current_j_gate)
double beta_j; // per_millisecond (in fast_sodium_current_j_gate)
//-----------
double j_inf; // dimensionless (in fast_sodium_current_j_gate)
double tau_j; // millisecond (in fast_sodium_current_j_gate)
//-----------
double alpha_m; // dimensionless (in fast_sodium_current_m_gate)
double beta_m; // dimensionless (in fast_sodium_current_m_gate)
//-----------
double m_inf; // dimensionless (in fast_sodium_current_m_gate)
double tau_m; // millisecond (in fast_sodium_current_m_gate)
//----------

double alpha_xr1; // dimensionless (in rapid_time_dependent_potassium_current_Xr1_gate)
double beta_xr1; // dimensionless (in rapid_time_dependent_potassium_current_Xr1_gate)
//-----------
double tau_xr1; // millisecond (in rapid_time_dependent_potassium_current_Xr1_gate)
double xr1_inf; // dimensionless (in rapid_time_dependent_potassium_current_Xr1_gate)
//----------
double alpha_xr2; // dimensionless (in rapid_time_dependent_potassium_current_Xr2_gate)
double beta_xr2; // dimensionless (in rapid_time_dependent_potassium_current_Xr2_gate)
//----------
double tau_xr2; // millisecond (in rapid_time_dependent_potassium_current_Xr2_gate)
double xr2_inf; // dimensionless (in rapid_time_dependent_potassium_current_Xr2_gate)
//------------
double E_Ca; // millivolt (in reversal_potentials)
double i_b_Ca; // picoA_per_picoF (in calcium_background_current)
double E_K; // millivolt (in reversal_potentials)
double alpha_K1; // dimensionless (in inward_rectifier_potassium_current)
double beta_K1; // dimensionless (in inward_rectifier_potassium_current)
//----------
double xK1_inf; // dimensionless (in inward_rectifier_potassium_current)
//----------
double i_K1; // picoA_per_picoF (in inward_rectifier_potassium_current)
double i_p_K; // picoA_per_picoF (in potassium_pump_current)
double i_Kr; // picoA_per_picoF (in rapid_time_dependent_potassium_current)

double alpha_xs; // dimensionless (in slow_time_dependent_potassium_current_Xs_gate)
double beta_xs; // dimensionless (in slow_time_dependent_potassium_current_Xs_gate)
//-----------
double tau_xs; // millisecond (in slow_time_dependent_potassium_current_Xs_gate)
double xs_inf; // dimensionless (in slow_time_dependent_potassium_current_Xs_gate)
//----------
double E_Ks; // millivolt (in reversal_potentials)
double E_Na; // millivolt (in reversal_potentials)
double i_Na; // picoA_per_picoF (in fast_sodium_current)
double i_Ks; // picoA_per_picoF (in slow_time_dependent_potassium_current)
double i_b_Na; // picoA_per_picoF (in sodium_background_current)
double i_NaCa; // picoA_per_picoF (in sodium_calcium_exchanger_current)
double i_NaK; // picoA_per_picoF (in sodium_potassium_pump_current)

double i_bck; // dimensionless (in i_bc)
double i_leak_comp; // dimensionless (in i_leak_com)
double i_inj; // dimensionless (in i_in)
double i_Stim; // picoA_per_picoF (in membrane)

//-----------
double r_inf; // dimensionless (in transient_outward_current_r_gate)
double tau_r; // millisecond (in transient_outward_current_r_gate)
//-----------
double i_to; // picoA_per_picoF (in transient_outward_current)
//
double s_inf; // dimensionless (in transient_outward_current_s_gate)
double tau_s; // millisecond (in transient_outward_current_s_gate)
//
    // Constants
    //------------------------------------------------------------------------------

double A0_bck = 1.0278;
double Ampl_gain = 1;
double Cext = 1;
double   E_l = 1;
double I = 1;
double R_seal = 1;
double Scale_bck = 1;
double Scaling = 1;
double k_bck = 0.098599999999999993;
double leak_comp_perc = 1;
double V_leak = 8.0000000000000007e-5;
double a_rel = 0.016463999999999999;
double Vmax_up = 0.00042499999999999998;
double conc_clamp = 1;
double Ca_o = 2;
double K_o = 5.4000000000000004;
double Na_o = 140;
double g_CaL = 0.000175;
double tau_fCa = 2;
double g_bca = 0.00059199999999999997;
double g_bna = 0.00029;
double g_pCa = 0.82499999999999996;
double Cm = 0.185;
double g_Na = 14.837999999999999;
double perc_reduced_inact_for_IpNa = 0;
double shift_INa_inact = 0;
double g_K1 = 5.4050000000000002;
double g_pK = 0.0146;
double g_Kr = 0.096000000000000002;
double g_Ks = 0.062;
double K_NaCa = 1000;
double P_NaK = 1.3620000000000001;
double stim_amplitude = -52;
double stim_duration = 1;
double stim_start = 100;
double stim_period = 1000;
double g_to = 0.29399999999999998;

double Buf_c = 0.14999999999999999;
double Buf_sr = 10;
double K_buf_c = 0.001;
double K_buf_sr = 0.29999999999999999;
double K_up = 0.00025000000000000001;
double V_sr = 0.0010939999999999999;
double b_rel = 0.25;
double c_rel = 0.0082319999999999997;
double K_pCa = 0.00050000000000000001; 
    
double F = 96485.341499999995;
double R = 8314.4719999999998;
double T = 310;
double V_c = 0.016403999999999998;
double P_kna = 0.029999999999999999;
double K_sat = 0.10000000000000001;
double Km_Ca = 1.3799999999999999;
double Km_Nai = 87.5;
double alpha = 2.5;
double gamm = 0.349;
double K_mNa = 40;
double K_mk = 1;

//------------------------------------------------------------------------------------------
double Xr1  = 0.0; //  2 Xr1 (dimensionless) (in rapid_time_dependent_potassium_current_Xr1_gate)
double Xr2  = 1.0; //  3 Xr2 (dimensionless) (in rapid_time_dependent_potassium_current_Xr2_gate)
double Xs   = 0.0; // 4 Xs (dimensionless) (in slow_time_dependent_potassium_current_Xs_gate)
double m    = 0.0; // 5 m (dimensionless) (in fast_sodium_current_m_gate)
double h    = 0.75; // 6 h (dimensionless) (in fast_sodium_current_h_gate)
double jj    = 0.75; // 7 j (dimensionless) (in fast_sodium_current_j_gate)
double f    = 1.0; //  9 f (dimensionless) (in L_type_Ca_current_f_gate)
double s    = 1.0; // 11 s (dimensionless) (in transient_outward_current_s_gate)
double r    = 0.0; // 12 r (dimensionless) (in transient_outward_current_r_gate)

//------------------------------------------------------------------------------
    
double d_dt_Ca_i;
double d_dt_g;
double d_dt_fCa;
double d_dt_Ca_SR;
double d_dt_Na_i;
double d_dt_K_i;
double dv_dt;

double dXr1_dt;
double dXr2_dt;
double dXs_dt;
double dm_dt;
double dh_dt;
double djj_dt;
double dd_dt;
double df_dt;
double dfCa_dt;
double ds_dt;
double dr_dt;  
//------------------------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------------------------
int main()
{
    std::vector<double> t(1000);
    double time;
    std::vector<double> v(t.size());
    std::vector<double> Xr1(t.size());
    std::vector<double> Xr2(t.size());
    std::vector<double> Xs(t.size());
    std::vector<double> m(t.size());
    std::vector<double> h(t.size());
    std::vector<double> jj(t.size());
    std::vector<double> d(t.size());
    std::vector<double> f(t.size());
    std::vector<double> fCa(t.size());
    std::vector<double> s(t.size());
    std::vector<double> r(t.size());
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
    std::vector<double> Ca_i(t.size());
    std::vector<double> Ca_SR(t.size());
    std::vector<double> g(t.size());
    std::vector<double> Na_i(t.size());
    std::vector<double> K_i(t.size());

//--------------------------------------------------------------------------
 double dt = 0.052;                                               // Integration steps  
//---------------------------------------------------------------------------------------

 for (size_t j = 0; j < t.size(); j++)
  {

    v[0]     = -86.2;    // 0 V (millivolt) (in membrane)
    Ca_i[0]  = 0.0002;  //1  Ca_i (millimolar) (in calcium_dynamics)
    fCa[0]   = 1;       // 10 fCa (dimensionless) (in L_type_Ca_current_fCa_gate)
    Ca_SR[0] = 0.2;     // 13 Ca_SR (millimolar) (in calcium_dynamics)
    g[0] = 1;           //g (dimensionless) (in calcium_dynamics)
    Na_i[0]  = 11.6;   //Na_i (millimolar) (in sodium_dynamics)
    K_i[0]   = 138.3;  //K_i (millimolar) (in potassium_dynamics)
//-----------------------------------------------------------------------------------------

 //--------------------------------------------------------------------------------------
    I_ext = I / (Ampl_gain * Cext);
    alpha_fCa = 1 / (1 + 8.034023767017109e+27 * pow(Ca_i[j], 8));
    beta_fCa = 0.10000000000000001 / (1 + exp(-5 + 10000 * Ca_i[j]));
    gama_fCa = 0.20000000000000001 / (1 + exp(-0.9375 + 1250 * Ca_i[j]));
//-----------------------------------------------------------------------------
    fCa_inf = 0.15753424657534248 + 0.68493150684931503 * alpha_fCa + 0.68493150684931503 * beta_fCa + 0.68493150684931503 * gama_fCa;
    //
    d_fCa = (-fCa[j] + fCa_inf) / tau_fCa;
//-----------------------------------------------------------------------------------

    Ca_i_bufc = 1 / (1 + Buf_c * K_buf_c / pow((Ca_i[j] + K_buf_c), 2));
    Ca_sr_bufsr = 1 / (1 + Buf_sr * K_buf_sr / pow((Ca_SR[j] + K_buf_sr), 2));
      //-----------------
      if (Ca_i[j] < 0.00035)
      
          {g_inf = 1 / (1 + 5.439910241481016e+20 * pow(Ca_i[j], 6));}
       else 
       
          {g_inf = 1 / (1 + 1.9720198874049176e+55 * pow(Ca_i[j], 16));}
       
     //-----------------

    i_leak = (-Ca_i[j] + Ca_SR[j]) * V_leak;
    i_rel = (pow(Ca_SR[j], 2) * a_rel / (pow(Ca_SR[j], 2) + pow(b_rel, 2)) + c_rel) * g[j] * d[j];
    i_up = Vmax_up / (1 + pow(K_up, 2) / pow(Ca_i[j], 2));
    //-------------------
    d_g = (-g[j] + g_inf) / tau_g;
    //------------------
//-------------------------------------------------------------------------------------
    i_p_Ca = Ca_i[j] * g_pCa / (Ca_i[j] + K_pCa);
    g_leak = 1 / R_seal;
    //-----------------
      if ((v[0] > -60) && (fCa[j] < fCa_inf)) 
          {d_dt_fCa = 0;}
       
       else
        
          {d_dt_fCa = d_fCa;}
       

      if ((v[0] > -60) && (g[j] < g_inf))
         {d_dt_g = 0;}
       
        else
       
           {d_dt_g = d_g;}
       
    //-----------------
    if ((v[0] >= -1.3356880329847825e-6) && (v[0] <= 1.3356880329847825e-6))
        {i_CaL = 374338.90822745475 * (1.3356880329847825e-6 + v[0]) * (0.019297068299999998 * (-0.34100000000000003 * Ca_o + Ca_i[j] * exp(9.9999999999999995e-8)) * g_CaL * fCa[j] * d[j] * f[j] / (-1 + exp(9.9999999999999995e-8)) + 0.019297068299999998 * (-0.34100000000000003 * Ca_o + Ca_i[j] * exp(-9.9999999999999995e-8)) * g_CaL * fCa[j] * d[j] * f[j] / (-1 + exp(-9.9999999999999995e-8))) - 0.019297068299999998 * (-0.34100000000000003 * Ca_o + Ca_i[j] * exp(-9.9999999999999995e-8)) * g_CaL * fCa[j] * d[j] * f[j] / (-1 + exp(-9.9999999999999995e-8));}
      
     else 
      
        {i_CaL = 14447.286958825251 * (-0.34100000000000003 * Ca_o + Ca_i[j] * exp(0.074867781645490947 * v[0])) * g_CaL * v[0] * fCa[j] * d[j] * f[j] / (-1 + exp(0.074867781645490947 * v[0]));}
      
   //----------------
//----------------------------------------------------------------------------------------------
    alpha_d = 0.25 + 1.3999999999999999 / (1 + exp(-2.6923076923076925 - 0.076923076923076927 * v[0]));
    beta_d = 1.3999999999999999 / (1 + exp(1 + 0.20000000000000001 * v[0]));
    //--------------
    d_inf = 1 / (1 + exp(-0.66666666666666663 - 0.13333333333333333 * v[0]));
    //--------------
    gamma_d = 1 / (1 + exp(2.5 - 0.050000000000000003 * v[0]));

    //-------------
    tau_d = alpha_d * beta_d + gamma_d;
    d[j] = (-d[j] + d_inf) / tau_d;
    //-------------
//----------------------------------------------------------------------------------------------
    
    f_inf = 1 / (1 + exp(2.8571428571428572 + 0.14285714285714285 * v[0]));
    tau_f = 80 + 165 / (1 + exp(2.5 - 0.10000000000000001 * v[0])) + 1125 * exp(-3.0375000000000001 * pow((1 + 0.037037037037037035 * v[0]), 2));
    f[j] = (-f[j] + f_inf) / tau_f;
    
//----------------------------------------------------------------------------------------------
    if (v[0] < -40)
        {alpha_h = 0.057000000000000002 * exp(-11.764705882352942 - 0.14705882352941177 * v[0]);}
    
     else
     
        {alpha_h = 0;}
     

    if (v[0] < -40)
        {beta_h = 310000 * exp(0.34849999999999998 * v[0]) + 2.7000000000000002 * exp(0.079000000000000001 * v[0]);}
    
     else
      
        {beta_h = 5.9230769230769234 / (1 + exp(-0.96036036036036043 - 0.0900900900900901 * v[0]));}
     
    //------------------
    h_inf = 0.01 * perc_reduced_inact_for_IpNa + (1 - 0.01 * perc_reduced_inact_for_IpNa) / pow((1 + exp(9.6298788694481825 + 0.13458950201884254 * v[0] - 0.13458950201884254 * shift_INa_inact)), 2);
    tau_h = 1 / (alpha_h + beta_h);
    h[j] = (-h[j] + h_inf) / tau_h;
    //------------------
//------------------------------------------------------------------------------------------------
    if (v[0] < -40)
        {alpha_j = (37.780000000000001 + v[0]) * (-25428 * exp(0.24440000000000001 * v[0]) - 6.9480000000000002e-6 * exp(-0.043909999999999998 * v[0])) / (1 + exp(24.640530000000002 + 0.311 * v[0]));}
    
     else  
        {alpha_j = 0;}
     
    if (v[0] < -40)
        {beta_j = 0.024240000000000001 * exp(-0.01052 * v[0]) / (1 + exp(-5.5312920000000005 - 0.13780000000000001 * v[0]));}
     
      else  
        {beta_j = 0.59999999999999998 * exp(0.057000000000000002 * v[0]) / (1 + exp(-3.2000000000000002 - 0.10000000000000001 * v[0]));}
     

    //----------------   
    j_inf = 0.01 * perc_reduced_inact_for_IpNa + (1 - 0.01 * perc_reduced_inact_for_IpNa) / pow((1 + exp(9.6298788694481825 + 0.13458950201884254 * v[0] - 0.13458950201884254 * shift_INa_inact)), 2);
    tau_j = 1 / (alpha_j + beta_j);
    jj[j] = (-jj[j] + j_inf) / tau_j;
    //----------------
//---------------------------------------------------------------------------------------------------
    alpha_m = 1 / (1 + exp(-12 - 0.20000000000000001 * v[0]));
    beta_m = 0.10000000000000001 / (1 + exp(7 + 0.20000000000000001 * v[0])) + 0.10000000000000001 / (1 + exp(-0.25 + 0.0050000000000000001 * v[0]));
    //-------------
    m_inf = 1 / pow((1 + exp(-6.2967884828349945 - 0.11074197120708749 * v[0])), 2);
    tau_m = alpha_m * beta_m;
    m[j] = (-m[j] + m_inf) / tau_m;
    //-------------
//--------------------------------------------------------------------------------------------------
    i_bck = A0_bck * Scale_bck / ((1 + exp(-v[0] * k_bck)) * Cext);
    i_leak_comp = (-E_l + v[0]) * g_leak * leak_comp_perc / (100 * Cext);
    i_inj = (-i_bck - i_leak_comp + I_ext) * Scaling;
//--------------------------------------------------------------------------------------------------
    d_dt_Ca_SR = (-i_leak - i_rel + i_up) * Ca_sr_bufsr * V_c / V_sr;

    if ((stim_start <= -stim_period * floor(j / stim_period) + j) && (stim_duration + stim_start >= -stim_period * floor(j / stim_period) + j))
        {i_Stim = stim_amplitude;}
    
     else 
        {i_Stim = 0;}
      

//--------------------------------------------------------------------------------------------------
    alpha_xr1 = 450 / (1 + exp(-4.5 - 0.10000000000000001 * v[0]));
    beta_xr1 = 6 / (1 + exp(2.6086956521739131 + 0.086956521739130432 * v[0]));
    //-------------
    tau_xr1 = alpha_xr1 * beta_xr1;
    xr1_inf = 1 / (1 + exp(-3.7142857142857144 - 0.14285714285714285 * v[0]));
    Xr1[j] = (-Xr1[j] + xr1_inf) / tau_xr1;
    //--------------
//--------------------------------------------------------------------------------------------------
    alpha_xr2 = 3 / (1 + exp(-3 - 0.050000000000000003 * v[0]));
    beta_xr2 = 1.1200000000000001 / (1 + exp(-3 + 0.050000000000000003 * v[0]));
    //-------------
    tau_xr2 = alpha_xr2 * beta_xr2;
    xr2_inf = 1 / (1 + exp(3.6666666666666665 + 0.041666666666666664 * v[0]));
    Xr2[j] = (-Xr2[j] + xr2_inf) / tau_xr2;
    //-------------
//--------------------------------------------------------------------------------------------------
    E_Ca = 0.5 * R * T * log(Ca_o / Ca_i[j]) / F;
    i_b_Ca = (-E_Ca + v[0]) * g_bca;
    E_K = R * T * log(K_o / K_i[j]) / F;
    alpha_K1 = 0.10000000000000001 / (1 + exp(-12 + 0.059999999999999998 * v[0] - 0.059999999999999998 * E_K));
    beta_K1 = (3 * exp(0.02 + 0.00020000000000000001 * v[0] - 0.00020000000000000001 * E_K) + exp(-1 + 0.10000000000000001 * v[0] - 0.10000000000000001 * E_K)) / (1 + exp(0.5 * E_K - 0.5 * v[0]));
    //-------------PROBLEMEMATIC BLOCK OF CODE ??-------------------

    //xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
    
    //-------------PROBLEMEMATIC BLOCK OF CODE ??-------------------
//---------------------------------------------------------------------------------------------------    
    i_K1 = 0.43033148291193518 * sqrt(K_o) * (-E_K + v[0]) * g_K1;
    i_p_K = (-E_K + v[0]) * g_pK / (1 + exp(4.1806020066889626 - 0.16722408026755853 * v[0]));
    i_Kr = 0.43033148291193518 * sqrt(K_o) * (-E_K + v[0]) * Xr1[j] * Xr2[j] * g_Kr;
    alpha_xs = 1100 / sqrt(1 + exp(-1.6666666666666667 - 0.16666666666666666 * v[0]));
    beta_xs = 1 / (1 + exp(-3 + 0.050000000000000003 * v[0]));
    //----------------
    tau_xs = alpha_xs * beta_xs;
    xs_inf = 1 / (1 + exp(-0.35714285714285715 - 0.071428571428571425 * v[0]));
    Xs[j] = (-Xs[j] + xs_inf) / tau_xs;
    //---------------
 //------------------------------------------------------------------------------------------------   
    
    E_Ks = R * T * log((P_kna * Na_o + K_o) / (Na_i[j] * P_kna + K_i[j])) / F;
    E_Na = R * T * log(Na_o / Na_i[j]) / F;
    i_Na = pow(m[j], 3) * (-E_Na + v[0]) * h[j] * jj[j] * g_Na;
    i_Ks = pow(Xs[j], 2) * (-E_Ks + v[0]) * g_Ks;
    i_b_Na = (-E_Na + v[0]) * g_bna;
    i_NaCa = (pow(Na_i[j], 3) * Ca_o * exp(v[0] * F * gamm / (R * T)) - pow(Na_o, 3) * Ca_i[j] * alpha * exp((-1 + gamm) * v[0] * F / (R * T))) * K_NaCa / ((1 + K_sat * exp((-1 + gamm) * v[0] * F / (R * T))) * (pow(Km_Nai, 3) + pow(Na_o, 3)) * (Ca_o + Km_Ca));
    
    d_dt_Ca_i = (-i_up - 0.5 * (-2 * i_NaCa + i_CaL + i_b_Ca + i_p_Ca) * Cm / (F * V_c) + i_leak + i_rel) * Ca_i_bufc;

//--------------------------------------------------------------------------------------------------
    i_NaK = Na_i[j] * K_o * P_NaK / ((Na_i[j] + K_mNa) * (K_o + K_mk) * (1 + 0.035299999999999998 * exp(-v[0] * F / (R * T)) + 0.1245 * exp(-0.10000000000000001 * v[0] * F / (R * T))));
    
    d_dt_Na_i = -(3 * i_NaCa + 3 * i_NaK + i_Na + i_b_Na) * Cm * conc_clamp / (F * V_c);

   //---------------
    r_inf = 1 / (1 + exp(3.3333333333333335 - 0.16666666666666666 * v[0]));
    tau_r = 0.80000000000000004 + 9.5 * exp(-0.88888888888888884 * pow((1 + 0.025000000000000001 * v[0]), 2));
    r[j] = (-r[j] + r_inf) / tau_r;
    //--------------
//---------------------------------------------------------------------------------------------------
    i_to = (-E_K + v[0]) * s[j] * r[j] * g_to;
    
    v[j] = -i_CaL - i_b_Ca - i_p_Ca - i_Na - i_inj - i_K1 - i_Stim - i_p_K - i_Kr - i_Ks - i_b_Na - i_NaCa - i_NaK - i_to;
    
     d_dt_K_i = -(-2 * i_NaK + i_K1 + i_Stim + i_p_K + i_Kr + i_Ks + i_to) * Cm * conc_clamp / (F * V_c);
    //-------------
    s_inf = 1 / (1 + exp(4 + 0.20000000000000001 * v[0]));
    tau_s = 3 + 5 / (1 + exp(-4 + 0.20000000000000001 * v[0])) + 85 * exp(-6.328125 * pow((1 + 0.022222222222222223 * v[0]), 2));
    s[j] = (-s[j] + s_inf) / tau_s;
    //-------------

//--------------------------------------------------------------------------------------   
    dv_dt =  v[j] / Cm;     
    
    Xr1[j+1] = (xr1_inf + ( Xr1[j] - xr1_inf )*exp( -dt*tau_xr1 ));
    
    Xr2[j+1] = (xr2_inf + ( Xr2[j] - xr2_inf )*exp(-dt*tau_xr2 ));
    
    Xs[j+1] = (xs_inf + ( Xs[j] - xs_inf )*exp( -dt*tau_xs ));
    
    m[j+1] = (m_inf + ( m[j] - m_inf )*exp( -dt*tau_m ));
    
    h[j+1] = (h_inf + ( h[j] - h_inf )*exp( -dt*tau_h ));
    
    jj[j+1] = (j_inf + ( jj[j] - j_inf )*exp( -dt*tau_j ));
    
    d[j+1] = (d_inf + ( d[j] - d_inf )*exp( -dt*tau_d ));
    
    f[j+1] = (f_inf + ( f[j] - f_inf )*exp( -dt*tau_f ));
    
    fCa[j+1] = (fCa_inf + ( fCa[j] - fCa_inf )*exp( -dt*tau_fCa));
    
    s[j+1] = (s_inf + ( s[j] - s_inf )*exp( -dt*tau_s ));
    
    r[j+1] = (r_inf + ( r[j] - r_inf )*exp( -dt*tau_r ));

// Integration Euler-----------------------------------------------------------------------
    v[j+1] = v[j] + dv_dt * dt;
//------------------------------------------------------------------------------------------
    Ca_i[j+1] = Ca_i[j] + d_dt_Ca_i * dt;
    Ca_SR[j+1] = Ca_SR[j] + d_dt_Ca_SR * dt;
    g[j+1] = g[j] + d_dt_g * dt;
    Na_i[j+1] = Na_i[j] + d_dt_Na_i * dt;
    K_i[j+1] = K_i[j] + d_dt_K_i * dt;

//-----------------------------------------------------------------------------------------        
    time = time + dt;
 
//-----------------------------------------------------------------------------------------

 
     
  }
 
//---------------------------------------------------------------------------------------
  
 plt::figure(); 
 plt::plot(v);
}
 
