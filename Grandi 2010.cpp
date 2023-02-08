// C++ code for model : Grandi et al., 2010
// Written by Faramarz Faghihi using Chaste_code from CellML model



#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <array>
#include "matplotlibcpp.h"
using namespace std;
#include <vector>
namespace plt = matplotlibcpp;
double stim_amplitude = -9.5;
double stim_duration = 10;
double stim_start = 100;
double stim_period =5000;

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------




//------------------------------------------------------------------

   double A0_bck = 1.0278;
   double I = 1;
   double SR_leak_max = 5.3480000000000003e-6;
   double ks = 25;
   double Vmax_SRCaP = 0.0053114;
   double conc_clamp = 1;
   double Cli = 15;
   double Cao = 1.8;
   double Clo = 150;
   double Ko = 5.4000000000000004;
   double Nao = 140;
   double Frdy = 96485;
   double R = 8314;
   double pCa_max = 0.5;
   double GNa = 23;
   double perc_reduced_inact_for_IpNa = 0;
   double shift_INa_inact = 0;
   double GtoFast_factor = 1;
   double g_ki_factor = 0.34999999999999998;
   double gkr_max = 0.035000000000000003;
   double GKs_factor = 1;
   double GtoSlow_factor = 1;
   double IbarNCX = 4.5;
   double IbarNaK = 1.8;
   double Cmem = 1.3809999999999999e-10;
   double Temp = 310;
//------------------------------------------------------------------------------
   
   double Ampl_gain = 1;
   double Cext = 1;
   double E_l = 1;
   double gks_junc = 0.0035000000000000001;
   double gks_sl = 0.0035000000000000001;
   double R_seal = 1;
   double MaxSR = 15;
   double MinSR = 1;
    double Scale_bck = 1;
   double Scaling = 1;
   double k_bck = 0.098599999999999993;
   double leak_comp_perc = 1;
   double Bmax_CaM = 0.024;
   double Bmax_Naj = 7.5609999999999999;
   double Bmax_Nasl = 1.6499999999999999;
   double Bmax_SR = 0.017100000000000001;
   double Bmax_TnChigh = 0.14000000000000001;
   double Bmax_TnClow = 0.070000000000000007;
   double Bmax_myosin = 0.14000000000000001;
   double Fjunc = 0.11;
   double Fjunc_CaL = 0.90000000000000002;
   double GCaB = 0.00055130000000000001;
   double GClB = 0.0089999999999999993;
   double GClCa = 0.0548125;
   double GNaB = 0.00059699999999999998;
   double IbarSLCaP = 0.067299999999999999;
   double J_ca_juncsl = 8.2413000000000004e-13;
   double J_ca_slmyo = 3.2742999999999999e-12;
   double J_na_juncsl = 1.8312999999999999e-14;
   double J_na_slmyo = 1.6385999999999999e-12;
   double KdClCa = 0.10000000000000001;
   double Kdact = 0.00014999999999999999;
   double KmCai = 0.0035899999999999999;
   double KmCao = 1.3;
   double KmKo = 1.5;
   double KmNai = 12.289999999999999;
   double KmNaip = 11;
   double KmNao = 87.5;
   double KmPCa = 0.00050000000000000001;
   double Kmf = 0.00024600000000000002;
   double Kmr = 1.7;
   double Mgi = 1;
   double Q10CaL = 1.8;
   double Q10NCX = 1.5700000000000001;
   double Q10SLCaP = 2.3500000000000001;
   double Q10SRCaP = 2.6000000000000001;
   double cellLength = 100;
   double cellRadius = 10.25;
   double ec50SR = 0.45000000000000001;
   double epi = 1;
   double gkp = 0.002;
   double hillSRCaP = 1.7869999999999999;
   double kiCa = 0.5;
   double kim = 0.0050000000000000001;
   double koCa = 10;
   double koff_cam = 0.23799999999999999;
   double koff_csqn = 65;
   double koff_myoca = 0.00046000000000000001;
   double koff_myomg = 5.7000000000000003e-5;
   double koff_na = 0.001;
   double koff_slh = 0.029999999999999999;
   double koff_sll = 1.3;
   double koff_sr = 0.059999999999999998;
   double koff_tnchca = 3.1999999999999999e-5;
   double koff_tnchmg = 0.0033300000000000001;
   double koff_tncl = 0.019599999999999999;
   double kom = 0.059999999999999998;
   double kon_cam = 34;
   double kon_csqn = 100;
   double kon_myoca = 13.800000000000001;
   double kon_myomg = 0.015699999999999999;
   double kon_na = 0.0001;
   double kon_slh = 100;
   double kon_sll = 100;
   double kon_sr = 100;
   double kon_tnchca = 2.3700000000000001;
   double kon_tnchmg = 0.0030000000000000001;
   double kon_tncl = 32.700000000000003;
   double ksat = 0.32000000000000001;
   double nu = 0.27000000000000002;
   double pNaK = 0.018329999999999999;

int main()
{

double dV;
double dCa_i;
double df_Ca_Bj;
double df_Ca_Bsl;
double dRy_Rr;
double dRy_Ro;
double dRy_Ri;
double dNa_Bj;
double dNa_Bsl;
double dTn_CL;
double dTn_CHm;
double dTn_CHc;
double dCaM;
double dMyo_c;
double dMyo_m;
double dSRB;
double dSLL_j;
double dSLL_sl;
double dSLH_j;
double dSLH_sl;
double dCsqn_b;
double dCa_sr;
double dNa_j;
double dNa_sl;
double dNa_i;
double dK_i;
double dCa_j;
double dCa_sl;
//------------------------------------------------------------------------------

double fcaCaMSL; // dimensionless (in I_Ca)
double fcaCaj; // dimensionless (in I_Ca)
double I_ext; // dimensionless (in I_ex)
double RI; // mM (in SR_Fluxes)
double J_SRleak; // mM_per_msec (in SR_Fluxes)
double g_leak; // one_over_ohm (in g_lea)
double dss; // dimensionless (in I_Ca)
double fss; // dimensionless (in I_Ca)
double taud; // msec (in I_Ca)
double tauf; // msec (in I_Ca)
double kp_kp; // dimensionless (in I_Kp)
double rkr; // dimensionless (in I_Kr)
double tauxr; // msec (in I_Kr)
double xrss; // dimensionless (in I_Kr)
double tauxs; // msec (in I_Ks)
double xsss; // dimensionless (in I_Ks)
double ah; // dimensionless (in I_Na)
double aj; // dimensionless (in I_Na)
double bh; // dimensionless (in I_Na)
double bj; // dimensionless (in I_Na)
double hss; // dimensionless (in I_Na)
double jss; // dimensionless (in I_Na)
double mss; // dimensionless (in I_Na)
double tauh; // msec (in I_Na)
double tauj; // msec (in I_Na)
double taum; // msec (in I_Na)
double tauxtof; // msec (in I_to)
double tauxtos; // msec (in I_to)
double tauytof; // msec (in I_to)
double tauytos; // msec (in I_to)
double xtoss; // dimensionless (in I_to)
double ytoss; // dimensionless (in I_to)
double i_bck; // dimensionless (in i_bc)
double i_leak_comp; // dimensionless (in i_leak_com)
//double i_inj; // dimensionless (in i_in)

double i_Stim=0; // uA_per_uF (in membrane_potential)
double Fsl; // dimensionless (in parameters)
double GKs_total; // mS_per_uF (in I_Ks)
double Fsl_CaL; // dimensionless (in parameters)
double Ka_junc; // dimensionless (in I_NCX)
double Ka_sl; // dimensionless (in I_NCX)
double g_K1; // mS_per_uF (in I_Ki)
double gkr; // mS_per_uF (in I_Kr)
double s3_junc; // mM4 (in I_NCX)
double s3_sl; // mM4 (in I_NCX)
double sigma; // dimensionless (in I_NaK)
double ibark; // uA_per_uF (in I_Ca)
double ibarna_j; // uA_per_uF (in I_Ca)
double ibarna_sl; // uA_per_uF (in I_Ca)
double FoRT; // per_mV (in parameters)
double fnak; // dimensionless (in I_NaK)
double I_nak_junc; // uA_per_uF (in I_NaK)
double I_nak_sl; // uA_per_uF (in I_NaK)
double I_nak; // uA_per_uF (in I_NaK)
double Qpow; // dimensionless (in parameters)
double I_CaK; // uA_per_uF (in I_Ca)
double I_CaNa_junc; // uA_per_uF (in I_Ca)
double I_CaNa_sl; // uA_per_uF (in I_Ca)
double I_pca_junc; // uA_per_uF (in I_PCa)
double I_pca_sl; // uA_per_uF (in I_PCa)
double Vcell; // liter (in parameters)
double Vjunc; // liter (in parameters)
double Vmyo; // liter (in parameters)
double Bmax_SLhighj; // mM (in parameters)
double Bmax_SLlowj; // mM (in parameters)
double Vsl; // liter (in parameters)
double Bmax_SLhighsl; // mM (in parameters)
double Bmax_SLlowsl; // mM (in parameters)
double Vsr; // liter (in parameters)
double Bmax_Csqn; // mM (in parameters)
double kCaSR; // dimensionless (in SR_Fluxes)
double eca_junc; // mV (in parameters)
double I_cabk_junc; // uA_per_uF (in I_CaBK)
double eca_sl; // mV (in parameters)
double I_cabk_sl; // uA_per_uF (in I_CaBK)
double ecl; // mV (in parameters)
double I_ClCa_junc; // uA_per_uF (in I_ClCa)
double I_ClCa_sl; // uA_per_uF (in I_ClCa)
double I_ClCa; // uA_per_uF (in I_ClCa)
double I_Clbk; // uA_per_uF (in I_ClCa)
//double I_Cl_tot; // uA_per_uF (in membrane_potential)
double ek; // mV (in parameters)
double aki; // dimensionless (in I_Ki)
double bki; // dimensionless (in I_Ki)
double kiss; // dimensionless (in I_Ki)
double I_ki; // uA_per_uF (in I_Ki)
double I_kr; // uA_per_uF (in I_Kr)
double ena_junc; // mV (in parameters)
double I_Na_junc; // uA_per_uF (in I_Na)
double I_nabk_junc; // uA_per_uF (in I_NaBK)
double ena_sl; // mV (in parameters)
double I_Na_sl; // uA_per_uF (in I_Na)
double I_nabk_sl; // uA_per_uF (in I_NaBK)
double GtoFast; // mS_per_uF (in I_to)
double GtoSlow; // mS_per_uF (in I_to)
double I_tof; // uA_per_uF (in I_to)
double I_tos; // uA_per_uF (in I_to)
double I_to; // uA_per_uF (in I_to)
double I_kp_junc; // uA_per_uF (in I_Kp)
double I_kp_sl; // uA_per_uF (in I_Kp)
double I_kp; // uA_per_uF (in I_Kp)
double J_serca; // mM_per_msec (in SR_Fluxes)
double kiSRCa; // per_mM_per_msec (in SR_Fluxes)
double koSRCa; // per_mM2_per_msec (in SR_Fluxes)
double dNa_Bj_dt; // mM_per_msec (in Na_Buffers)
double dNa_Bsl_dt; // mM_per_msec (in Na_Buffers)
double J_CaB_junction; // mM_per_msec (in Junctional_and_SL_Ca_Buffers)
double J_CaB_sl; // mM_per_msec (in Junctional_and_SL_Ca_Buffers)
double J_CaB_cytosol; // mM_per_msec (in Cytosolic_Ca_Buffers)
double J_SRCarel; // mM_per_msec (in SR_Fluxes)
double s1_junc; // mM4 (in I_NCX)
double s1_sl; // mM4 (in I_NCX)
double s2_junc; // mM4 (in I_NCX)
double I_ncx_junc; // uA_per_uF (in I_NCX)
double s2_sl; // mM4 (in I_NCX)
double I_ncx_sl; // uA_per_uF (in I_NCX)
double I_Na_tot_junc; // uA_per_uF (in Na_Concentrations)
double I_Na_tot_sl; // uA_per_uF (in Na_Concentrations)
//double I_Na_tot; // uA_per_uF (in membrane_potential)
double ibarca_j; // uA_per_uF (in I_Ca)
double I_Ca_junc; // uA_per_uF (in I_Ca)
double I_Ca_tot_junc; // uA_per_uF (in Ca_Concentrations)
double ibarca_sl; // uA_per_uF (in I_Ca)
double I_Ca_sl; // uA_per_uF (in I_Ca)
double I_Ca_tot_sl; // uA_per_uF (in Ca_Concentrations)
//double I_Ca_tot; // uA_per_uF (in membrane_potential)
double eks; // mV (in I_Ks)
double I_ks; // uA_per_uF (in I_Ks)
//double I_K_tot; // uA_per_uF (in K_Concentration)
//double I_tot; // uA_per_uF (in membrane_potential)

//---------------------------------------------------------------------

std::vector<double> t(5000);
std::vector<double> V(t.size());

std::vector<double> I_tot(t.size());
std::vector<double> i_inj(t.size());
std::vector<double> I_K_tot(t.size());
std::vector<double> I_Ca_tot(t.size());
std::vector<double> I_Cl_tot(t.size());
std::vector<double> I_Na_tot(t.size());

std::vector<double> Ca_sl(t.size());
std::vector<double> Na_i(t.size());
std::vector<double> Na_j(t.size());
std::vector<double> Ca_sr(t.size());
std::vector<double> SLH_sl(t.size());
std::vector<double> SLL_sl(t.size());   
std::vector<double> SRB(t.size());
std::vector<double> CaM(t.size());
std::vector<double> Tn_CHm(t.size());
std::vector<double> Tn_CL(t.size());
std::vector<double> Na_Bj(t.size());
std::vector<double> f(t.size());
std::vector<double> d(t.size());
std::vector<double> y_to_f(t.size());
std::vector<double> x_to_f(t.size());
std::vector<double> y_to_s(t.size());
std::vector<double> x_to_s(t.size());
std::vector<double> m(t.size());
std::vector<double> h(t.size());
std::vector<double> jj(t.size());
std::vector<double> x_ks(t.size());
std::vector<double> x_kr(t.size());
std::vector<double> Ca_i(t.size());
 
std::vector<double> f_Ca_Bj(t.size()); 
std::vector<double> f_Ca_Bsl(t.size()); 
std::vector<double> Ry_Rr(t.size()); 
std::vector<double> Ry_Ro(t.size()); 
std::vector<double> Ry_Ri(t.size()); 
std::vector<double> Na_Bsl(t.size()); 
std::vector<double> Tn_CHc(t.size());  
std::vector<double> Myo_c(t.size()); 
std::vector<double> Myo_m(t.size()); 
std::vector<double> SLL_j(t.size());
std::vector<double> SLH_j(t.size());
std::vector<double> Csqn_b(t.size()); 
std::vector<double> Na_sl(t.size()); 
std::vector<double> K_i(t.size()); 
std::vector<double> Ca_j(t.size());   

//---------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // State variables
    //---------------------------------------------------------------------------
   

    V[0] = -61.4229700631461; // V_m (mV) (in membrane_potential)
    
    Ca_i[0] = 8.72745589849657e-05; // Ca_i (mM) (in Ca_Concentrations)

    m[0] = 0.00381858135062259; // m (dimensionless) (in I_Na)
    
    h[0] = 0.625086621635645; // h (dimensionless) (in I_Na)

    jj[0] = 0.62332507235506; // j (dimensionless) (in I_Na)

    x_kr[0] = 0.021733554982663; // x_kr (dimensionless) (in I_Kr)
    
    x_ks[0] = 0.00428981903391989; // x_ks (dimensionless) (in I_Ks)

    x_to_s[0] = 0.000441539203559411; // x_to_s (dimensionless) (in I_to)
 
    y_to_s[0] = 0.784875334693892; // y_to_s (dimensionless) (in I_to)

    x_to_f[0] = 0.000441531250866821; // x_to_f (dimensionless) (in I_to)

    y_to_f[0] = 0.999995817153572; // y_to_f (dimensionless) (in I_to)

    d[0] = 2.93982287251302e-06; // d (dimensionless) (in I_Ca)

    f[0] = 0.99511673495949; // f (dimensionless) (in I_Ca)

    f_Ca_Bj[0] = 0.0246142436477748; // f_Ca_Bj (dimensionless) (in I_Ca)

    f_Ca_Bsl[0] = 0.0152416826209301; // f_Ca_Bsl (dimensionless) (in I_Ca)

    Ry_Rr[0] = 0.891022230597263; // Ry_Rr (mM) (in SR_Fluxes)

    Ry_Ro[0] = 7.37484660389498e-07; // Ry_Ro (mM) (in SR_Fluxes)

    Ry_Ri[0] = 9.01984485847386e-08; // Ry_Ri (mM) (in SR_Fluxes)

    Na_Bj[0] = 3.43545459048316; // Na_Bj (mM) (in Na_Buffers)

    Na_Bsl[0] = 0.749601264899653; // Na_Bsl (mM) (in Na_Buffers)
    
    Tn_CL[0] = 0.00893708435270205; // Tn_CL (mM) (in Cytosolic_Ca_Buffers)

    Tn_CHc[0] = 0.117445983314504; // Tn_CHc (mM) (in Cytosolic_Ca_Buffers)

    Tn_CHm[0] = 0.0105996734077994; // Tn_CHm (mM) (in Cytosolic_Ca_Buffers)

    CaM[0] = 0.000295653619580701; // CaM (mM) (in Cytosolic_Ca_Buffers)

    Myo_c[0] = 0.00192645052472679; // Myo_c (mM) (in Cytosolic_Ca_Buffers)

    Myo_m[0] = 0.137557201546068; // Myo_m (mM) (in Cytosolic_Ca_Buffers)

    SRB[0] = 0.00217414510791738; // SRB (mM) (in Cytosolic_Ca_Buffers)
    
    SLL_j[0] = 0.00738583890572642; // SLL_j (mM) (in Junctional_and_SL_Ca_Buffers)
     
    SLL_sl[0] = 0.00988178900584875; // SLL_sl (mM) (in Junctional_and_SL_Ca_Buffers)
       
    SLH_j[0] = 0.0734662466011574; // SLH_j (mM) (in Junctional_and_SL_Ca_Buffers)
       
    SLH_sl[0] = 0.114400081504523; // SLH_sl (mM) (in Junctional_and_SL_Ca_Buffers)
        
    Csqn_b[0] = 1.19772047585784; // Csqn_b (mM) (in SR_Ca_Concentrations)
    
    Ca_sr[0] = 0.555180633859957; // Ca_sr (mM) (in SR_Ca_Concentrations)

    Na_j[0] = 8.3215690202059; // Na_j (mM) (in Na_Concentrations)

    Na_sl[0] = 8.32094589677861; // Na_sl (mM) (in Na_Concentrations)

    Na_i[0] = 8.32114502072456; // Na_i (mM) (in Na_Concentrations)
 
    K_i[0] = 120.0; // K_i (mM) (in K_Concentration)
 
    Ca_j[0] = 0.000175415190830688; // Ca_j (mM) (in Ca_Concentrations)
     
    Ca_sl[0]  = 0.000106544589194246; // Ca_sl (mM) (in Ca_Concentrations)

  double dt = 0.001;
  double stim_end = 100000;// millisecond (in cell)




  for (size_t j = 0; j<t.size(); j++)
   {  
    //    if ((stim_start <= -stim_period * floor(j / stim_period) + j) && (stim_duration + stim_start >= -stim_period * floor(j / stim_period) + j))
    //         {i_Stim = 35*stim_amplitude;}
    //        else 
    //         {i_Stim = 0;}

       if ((j >= stim_start) && (j <= stim_end) && (stim_duration >= -stim_start - stim_period * floor((-stim_start + j) / stim_period) + j))
           i_Stim = 68*stim_amplitude;
         else 
          i_Stim = 0;
       
    //---------------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------------


    dK_i = 0;
    df_Ca_Bj = -0.011900000000000001 * f_Ca_Bj[j] + 1.7 * (1 - f_Ca_Bj[j]) * Ca_j[j];
    df_Ca_Bsl = -0.011900000000000001 * f_Ca_Bsl[j] + 1.7 * (1 - f_Ca_Bsl[j]) * Ca_sl[j];
    fcaCaMSL = 0;
    fcaCaj = 0;
    I_ext = I / (Ampl_gain * Cext);
    RI = 1 - Ry_Rr[j] - Ry_Ro[j] - Ry_Ri[j];
    J_SRleak = (-Ca_j[j] + Ca_sr[j]) * SR_leak_max;
    g_leak = 1 / R_seal;
    dss = 1 / (1 + exp(-0.83333333333333337 - 0.16666666666666666 * V[j]));

    fss = 1 / (1 + exp(3.8888888888888888 + 0.1111111111111111 * V[j])) + 0.59999999999999998 / (1 + exp(2.5 - 0.050000000000000003 * V[j]));
    
    if ((V[j] >= -5.0000005999999999) && (V[j] <= -4.9999994000000001))
        taud = 833333.33332864498 * (5.0000005999999999 + V[j]) * (47619047.618779711 * (1 - exp(1.000000000005626e-7)) / (1 + exp(1.000000000005626e-7)) + 47619047.618779711 * (1 - exp(-1.000000000005626e-7)) / (1 + exp(-1.000000000005626e-7))) - 47619047.618779711 * (1 - exp(1.000000000005626e-7)) / (1 + exp(1.000000000005626e-7));
    else 
        taud = 28.571428571428569 * (1 - exp(-0.83333333333333337 - 0.16666666666666666 * V[j])) / ((1 + exp(-0.83333333333333337 - 0.16666666666666666 * V[j])) * (5 + V[j]));


    
    tauf = 1 / (0.02 + 0.019699999999999999 * exp(-0.23877882250000002 * pow((1 + 0.068965517241379309 * V[j]), 2)));
    
    
    kp_kp = 1 / (1 + exp(7.4880000000000004 - 0.16722408026755853 * V[j]));
    
    rkr = 1 / (1 + exp(3.0833333333333335 + 0.041666666666666664 * V[j]));
    
    tauxr = 230 / (1 + exp(2 + 0.050000000000000003 * V[j])) + 3300 / ((1 + exp(1.2222222222222223 + 0.1111111111111111 * V[j])) * (1 + exp(-2.4444444444444446 - 0.1111111111111111 * V[j])));
    
    xrss = 1 / (1 + exp(-2 - 0.20000000000000001 * V[j]));
    
    
    tauxs = 990.10000000000002 / (1 + exp(-0.17252124645892353 - 0.070821529745042494 * V[j]));
    
    xsss = 1 / (1 + exp(-0.26666666666666666 - 0.070175438596491224 * V[j]));
    
  
    
    if (V[j] >= -40)
        {ah = 0;}
       else 
         {ah = 0.057000000000000002 * exp(-11.764705882352942 + 0.14705882352941177 * shift_INa_inact - 0.14705882352941177 * V[j]);}


    if (V[j] >= -40)
        {aj = 0;}
       else 
         {aj = (37.780000000000001 + V[j]) * (-25428 * exp(0.24440000000000001 * V[j] - 0.24440000000000001 * shift_INa_inact) - 6.9480000000000002e-6 * exp(0.043909999999999998 * shift_INa_inact - 0.043909999999999998 * V[j])) / (1 + exp(24.640530000000002 + 0.311 * V[j] - 0.311 * shift_INa_inact));}


    if (V[j] >= -40)
        {bh = 5.9230769230769234 / (1 + exp(-0.96036036036036043 + 0.0900900900900901 * shift_INa_inact - 0.0900900900900901 * V[j]));}
       else 
         {bh = 310000 * exp(0.34849999999999998 * V[j] - 0.34849999999999998 * shift_INa_inact) + 2.7000000000000002 * exp(0.079000000000000001 * V[j] - 0.079000000000000001 * shift_INa_inact);}


    if (V[j] >= -40)
        {bj = 0.59999999999999998 * exp(0.057000000000000002 * V[j] - 0.057000000000000002 * shift_INa_inact) / (1 + exp(-3.2000000000000002 + 0.10000000000000001 * shift_INa_inact - 0.10000000000000001 * V[j]));}
       else 
         {bj = 0.024240000000000001 * exp(0.01052 * shift_INa_inact - 0.01052 * V[j]) / (1 + exp(-5.5312920000000005 + 0.13780000000000001 * shift_INa_inact - 0.13780000000000001 * V[j]));}


    hss = 0.01 * perc_reduced_inact_for_IpNa + (1 - 0.01 * perc_reduced_inact_for_IpNa) / pow((1 + exp(9.6298788694481825 + 0.13458950201884254 * V[j] - 0.13458950201884254 * shift_INa_inact)), 2);
    
    jss = 0.01 * perc_reduced_inact_for_IpNa + (1 - 0.01 * perc_reduced_inact_for_IpNa) / pow((1 + exp(9.6298788694481825 + 0.13458950201884254 * V[j] - 0.13458950201884254 * shift_INa_inact)), 2);
    
    mss = 1 / pow((1 + exp(-6.2967884828349945 - 0.11074197120708749 * V[j])), 2);
    
    
    taum = 0.12920000000000001 * exp(-8.682389366752302 * pow((1 + 0.021838829438742085 * V[j]), 2)) + 0.064869999999999997 * exp(-0.0089012876052174655 * pow((-1 + 0.2073398299813394 * V[j]), 2));
    
    
    tauxtof = 0.5 + 8.5 * exp(-0.81000000000000005 * pow((1 + 0.022222222222222223 * V[j]), 2));
    
    tauxtos = 0.5 + 9 / (1 + exp(0.20000000000000001 + 0.066666666666666666 * V[j]));
    
    tauytof = 7 + 85 * exp(-7.2727272727272725 * pow((1 + 0.025000000000000001 * V[j]), 2));
    
    tauytos = 30 + 800 / (1 + exp(6 + 0.10000000000000001 * V[j]));
    
    xtoss = 1 / (1 + exp(1.4615384615384615 - 0.076923076923076927 * V[j]));
    

    ytoss = 1 / (1 + exp(3.8999999999999999 + 0.20000000000000001 * V[j]));
    

    i_bck = A0_bck * Scale_bck / ((1 + exp(-V[j] * k_bck)) * Cext);

    i_leak_comp = (-E_l + V[j]) * g_leak * leak_comp_perc / (100 * Cext);

    i_inj[j] = (-i_bck - i_leak_comp + I_ext) * Scaling;

//-----------------------------------------------------
 
//-----------------------------------------------------
    Fsl = 1 - Fjunc;

    GKs_total = (gks_junc * Fjunc + gks_sl * Fsl) * GKs_factor;

    Fsl_CaL = 1 - Fjunc_CaL;
    
    Ka_junc = 1 / (1 + pow(Kdact, 2) / pow(Ca_j[j], 2));
    
    Ka_sl = 1 / (1 + pow(Kdact, 2) / pow(Ca_sl[j], 2));
    
    g_K1 = 0.43033148291193518 * sqrt(Ko) * g_ki_factor;
    
    gkr = 0.43033148291193518 * sqrt(Ko) * gkr_max;
    
    s3_junc = pow(Na_j[j], 3) * Cao + pow(Na_j[j], 3) * KmCao + pow(Nao, 3) * Ca_j[j] + pow(KmNao, 3) * (1 + Ca_j[j] / KmCai) * Ca_j[j] + pow(Nao, 3) * (1 + pow(Na_j[j], 3) / pow(KmNai, 3)) * KmCai;
    
    s3_sl = pow(Na_sl[j], 3) * Cao + pow(Na_sl[j], 3) * KmCao + pow(Nao, 3) * Ca_sl[j] + pow(KmNao, 3) * (1 + Ca_sl[j] / KmCai) * Ca_sl[j] + pow(Nao, 3) * (1 + pow(Na_sl[j], 3) / pow(KmNai, 3)) * KmCai;
    
    sigma = -0.14285714285714285 + 0.14285714285714285 * exp(0.01485884101040119 * Nao);
    
    if (((V[j] >= -9.9999999999999995e-8 * R * Temp / Frdy) && (V[j] <= 9.9999999999999995e-8 * R * Temp / Frdy)) || ((V[j] >= 9.9999999999999995e-8 * R * Temp / Frdy) && (V[j] <= -9.9999999999999995e-8 * R * Temp / Frdy)))
        {ibark = 1.3499999999999999e-14 * (-0.75 * Ko + 0.75 * K_i[j] * exp(9.9999999999999995e-8)) * Frdy / (-1 + exp(9.9999999999999995e-8)) - 5000000 * (-9.9999999999999995e-8 * R * Temp / Frdy + V[j]) * (-1.3499999999999999e-14 * (-0.75 * Ko + 0.75 * K_i[j] * exp(9.9999999999999995e-8)) * Frdy / (-1 + exp(9.9999999999999995e-8)) - 1.3499999999999999e-14 * (-0.75 * Ko + 0.75 * K_i[j] * exp(-9.9999999999999995e-8)) * Frdy / (-1 + exp(-9.9999999999999995e-8))) * Frdy / (R * Temp);}
       else 
        {ibark = 1.35e-7 * pow(Frdy, 2) * (-0.75 * Ko + 0.75 * K_i[j] * exp(V[j] * Frdy / (R * Temp))) * V[j] / ((-1 + exp(V[j] * Frdy / (R * Temp))) * R * Temp);}


    if (((V[j] >= -9.9999999999999995e-8 * R * Temp / Frdy) && (V[j] <= 9.9999999999999995e-8 * R * Temp / Frdy)) || ((V[j] >= 9.9999999999999995e-8 * R * Temp / Frdy) && (V[j] <= -9.9999999999999995e-8 * R * Temp / Frdy)))
        {ibarna_j = 7.4999999999999986e-16 * (-0.75 * Nao + 0.75 * Na_j[j] * exp(9.9999999999999995e-8)) * Frdy / (-1 + exp(9.9999999999999995e-8)) - 5000000 * (-9.9999999999999995e-8 * R * Temp / Frdy + V[j]) * (-7.4999999999999986e-16 * (-0.75 * Nao + 0.75 * Na_j[j] * exp(9.9999999999999995e-8)) * Frdy / (-1 + exp(9.9999999999999995e-8)) - 7.4999999999999986e-16 * (-0.75 * Nao + 0.75 * Na_j[j] * exp(-9.9999999999999995e-8)) * Frdy / (-1 + exp(-9.9999999999999995e-8))) * Frdy / (R * Temp);}
        else 
         {ibarna_j = 7.4999999999999993e-9 * pow(Frdy, 2) * (-0.75 * Nao + 0.75 * Na_j[j] * exp(V[j] * Frdy / (R * Temp))) * V[j] / ((-1 + exp(V[j] * Frdy / (R * Temp))) * R * Temp);}


    if (((V[j] >= -9.9999999999999995e-8 * R * Temp / Frdy) && (V[j] <= 9.9999999999999995e-8 * R * Temp / Frdy)) || ((V[j] >= 9.9999999999999995e-8 * R * Temp / Frdy) && (V[j] <= -9.9999999999999995e-8 * R * Temp / Frdy)))
        {ibarna_sl = 7.4999999999999986e-16 * (-0.75 * Nao + 0.75 * Na_sl[j] * exp(9.9999999999999995e-8)) * Frdy / (-1 + exp(9.9999999999999995e-8)) - 5000000 * (-9.9999999999999995e-8 * R * Temp / Frdy + V[j]) * (-7.4999999999999986e-16 * (-0.75 * Nao + 0.75 * Na_sl[j] * exp(9.9999999999999995e-8)) * Frdy / (-1 + exp(9.9999999999999995e-8)) - 7.4999999999999986e-16 * (-0.75 * Nao + 0.75 * Na_sl[j] * exp(-9.9999999999999995e-8)) * Frdy / (-1 + exp(-9.9999999999999995e-8))) * Frdy / (R * Temp);}
       else 
         {ibarna_sl = 7.4999999999999993e-9 * pow(Frdy, 2) * (-0.75 * Nao + 0.75 * Na_sl[j] * exp(V[j] * Frdy / (R * Temp))) * V[j] / ((-1 + exp(V[j] * Frdy / (R * Temp))) * R * Temp);}


    FoRT = Frdy / (R * Temp);
    
    fnak = 1 / (1 + 0.1245 * exp(-0.10000000000000001 * V[j] * FoRT) + 0.036499999999999998 * sigma * exp(-V[j] * FoRT));
    
    I_nak_junc = fnak * Fjunc * IbarNaK * Ko / ((1 + pow(KmNaip, 4) / pow(Na_j[j], 4)) * (KmKo + Ko));
    
    I_nak_sl = fnak * Fsl * IbarNaK * Ko / ((1 + pow(KmNaip, 4) / pow(Na_sl[j], 4)) * (KmKo + Ko));
    
    I_nak = I_nak_junc + I_nak_sl;
    
    Qpow = -31 + 0.10000000000000001 * Temp;
    
    I_CaK = 0.45000000000000001 * pow(Q10CaL, Qpow) * ((1 - f_Ca_Bj[j] + fcaCaj) * Fjunc_CaL + (1 - f_Ca_Bsl[j] + fcaCaMSL) * Fsl_CaL) * ibark * d[j] * f[j];
    
    I_CaNa_junc = 0.45000000000000001 * pow(Q10CaL, Qpow) * (1 - f_Ca_Bj[j] + fcaCaj) * ibarna_j * d[j] * f[j] * Fjunc_CaL;
    
    I_CaNa_sl = 0.45000000000000001 * pow(Q10CaL, Qpow) * (1 - f_Ca_Bsl[j] + fcaCaMSL) * ibarna_sl * d[j] * f[j] * Fsl_CaL;
    
    I_pca_junc = pow(Ca_j[j], 1.6000000000000001) * pow(Q10SLCaP, Qpow) * Fjunc * IbarSLCaP / (pow(Ca_j[j], 1.6000000000000001) + pow(KmPCa, 1.6000000000000001));
    
    I_pca_sl = pow(Ca_sl[j], 1.6000000000000001) * pow(Q10SLCaP, Qpow) * Fsl * IbarSLCaP / (pow(Ca_sl[j], 1.6000000000000001) + pow(KmPCa, 1.6000000000000001));
    
    Vcell = 1.0000000000000001e-15 * M_PI * pow(cellRadius, 2) * cellLength;
    
    Vjunc = 0.00053900000000000009 * Vcell;
    
    Vmyo = 0.65000000000000002 * Vcell;
    
    dNa_i = (-Na_i[j] + Na_sl[j]) * conc_clamp * J_na_slmyo / Vmyo;
    
    Bmax_SLhighj = 0.000165 * Vmyo / Vjunc;
    
    Bmax_SLlowj = 0.00046000000000000001 * Vmyo / Vjunc;
    
    Vsl = 0.02 * Vcell;
    
    Bmax_SLhighsl = 0.0134 * Vmyo / Vsl;
    
    Bmax_SLlowsl = 0.037400000000000003 * Vmyo / Vsl;
    
    Vsr = 0.035000000000000003 * Vcell;
    
    Bmax_Csqn = 0.14000000000000001 * Vmyo / Vsr;
    
    kCaSR = -(-MinSR + MaxSR) / (1 + pow((ec50SR / Ca_sr[j]), 2.5)) + MaxSR;
    
    eca_junc = 0.5 * log(Cao / Ca_j[j]) / FoRT;
    
    I_cabk_junc = (-eca_junc + V[j]) * Fjunc * GCaB;
    
    eca_sl = 0.5 * log(Cao / Ca_sl[j]) / FoRT;
    
    I_cabk_sl = (-eca_sl + V[j]) * Fsl * GCaB;
    
    ecl = log(Cli / Clo) / FoRT;
    
    I_ClCa_junc = (-ecl + V[j]) * Fjunc * GClCa / (1 + KdClCa / Ca_j[j]);
    
    I_ClCa_sl = (-ecl + V[j]) * Fsl * GClCa / (1 + KdClCa / Ca_sl[j]);
    
    I_ClCa = I_ClCa_junc + I_ClCa_sl;
    
    I_Clbk = (-ecl + V[j]) * GClB;
    
    I_Cl_tot[j] = I_ClCa + I_Clbk;
    
    ek = log(Ko / K_i[j]) / FoRT;
    
    aki = 1.02 / (1 + exp(-14.1227775 + 0.23849999999999999 * V[j] - 0.23849999999999999 * ek));
    
    bki = (0.49124000000000001 * exp(0.43983232 + 0.080320000000000003 * V[j] - 0.080320000000000003 * ek) + exp(-36.698642499999998 + 0.061749999999999999 * V[j] - 0.061749999999999999 * ek)) / (1 + exp(-2.4444678999999998 + 0.51429999999999998 * ek - 0.51429999999999998 * V[j]));
    
    kiss = aki / (aki + bki);
    
    I_ki = (-ek + V[j]) * g_K1 * kiss;
    
    I_kr = (-ek + V[j]) * gkr * rkr * x_kr[j];
    
    ena_junc = log(Nao / Na_j[j]) / FoRT;
    
    I_Na_junc = pow(m[j], 3) * (-ena_junc + V[j]) * h[j] * jj[j] * Fjunc * GNa;
    
    I_nabk_junc = (-ena_junc + V[j]) * Fjunc * GNaB;
    
    ena_sl = log(Nao / Na_sl[j]) / FoRT;
    
    I_Na_sl = pow(m[j], 3) * (-ena_sl + V[j]) * h[j] * jj[j] * Fsl * GNa;
    
    I_nabk_sl = (-ena_sl + V[j]) * Fsl * GNaB;
    
    if (epi == 1)
        {GtoFast = 0.1144 * GtoFast_factor;}
       else 
         {GtoFast = 0.0014039999999999999 * GtoFast_factor;}


    if (epi == 1)
        {GtoSlow = 0.015599999999999999 * GtoSlow_factor;}
       else 
         {GtoSlow = 0.037595999999999997 * GtoSlow_factor;}


    I_tof = (-ek + V[j]) * GtoFast * y_to_f[j] * x_to_f[j] ;
    
    I_tos = (-ek + V[j]) * GtoSlow * x_to_s[j] * y_to_s[j];
    
    I_to = I_tof + I_tos;
    
    I_kp_junc = (-ek + V[j]) * kp_kp * Fjunc * gkp;
    
    I_kp_sl = (-ek + V[j]) * kp_kp * Fsl * gkp;
    
    I_kp = I_kp_junc + I_kp_sl;
    
    J_serca = pow(Q10SRCaP, Qpow) * (pow((Ca_i[j] / Kmf), hillSRCaP) - pow((Ca_sr[j] / Kmr), hillSRCaP)) * Vmax_SRCaP / (1 + pow((Ca_i[j] / Kmf), hillSRCaP) + pow((Ca_sr[j] / Kmr), hillSRCaP));
    
    kiSRCa = kCaSR * kiCa;
    
    koSRCa = koCa / kCaSR;
    
    dRy_Ri = -Ry_Ri[j] * kim - Ry_Ri[j] * kom + pow(Ca_j[j], 2) * RI * koSRCa + kiSRCa * Ry_Ro[j] * Ca_j[j];
    dRy_Ro = Ry_Ri[j] * kim - Ry_Ro[j] * kom + pow(Ca_j[j], 2) * koSRCa * Ry_Rr[j] - kiSRCa * Ry_Ro[j] * Ca_j[j];
    dRy_Rr = RI * kim + Ry_Ro[j] * kom - pow(Ca_j[j], 2) * koSRCa * Ry_Rr[j] - kiSRCa * Ry_Rr[j] * Ca_j[j];
    dCaM = -CaM[j] * koff_cam + (-CaM[j] + Bmax_CaM) * Ca_i[j] * kon_cam;
    dCsqn_b = -Csqn_b[j] * koff_csqn + (-Csqn_b[j] + Bmax_Csqn) * Ca_sr[j] * kon_csqn;
    dMyo_c = -Myo_c[j] * koff_myoca + (-Myo_c[j] - Myo_m[j] + Bmax_myosin) * Ca_i[j] * kon_myoca;
    dMyo_m = -Myo_m[j] * koff_myomg + (-Myo_c[j] - Myo_m[j] + Bmax_myosin) * Mgi * kon_myomg;
    dNa_Bj_dt = -Na_Bj[j] * koff_na + (-Na_Bj[j] + Bmax_Naj) * Na_j[j] * kon_na;
    dNa_Bj = dNa_Bj_dt;
    dNa_Bsl_dt = -Na_Bsl[j] * koff_na + (-Na_Bsl[j] + Bmax_Nasl) * Na_sl[j] * kon_na;
    dNa_Bsl = dNa_Bsl_dt;
    dSLH_j = -SLH_j[j] * koff_slh + (-SLH_j[j] + Bmax_SLhighj) * Ca_j[j] * kon_slh;
    dSLH_sl = -SLH_sl[j] * koff_slh + (-SLH_sl[j] + Bmax_SLhighsl) * Ca_sl[j] * kon_slh;
    dSLL_j = -SLL_j[j] * koff_sll + (-SLL_j[j] + Bmax_SLlowj) * Ca_j[j] * kon_sll;
    dSLL_sl = -SLL_sl[j] * koff_sll + (-SLL_sl[j] + Bmax_SLlowsl) * Ca_sl[j] * kon_sll;
    
    J_CaB_junction = -SLL_j[j] * koff_sll - SLH_j[j] * koff_slh + (-SLL_j[j] + Bmax_SLlowj) * Ca_j[j] * kon_sll + (-SLH_j[j] + Bmax_SLhighj) * Ca_j[j] * kon_slh;
    
    J_CaB_sl = -SLL_sl[j] * koff_sll - SLH_sl[j] * koff_slh + (-SLL_sl[j] + Bmax_SLlowsl) * Ca_sl[j] * kon_sll + (-SLH_sl[j] + Bmax_SLhighsl) * Ca_sl[j] * kon_slh;
    
    dSRB = -SRB[j] * koff_sr + (-SRB[j] + Bmax_SR) * Ca_i[j] * kon_sr;
    dTn_CHc = -Tn_CHc[j] * koff_tnchca + (-Tn_CHc[j] - Tn_CHm[j] + Bmax_TnChigh) * Ca_i[j] * kon_tnchca;
    dTn_CHm = -Tn_CHm[j] * koff_tnchmg + (-Tn_CHc[j] - Tn_CHm[j] + Bmax_TnChigh) * Mgi * kon_tnchmg;
    
    J_CaB_cytosol = -Tn_CL[j] * koff_tncl - Tn_CHc[j] * koff_tnchca - Tn_CHm[j] * koff_tnchmg - CaM[j] * koff_cam - Myo_c[j] * koff_myoca - Myo_m[j] * koff_myomg - SRB[j] * koff_sr + (-Tn_CL[j] + Bmax_TnClow) * Ca_i[j] * kon_tncl + (-CaM[j] + Bmax_CaM) * Ca_i[j] * kon_cam + (-SRB[j] + Bmax_SR) * Ca_i[j] * kon_sr + (-Tn_CHc[j] - Tn_CHm[j] + Bmax_TnChigh) * Ca_i[j] * kon_tnchca + (-Tn_CHc[j] - Tn_CHm[j] + Bmax_TnChigh) * Mgi * kon_tnchmg + (-Myo_c[j] - Myo_m[j] + Bmax_myosin) * Ca_i[j] * kon_myoca + (-Myo_c[j] - Myo_m[j] + Bmax_myosin) * Mgi * kon_myomg;
    
    dCa_i = -J_CaB_cytosol + (-Ca_i[j] + Ca_sl[j]) * J_ca_slmyo / Vmyo - J_serca * Vsr / Vmyo;
    dTn_CL = -Tn_CL[j] * koff_tncl + (-Tn_CL[j] + Bmax_TnClow) * Ca_i[j] * kon_tncl;
    
    J_SRCarel = (-Ca_j[j] + Ca_sr[j]) * Ry_Ro[j] * ks;
   
    dCa_sr = -J_SRCarel + Csqn_b[j] * koff_csqn - J_SRleak * Vmyo / Vsr - (-Csqn_b[j] + Bmax_Csqn) * Ca_sr[j] * kon_csqn + J_serca;
    
    s1_junc = pow(Na_j[j], 3) * Cao * exp(V[j] * FoRT * nu);
    
    s1_sl = pow(Na_sl[j], 3) * Cao * exp(V[j] * FoRT * nu);
    
    s2_junc = pow(Nao, 3) * Ca_j[j] * exp((-1 + nu) * V[j] * FoRT);
    
    I_ncx_junc = pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc * Fjunc * IbarNCX / ((1 + ksat * exp((-1 + nu) * V[j] * FoRT)) * s3_junc);
    
    s2_sl = pow(Nao, 3) * Ca_sl[j] * exp((-1 + nu) * V[j] * FoRT);
    
    I_ncx_sl = pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Ka_sl * Fsl * IbarNCX / ((1 + ksat * exp((-1 + nu) * V[j] * FoRT)) * s3_sl);
    
    I_Na_tot_junc = 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc + I_Na_junc + I_nabk_junc;
    
    dNa_j = (-dNa_Bj_dt + (-Na_j[j] + Na_sl[j]) * J_na_juncsl / Vjunc - I_Na_tot_junc * Cmem / (Frdy * Vjunc)) * conc_clamp;
    
    I_Na_tot_sl = 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl;
    
    dNa_sl = (-dNa_Bsl_dt + (-Na_sl[j] + Na_j[j]) * J_na_juncsl / Vsl + (-Na_sl[j] + Na_i[j]) * J_na_slmyo / Vsl - I_Na_tot_sl * Cmem / (Frdy * Vsl)) * conc_clamp;
    
    I_Na_tot[j] = I_Na_tot_junc + I_Na_tot_sl;
    
    if (((V[j] >= -4.9999999999999998e-8 * R * Temp / Frdy) && (V[j] <= 4.9999999999999998e-8 * R * Temp / Frdy)) || ((V[j] >= 4.9999999999999998e-8 * R * Temp / Frdy) && (V[j] <= -4.9999999999999998e-8 * R * Temp / Frdy)))
        {ibarca_j = 1.08e-10 * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_j[j] * exp(9.9999999999999995e-8)) * Frdy * pCa_max / (-1 + exp(9.9999999999999995e-8)) - 10000000 * (-4.9999999999999998e-8 * R * Temp / Frdy + V[j]) * (-1.08e-10 * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_j[j] * exp(9.9999999999999995e-8)) * Frdy * pCa_max / (-1 + exp(9.9999999999999995e-8)) - 1.08e-10 * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_j[j] * exp(-9.9999999999999995e-8)) * Frdy * pCa_max / (-1 + exp(-9.9999999999999995e-8))) * Frdy / (R * Temp);}
       else 
         {ibarca_j = 0.00216 * pow(Frdy, 2) * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_j[j] * exp(2 * V[j] * Frdy / (R * Temp))) * V[j] * pCa_max / ((-1 + exp(2 * V[j] * Frdy / (R * Temp))) * R * Temp);}


    I_Ca_junc = 0.45000000000000001 * pow(Q10CaL, Qpow) * (1 - f_Ca_Bj[j] + fcaCaj) * ibarca_j * d[j] * f[j] * Fjunc_CaL;
    
    I_Ca_tot_junc = -2 * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc;
    
    dCa_j = -J_CaB_junction + (-Ca_j[j] + Ca_sl[j]) * J_ca_juncsl / Vjunc + J_SRCarel * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc - 0.5 * I_Ca_tot_junc * Cmem / (Frdy * Vjunc);
    
    if (((V[j] >= -4.9999999999999998e-8 * R * Temp / Frdy) && (V[j] <= 4.9999999999999998e-8 * R * Temp / Frdy)) || ((V[j] >= 4.9999999999999998e-8 * R * Temp / Frdy) && (V[j] <= -4.9999999999999998e-8 * R * Temp / Frdy)))
        {ibarca_sl = 1.08e-10 * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_sl[j] * exp(9.9999999999999995e-8)) * Frdy * pCa_max / (-1 + exp(9.9999999999999995e-8)) - 10000000 * (-4.9999999999999998e-8 * R * Temp / Frdy + V[j]) * (-1.08e-10 * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_sl[j] * exp(9.9999999999999995e-8)) * Frdy * pCa_max / (-1 + exp(9.9999999999999995e-8)) - 1.08e-10 * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_sl[j] * exp(-9.9999999999999995e-8)) * Frdy * pCa_max / (-1 + exp(-9.9999999999999995e-8))) * Frdy / (R * Temp);}
       else 
         {ibarca_sl = 0.00216 * pow(Frdy, 2) * (-0.34100000000000003 * Cao + 0.34100000000000003 * Ca_sl[j] * exp(2 * V[j] * Frdy / (R * Temp))) * V[j] * pCa_max / ((-1 + exp(2 * V[j] * Frdy / (R * Temp))) * R * Temp);}


    I_Ca_sl = 0.45000000000000001 * pow(Q10CaL, Qpow) * (1 - f_Ca_Bsl[j] + fcaCaMSL) * ibarca_sl * d[j] * f[j] * Fsl_CaL;
    
    I_Ca_tot_sl = -2 * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl;
    
    dCa_sl = -J_CaB_sl + (-Ca_sl[j] + Ca_i[j]) * J_ca_slmyo / Vsl + (-Ca_sl[j] + Ca_j[j]) * J_ca_juncsl / Vsl - 0.5 * I_Ca_tot_sl * Cmem / (Frdy * Vsl);
    
    I_Ca_tot[j] = I_Ca_tot_junc + I_Ca_tot_sl;
    
    eks = log((Nao * pNaK + Ko) / (Na_i[j] * pNaK + K_i[j])) / FoRT;
    
    I_ks = pow(x_ks[j], 2) * (-eks + V[j]) * GKs_total;
    
    I_K_tot[j] = -2 * I_nak + I_CaK + I_ki + I_kp + I_kr + I_ks + I_to;
    
    I_tot[j] = I_K_tot[j] + I_Ca_tot[j] + I_Cl_tot[j] + I_Na_tot[j];
    
    dV = -15*i_inj[j] - 15*I_tot[j] - 8*i_Stim;

   


    //------------------------------------------------------------------------------
    // Integration & Output
    //------------------------------------------------------------------------------
    // Rush-Larsen method
  
    m[j+1] = mss + (m[j] - mss)*exp(-dt/taum);

    h[j+1] = hss + (h[j] - hss)*exp(-dt/tauh);

    jj[j+1] = jss + (jj[j] - jss)*exp(-dt/tauj);

    x_kr[j+1] = xrss + (x_kr[j] - xrss)*exp(-dt/tauxr);

    x_ks[j+1] = xsss + (x_ks[j] - xsss)*exp(-dt/tauxs);

    x_to_s[j+1] = xtoss + (x_to_s[j] - xtoss)*exp(-dt/tauxtos);

    y_to_s[j+1] = ytoss + (y_to_s[j] - ytoss)*exp(-dt/tauytos);

    x_to_f[j+1] = xtoss + (x_to_f[j]  - xtoss)*exp(-dt/tauxtof);

    y_to_f[j+1] = ytoss + (y_to_f[j] - ytoss)*exp(-dt/tauytof);

    d[j+1]  = dss + (d[j] - dss)*exp(-dt/taud);

    f[j+1] = fss + (f[j] - fss)*exp(-dt/tauf);


// Remainder: Forward Euler

 
    V[j+1] = V[j] + dt * dV;

  
    Ca_i[j+1] = Ca_i[j] + dt * dCa_i;


    f_Ca_Bj[j+1] = f_Ca_Bj[j] + dt * df_Ca_Bj;


    f_Ca_Bsl[j+1] = f_Ca_Bsl[j] + dt * df_Ca_Bsl;


    Ry_Rr[j+1] = Ry_Rr[j] + dt * dRy_Rr;


    Ry_Ro[j+1] = Ry_Ro[j] + dt * dRy_Ro;


    Ry_Ri[j+1] = Ry_Ri[j] + dt * dRy_Ri;


    Na_Bj[j+1] = Na_Bj[j] + dt * dNa_Bj;


    Na_Bsl[j+1] = Na_Bsl[j] + dt * dNa_Bsl;


    Tn_CL[j+1] = Tn_CL[j] + dt * dTn_CL;


    Tn_CHc[j+1] = Tn_CHc[j] + dt * dTn_CHc;


    Tn_CHm[j+1] = Tn_CHm[j] + dt * dTn_CHm;


    CaM[j+1] = CaM[j] + dt * dCaM;


    Myo_c[j+1] = Myo_c[j] + dt * dMyo_c;


    Myo_m[j+1] = Myo_m[j] + dt * dMyo_m;

 
    SRB[j+1] = SRB[j] + dt * dSRB;

 
    SLL_j[j+1] = SLL_j[j] + dt * dSLL_j;


    SLL_sl[j+1] = SLL_sl[j] + dt * dSLL_sl;

 
    SLH_j[j+1] = SLH_j[j] + dt * dSLH_j;


    SLH_sl[j+1] = SLH_sl[j] + dt * dSLH_sl;


    Csqn_b[j+1] = Csqn_b[j] + dt * dCsqn_b;


    Ca_sr[j+1] = Ca_sr[j] + dt * dCa_sr;


    Na_j[j+1] = Na_j[j] + dt * dNa_j;


    Na_sl[j+1] = Na_sl[j] + dt * dNa_sl;


    Na_i[j+1] = Na_i[j] + dt * dNa_i;


    K_i[j+1] = K_i[j] + dt * dK_i;


    Ca_j[j+1] = Ca_j[j] + dt * dCa_j;


    Ca_sl[j+1] = Ca_sl[j] + dt * dCa_sl;
   } 
plt::figure();
plt::plot(V);
     
// plt::figure();
// plt::plot(i_inj); 

// plt::figure();
// plt::plot(I_tot); 

// plt::figure();
// plt::plot(I_K_tot);

// plt::figure();
// plt::plot(I_Ca_tot);

// plt::figure();
// plt::plot(I_Cl_tot);

// plt::figure();
// plt::plot(I_Na_tot);


plt::show();
}

 






//==============================================================================
// End of file
//====================
