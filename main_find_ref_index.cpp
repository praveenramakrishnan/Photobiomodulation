/*This program calculates the amplitudes of the reflected and transmitted wave when 
a multi layered dielectric is illuminated by an obliquely incident transverse electric source.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>

using namespace std;
typedef complex<double> dcomplex;
#define j dcomplex(0.0,1.0)
const double pi = 3.1415926535897932;
const double C0 = 299792458.0;//speed of light in vacuum
const double MU0 = 4.0*pi*1.0e-7; 
const double EPS0 = 1.0/(MU0*C0*C0);
const double planckconstantbar = 1.054571800e-34;

void matrixMultiply(dcomplex M1[][2], dcomplex M2[][2], dcomplex M_ans[0][2])
{
    //This function multiplies matrices M1 and M2 and store the answer in M_ans. ie M_ans = M1*M2
    int iloop1,iloop2,iloop3;
    dcomplex temp[2][2]; //matrix to store temporary values

   for(iloop1=0; iloop1<2; ++iloop1)
       for(iloop2=0; iloop2<2; ++iloop2)
       {
           temp[iloop1][iloop2] = 0.0;
       }

  for(iloop1=0; iloop1<2; ++iloop1)
     for(iloop2=0; iloop2<2; ++iloop2)
        for(iloop3=0; iloop3<2; ++iloop3)
           {
              temp[iloop1][iloop2]+=M1[iloop1][iloop3]*M2[iloop3][iloop2];
           }
   for(iloop1=0; iloop1<2; ++iloop1)
      for(iloop2=0; iloop2<2; ++iloop2)
        {
           M_ans[iloop1][iloop2] =temp[iloop1][iloop2];
        }
}

void calculations_TE_wave(int ni, double inputpowerdensity, vector<dcomplex>& vMu, vector<dcomplex>& vEpsilon, vector<dcomplex>& vKprop_z, vector<double>& vThickness, double&  active_power_reflected_TE, double& active_power_transmitted_TE );

void calculations_TM_wave(int ni, double inputpowerdensity, vector<dcomplex>& vMu, vector<dcomplex>& vEpsilon, vector<dcomplex>& vKprop_z, vector<double>& vThickness, double&  active_power_reflected_TM, double& active_power_transmitted_TM );

int main()
{
	
   //Input data
   int ni;
   cin >> ni;

   vector <double> vThickness(ni-1); //vector containing thickness of layers
   vector <double> vZcoord(ni); //vector containing z-coordinate of intefaces
   vector <dcomplex> vEpsilon(ni+1); //vector containing permittivity of layers
   vector <dcomplex> vMu(ni+1); //vector containing permeability of layers
   vector <dcomplex> vImp(ni+1); //vector containing impedence of layers
   vector <dcomplex> vKprop_z(ni+1); //vector containg z-component of propagation constant of layers

   // Input the thickness of layers
   for(int iloop = 0; iloop < ni-1;  iloop++)
      {
         //cout << "input the thickness of layer #" << iloop+2 <<  endl;
	 cin >> vThickness[iloop];
      }

   //Input wavelength in air of incident wave
   //cout << "input the wavelength in air of incident wave" << endl; 
   double lambda_input, lambda0;
   cin >> lambda_input;

   //Input power density of incident wave
   //cout << "input power intensity in W/m^2 of incident wave" << endl;
   double inputpowerdensity;
   cin >> inputpowerdensity;

   //Input the angle of incidence of the source
   //cout << "input the angle of incident wave in radians" << endl;
   double theta_in;
   cin >> theta_in;


   //Read the values of refractive index of layers

    dcomplex refractiveindex, Kprop, Kprop_x_1;

    double realpartofrefindex[ni], imagpartofrefindex, extinctioncoefficient[ni];

   for(int iloop = 0; iloop < ni+1; iloop++)
   {
	   //cout << "Input real part of refractive index of layer " << iloop+1 << endl;
	   cin >> realpartofrefindex[iloop];
	   //cout << "Extinction coefficient of layer " << iloop+1 << endl;
	   cin >> extinctioncoefficient[iloop];
   }
    //Read measured values
    double Power_out_measured, Power_reflected_measured;
    cin >> Power_out_measured;
    cin >> Power_reflected_measured;

    //Read output file name
    string outfilename;
    cin >>outfilename; 


   //Calculations

   double omega, photonenergy;
  
   double factor_of_polarization = 0.5;  //The relative strength of TE and TM polarization

   double active_power_reflected_TE, active_power_reflected_TM, active_power_transmitted_TE, active_power_transmitted_TM;
   
   double active_power_reflected_avg=0.0, active_power_transmitted_avg=0.0;

   int n_loop=0, index_unknown_layer=1;
 
   double real_n, k_ext;
   double active_power_reflected,active_power_transmitted, reflection_coefficient,relative_error ; 
   double thickness_layer;   
   for(real_n=1.0001; real_n<=2.0; real_n+=1.0e-4)
   {
   for(k_ext=0.0; k_ext<=100.0; k_ext+=0.1)
   {

     n_loop = 0;
     active_power_reflected_avg = 0;
     active_power_transmitted_avg = 0;

   for(thickness_layer=0.88e-3; thickness_layer<=1.06e-3; thickness_layer+=1.0e-6)

   {

   for(lambda0=lambda_input-5.0e-9; lambda0 <= lambda_input+5e-9; lambda0+=1e-9)
   {
 
   omega = (2*pi*C0) / lambda0; //angular frequency of incident wave
   photonenergy = planckconstantbar*omega; //energy of one photon of incident wave
   
   realpartofrefindex[index_unknown_layer] = real_n;
   extinctioncoefficient[index_unknown_layer] = k_ext;
   vThickness[index_unknown_layer-1]=thickness_layer;

   for(int iloop = 0; iloop < ni+1; iloop++)
   {
	   imagpartofrefindex = -extinctioncoefficient[iloop]*lambda0/(4.0*pi);
	   refractiveindex = realpartofrefindex[iloop] + (1.0*j)*imagpartofrefindex;
	   vEpsilon[iloop] = refractiveindex*refractiveindex*EPS0;
	   vMu[iloop] = MU0;
	   Kprop = omega*sqrt(vEpsilon[iloop]*vMu[iloop]);
	   Kprop_x_1 = omega*sqrt(vEpsilon[0]*vMu[0])*sin(theta_in);
	   vKprop_z[iloop] = sqrt(Kprop* Kprop - Kprop_x_1*Kprop_x_1);
	   vImp[iloop] = sqrt(vMu[iloop]/vEpsilon[iloop]);  
   }
     

   //The z-coordinates of interfaces 
   
   vZcoord[0] = 0;
   
   for(int iloop = 1; iloop < ni; iloop++)
   {
       vZcoord[iloop] = vZcoord[iloop-1] + vThickness[iloop-1];
   }


  
   //Calculations for TE wave
	calculations_TE_wave(ni, inputpowerdensity*factor_of_polarization, vMu, vEpsilon, vKprop_z, vThickness, active_power_reflected_TE, active_power_transmitted_TE );

   //Calculations for TM wave
	calculations_TM_wave(ni, inputpowerdensity*(1.0-factor_of_polarization), vMu, vEpsilon, vKprop_z, vThickness, active_power_reflected_TM, active_power_transmitted_TM );

   //Calculate average
   active_power_reflected_avg += (active_power_reflected_TE + active_power_reflected_TM);
   active_power_transmitted_avg += (active_power_transmitted_TE + active_power_transmitted_TM);
   n_loop+=1;
   }
   }

   //Print output
   active_power_reflected = active_power_reflected_avg/n_loop;
   active_power_transmitted = active_power_transmitted_avg/n_loop;
   reflection_coefficient = active_power_reflected/inputpowerdensity;
   relative_error = (active_power_reflected - Power_reflected_measured)/active_power_reflected;
   cout << realpartofrefindex[index_unknown_layer] << "\t" << extinctioncoefficient[index_unknown_layer] << "\t" << vThickness[0] << "\t"<< lambda_input << "\t" <<  reflection_coefficient << "\t" << relative_error << endl; 
   }
   }
   return 0;
}

void calculations_TE_wave(int ni, double inputpowerdensity, vector<dcomplex>& vMu, vector<dcomplex>& vEpsilon, vector<dcomplex>& vKprop_z, vector<double>& vThickness,  double&  active_power_reflected_TE, double& active_power_transmitted_TE )
{

 //Calculate amplitude of incident wave a1
   dcomplex a1 = sqrt(2* inputpowerdensity *sqrt(vMu[0] / vEpsilon[0]));
   
   //Calculate the matrix M
   
   dcomplex Ti[2][2][ni-1]; //Ti matrix encodes effect of transmission through layer i+1. 

   for(int iloop=0; iloop < ni-1; iloop++)
   {
      Ti[0][0][iloop] = cos(vKprop_z[iloop+1]*vThickness[iloop])-(1.0*j)*sin(vKprop_z[iloop+1]*vThickness[iloop]);
      Ti[0][1][iloop] = 0.0;
      Ti[1][0][iloop] = 0.0;
      Ti[1][1][iloop] = cos(vKprop_z[iloop+1]*vThickness[iloop])+(1.0*j)*sin(vKprop_z[iloop+1]*vThickness[iloop]);
   }

   dcomplex  Mi[2][2][ni]; //Mi matrix encodes effect of interface i. 

   for(int iloop=0; iloop < ni; iloop++)
   {
      Mi[0][0][iloop] = 0.5*(1.0+(vMu[iloop+1]/vMu[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
      Mi[0][1][iloop] = 0.5*(1.0-(vMu[iloop+1]/vMu[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
      Mi[1][0][iloop] = 0.5*(1.0-(vMu[iloop+1]/vMu[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
      Mi[1][1][iloop] = 0.5*(1.0+(vMu[iloop+1]/vMu[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
   }

   dcomplex Mtot[2][2]; //array containg the total transmission matrix
   Mtot[0][0] = Mi[0][0][0];
   Mtot[0][1] = Mi[0][1][0];
   Mtot[1][0] = Mi[1][0][0];
   Mtot[1][1] = Mi[1][1][0];

   dcomplex Ttemp[2][2], Mtemp[2][2]; //matrices storing temporary values

   for(int iloop = 0; iloop < ni-1; iloop++)
   {
     Ttemp[0][0] = Ti[0][0][iloop];
     Ttemp[0][1] = Ti[0][1][iloop];
     Ttemp[1][0] = Ti[1][0][iloop];
     Ttemp[1][1] = Ti[1][1][iloop];

     Mtemp[0][0] = Mi[0][0][iloop+1];
     Mtemp[0][1] = Mi[0][1][iloop+1];
     Mtemp[1][0] = Mi[1][0][iloop+1];
     Mtemp[1][1] = Mi[1][1][iloop+1];

     matrixMultiply(Ttemp,Mtot,Mtot); 
     matrixMultiply(Mtemp,Mtot,Mtot); 

   }

   dcomplex b1 = - Mtot[1][0]*a1/(Mtot[1][1]);
   dcomplex anplus1 = (Mtot[0][0]*Mtot[1][1]-Mtot[0][1]*Mtot[1][0])/Mtot[1][1]*a1;

   active_power_reflected_TE = real(0.5*b1*conj(b1*sqrt(vEpsilon[0]/vMu[0])));
   active_power_transmitted_TE = real(0.5*anplus1*conj(anplus1*sqrt(vEpsilon[ni]/vMu[ni])));

}

void calculations_TM_wave(int ni, double inputpowerdensity, vector<dcomplex>& vMu, vector<dcomplex>& vEpsilon, vector<dcomplex>& vKprop_z, vector<double>& vThickness,  double&  active_power_reflected_TM, double& active_power_transmitted_TM )
{
 //Calculate amplitude of incident wave a1
   dcomplex a1 = sqrt(2* inputpowerdensity *sqrt(vEpsilon[0] / vMu[0]));
   
   //Calculate the matrix M
   
   dcomplex Ti[2][2][ni-1]; //Ti matrix encodes effect of transmission through layer i+1. 

   for(int iloop=0; iloop < ni-1; iloop++)
   {
      Ti[0][0][iloop] = cos(vKprop_z[iloop+1]*vThickness[iloop])-(1.0*j)*sin(vKprop_z[iloop+1]*vThickness[iloop]);
      Ti[0][1][iloop] = 0.0;
      Ti[1][0][iloop] = 0.0;
      Ti[1][1][iloop] = cos(vKprop_z[iloop+1]*vThickness[iloop])+(1.0*j)*sin(vKprop_z[iloop+1]*vThickness[iloop]);
   }

   dcomplex  Mi[2][2][ni]; //Mi matrix encodes effect of interface i. 

   for(int iloop=0; iloop < ni; iloop++)
   {
      Mi[0][0][iloop] = 0.5*(1.0+(vEpsilon[iloop+1]/vEpsilon[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
      Mi[0][1][iloop] = 0.5*(1.0-(vEpsilon[iloop+1]/vEpsilon[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
      Mi[1][0][iloop] = 0.5*(1.0-(vEpsilon[iloop+1]/vEpsilon[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
      Mi[1][1][iloop] = 0.5*(1.0+(vEpsilon[iloop+1]/vEpsilon[iloop])*(vKprop_z[iloop]/vKprop_z[iloop+1]));
   }

   dcomplex Mtot[2][2]; //array containg the total transmission matrix
   Mtot[0][0] = Mi[0][0][0];
   Mtot[0][1] = Mi[0][1][0];
   Mtot[1][0] = Mi[1][0][0];
   Mtot[1][1] = Mi[1][1][0];

   dcomplex Ttemp[2][2], Mtemp[2][2]; //matrices storing temporary values

   for(int iloop = 0; iloop < ni-1; iloop++)
   {
     Ttemp[0][0] = Ti[0][0][iloop];
     Ttemp[0][1] = Ti[0][1][iloop];
     Ttemp[1][0] = Ti[1][0][iloop];
     Ttemp[1][1] = Ti[1][1][iloop];

     Mtemp[0][0] = Mi[0][0][iloop+1];
     Mtemp[0][1] = Mi[0][1][iloop+1];
     Mtemp[1][0] = Mi[1][0][iloop+1];
     Mtemp[1][1] = Mi[1][1][iloop+1];

     matrixMultiply(Ttemp,Mtot,Mtot); 
     matrixMultiply(Mtemp,Mtot,Mtot); 

   }

   dcomplex b1 = - Mtot[1][0]*a1/(Mtot[1][1]);
   dcomplex anplus1 = (Mtot[0][0]*Mtot[1][1]-Mtot[0][1]*Mtot[1][0])/Mtot[1][1]*a1;

   dcomplex A1 = a1*sqrt(vMu[0]/vEpsilon[0]);
   dcomplex B1 = b1*sqrt(vMu[0]/vEpsilon[0]);
   dcomplex Anplus1 = anplus1*sqrt(vMu[ni]/vEpsilon[ni]);

   active_power_reflected_TM = real(0.5*B1*conj(B1*sqrt(vEpsilon[0]/vMu[0])));
   active_power_transmitted_TM = real(0.5*Anplus1*conj(Anplus1*sqrt(vEpsilon[ni]/vMu[ni])));

}
