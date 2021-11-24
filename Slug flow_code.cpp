#include<stdio.h>
#include <sstream>
#include <fstream>
#include<math.h>
#include<dir.h>
using namespace std;

#define NI 5
#define NJ 5
#define Tmax .4
#define L 1
#define H 1
#define alpha 0.01
#define q_gen 0
#define Rho 1
#define Cp 1
#define DELT 0.1
#define pi 3.14

double PHI[NJ][NI];
double PHI_OLD[NJ][NI];
double PHI_previ[NJ][NI];
double PHI_new[NI-1];
double X0[NJ][NI];
double Y0[NJ][NI];
double X1[NJ][NI];
double Y1[NJ][NI];
double X2[NJ][NI];
double Y2[NJ][NI];
double X3[NJ][NI];
double Y3[NJ][NI];
double XCELL_X[NJ][NI];
double XCELL_Y[NJ][NI];
double S0[NJ][NI];
double S1[NJ][NI];
double S2[NJ][NI];
double S3[NJ][NI];
double DXF0[NJ][NI];
double DXF1[NJ][NI];
double DXF2[NJ][NI];
double DXF3[NJ][NI];
double Vol[NJ][NI];
double aw[NJ][NI];
double as[NJ][NI];
double ae[NJ][NI];
double an[NJ][NI];
double ap[NJ][NI];
double a[NI];
double b[NI];
double c[NI];
double d[NI];
double Sp[NJ][NI];
double DELX=L/(NI-2.0);
double DELY=H/(NJ-2.0);
int i;
int j;
int k;
int NCELLI = NI-1;
int NCELLJ = NJ-1;
int NSTEP=1;
int kmax=NCELLI;
double S=(q_gen)/(Rho*Cp);
double vel=1.0;
double theta=pi/4;
double u=-1;
double v=1;
//Rho*vel*cos(theta)
// Rho*vel*sin(theta)
void grid_generation(){
	for(j=1;j<NCELLJ;j++){
		for(i=1;i<NCELLI;i++){
			X0[j][i]=DELX*(i-1);
			Y0[j][i]=DELY*(j-1);
			X1[j][i]=DELX*i;
			Y1[j][i]=DELY*(j-1);
			X2[j][i]=DELX*i;
			Y2[j][i]=DELY*j;
			X3[j][i]=DELX*(i-1);
			Y3[j][i]=DELY*j;

			XCELL_X[j][i]=0.25*(X0[j][i]+X1[j][i]+X2[j][i]+X3[j][i]);
			XCELL_Y[j][i]=0.25*(Y0[j][i]+Y1[j][i]+Y2[j][i]+Y3[j][i]);
		}
	}
	for(j=1;j<NCELLJ;j++){
		XCELL_X[j][0]=0;
		XCELL_Y[j][0]=XCELL_Y[j][1];

		XCELL_X[j][NCELLI]=L;
		XCELL_Y[j][NCELLI]=XCELL_Y[j][NCELLI-1];

	}
	for(i=1;i<NCELLI;i++){
		XCELL_X[0][i]=XCELL_X[1][i];
		XCELL_Y[0][i]=0;

		XCELL_X[NCELLJ][i]=XCELL_X[NCELLJ-1][i];
		XCELL_Y[NCELLJ][i]=H;
	}
	for(j=1;j<NCELLJ;j++){
		for(i=1;i<NCELLI;i++){
			S0[j][i]=-DELY;
			S1[j][i]=-DELX;
			S2[j][i]=DELY;
			S3[j][i]=DELX;

			DXF0[j][i]=DELX;
			DXF1[j][i]=DELY;
			DXF2[j][i]=DELX;
			DXF3[j][i]=DELY;

			Vol[j][i]=DELX*DELY;
		}
	}
	XCELL_X[NCELLJ][0] = 0;
	XCELL_Y[NCELLJ][0] = H;
	XCELL_X[0][NCELLI] = L;
	XCELL_Y[0][NCELLI] = 0;
	XCELL_X[NCELLJ][NCELLI] = L;
	XCELL_Y[NCELLJ][NCELLJ] = H;
}
void initial_condition(){
	for(j=0;j<=NCELLJ;j++){
		for(i=0;i<=NCELLI;i++){
			PHI_OLD[j][i]=0;
			PHI[j][i]=0;
		}
	}
}
void boundary_condition(){
	for(j=1;j<NCELLJ;j++){
		PHI[j][0]=PHI[j][1] + 3.33;
		PHI_OLD[j][0]=PHI_OLD[j][1] + 3.33;
		PHI[j][NCELLI]=300;
		PHI_OLD[j][NCELLI]=300;
	}
	for(i=1;i<NCELLI;i++){
		PHI[0][i]=100;
		PHI_OLD[0][i]=100;
		PHI[NCELLJ][i]= PHI[NCELLJ-1][i];
		PHI_OLD[NCELLJ][i]=PHI_OLD[NCELLJ-1][i];
	}
}
void coefficient(){
	for(j=1;j<NCELLJ;j++){
		for(i=1;i<NCELLI;i++){
			aw[j][i]=fabs(alpha*S0[j][i]/DXF0[j][i]);
			as[j][i]=fabs(alpha*S1[j][i]/DXF1[j][i]);
			ae[j][i]=fabs(alpha*S2[j][i]/DXF2[j][i]);
			an[j][i]=fabs(alpha*S3[j][i]/DXF3[j][i]);
			ap[j][i]=(aw[j][i]+as[j][i]+ae[j][i]+an[j][i]);
		}
	}
}
void solve_expl(){
	for(j=1;j<NCELLJ;j++){  //Advection(FOU) and diffusion PHI_OLD[j][i]+
		for(i=1;i<NCELLI;i++){
			PHI[j][i] = PHI_OLD[j][i] + (DELT/Vol[j][i])*(aw[j][i]*PHI_OLD[j][i-1] + ae[j][i]*PHI_OLD[j][i+1] + as[j][i]*PHI_OLD[j-1][i] + an[j][i]*PHI_OLD[j+1][i] - ap[j][i]*PHI_OLD[j][i])+ (0.125*u*aw[j][i]*((3*PHI_OLD[j][i-1])+6*PHI_OLD[j][i]-PHI_OLD[j][i+1]))-(u*ae[j][i]*(PHI_OLD[j][i+1]))+(v*as[j][i]*((PHI_OLD[j-1][i])))-(0.125*v*an[j][i]*((3*PHI_OLD[j+1][i])+6*PHI_OLD[j][i]-PHI_OLD[j-1][i])) ;
			//+(u*aw[j][i]*PHI_OLD[j][i-1])-(u*ae[j][i]*PHI_OLD[j][i]) (v*as[j][i]*PHI_OLD[j-1][i])-(v*an[j][i]*PHI_OLD[j][i])+((q_gen*DELT)/(Rho*Cp))
		}
	}
}

void update(){
	for(j=1;j<NCELLJ;j++){
		for(i=1;i<NCELLI;i++){
			PHI_OLD[j][i]=PHI[j][i];
		}
	}
}
void write_file()
{
	int i,j;
	ofstream file;
	stringstream ss;
	ss << "berger/TEMP_PROFILE_TIME_" << NSTEP*DELT << ".dat";
	file.open(ss.str().c_str());
	file<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<", "<<'"'<<"TEMPERATURE"<<'"'<<endl;
	file<<"ZONE I="<<NI<<", J="<<NJ<<", F=POINT"<< endl;
	file<<"SOLUTIONTIME="<<NSTEP*DELT<< endl;
	for(j=0;j<=NCELLJ;j++)
	{
		for(i=0;i<=NCELLI;i++)
		{
			file<<XCELL_X[j][i]<<" "<<XCELL_Y[j][i]<<" "<<PHI[j][i]<<endl;
		}
	}
}

int main(){
	mkdir("berger");
	grid_generation();
	initial_condition();
	coefficient();
	for(NSTEP=1;NSTEP<=(Tmax/DELT);NSTEP++){
		boundary_condition();
		solve_expl();
		//if(NSTEP%200==0){
			write_file();
		//}//
		update();
	}
	return 0;
}
