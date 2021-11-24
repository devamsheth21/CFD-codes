#include<stdio.h>
#include <sstream>
#include <fstream>
#include<math.h>
#include<dir.h>
#include<math.h>

using namespace std;

#define NI 102
#define NJ 102
#define Tmax 10
#define L 1
#define H 1
#define alpha 0.01
#define q_gen 0
#define Rho 1
#define Cp 1
#define DELT 0.0001
#define pi 3.14

double PHI[NJ][NI];
double PHI_OLD[NJ][NI];
double PHI_previ[NJ][NI];
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
double DELX=L/(NI-2.0);
double DELY=H/(NJ-2.0);
int i;
int j;
int NCELLI = NI-1;
int NCELLJ = NJ-1;
int NSTEP=1;
double vel=1.0;
double theta=pi/4;
double u=vel*cos(theta);
double v=vel*sin(theta);
int count=0;

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
			PHI_previ[j][i]=0;
		}
	}
}
void boundary_condition(){
	for(j=1;j<NCELLJ;j++){
		PHI[j][0]=1;
		PHI_OLD[j][0]=1;
		PHI_previ[j][0]=1;
	}
	for(i=1;i<NCELLI;i++){
		PHI[0][i]=0;
		PHI_OLD[0][i]=0;
		PHI_previ[0][i]=0;
	}
}
void coefficient(){
	for(j=1;j<NCELLJ;j++){
		for(i=1;i<NCELLI;i++){
			aw[j][i]=fabs(DELT*S0[j][i]/Vol[j][i]);
			as[j][i]=fabs(DELT*S1[j][i]/Vol[j][i]);
			ae[j][i]=fabs(DELT*S2[j][i]/Vol[j][i]);
			an[j][i]=fabs(DELT*S3[j][i]/Vol[j][i]);
		}
	}
}

void FOU_SEIDEL(){
	double cc= 0.00000001;
	double rms=1;
	while(rms>cc){
		double sqr=0;
		for(j=1;j<NCELLJ;j++){
			for(i=1;i<NCELLI;i++){
				PHI[j][i] = PHI_OLD[j][i]+(u*aw[j][i]*PHI[j][i-1])-(u*ae[j][i]*PHI[j][i])+(v*as[j][i]*PHI[j-1][i])-(v*an[j][i]*PHI[j][i]);
			}
		}
		for(j=1;j<NCELLJ;j++){
			for(i=1;i<NCELLI;i++){
				sqr=sqr+(PHI[j][i] - PHI_previ[j][i])*(PHI[j][i] - PHI_previ[j][i]);
			}
		}
		rms=sqrt(sqr/(NCELLI*NCELLJ));
		for(j=1;j<NCELLJ;j++){
			for(i=1;i<NCELLI;i++){
				PHI_previ[j][i]=PHI[j][i];
			}
		}
		count=count+1;
	}
}

void SOU_SEIDEL(){
	double cc= 0.00000001;
	double rms=1;
	for(j=1;j<NCELLJ;j++){
		PHI[j][1]=1;
		PHI_OLD[j][1]=1;
		PHI_previ[j][1]=1;
	}
	for(i=1;i<NCELLI;i++){
		PHI[1][i]=0;
		PHI_OLD[1][i]=0;
		PHI_previ[1][i]=0;
	}
	while(rms>cc){
		double sqr=0;
		for(j=2;j<NCELLJ;j++){
			for(i=2;i<NCELLI;i++){
			PHI[j][i] = PHI_OLD[j][i]+(0.5*u*aw[j][i]*((3*PHI[j][i-1])-PHI[j][i-2]))-(0.5*u*ae[j][i]*((3*PHI[j][i])-PHI[j][i-1]))+(0.5*v*as[j][i]*((3*PHI[j-1][i])-PHI[j-2][i]))-(0.5*v*an[j][i]*((3*PHI[j][i])-PHI[j-1][i]));
			}
		}
		for(j=2;j<NCELLJ;j++){
			for(i=2;i<NCELLI;i++){
				sqr=sqr+(PHI[j][i] - PHI_previ[j][i])*(PHI[j][i] - PHI_previ[j][i]);
			}
		}
		rms=sqrt(sqr/(NCELLI*NCELLJ));
		for(j=2;j<NCELLJ;j++){
			for(i=2;i<NCELLI;i++){
				PHI_previ[j][i]=PHI[j][i];
			}
		}
	}
}

void QUICK_SEIDEL(){
	double cc= 0.00000001;
	double rms=1;
	for(j=1;j<NCELLJ;j++){
		PHI[j][1]=1;
		PHI_OLD[j][1]=1;
		PHI_previ[j][1]=1;
	}
	for(i=1;i<NCELLI;i++){
		PHI[1][i]=0;
		PHI_OLD[1][i]=0;
		PHI_previ[1][i]=0;
	}
	while(rms>cc){
		double sqr=0;
		for(j=2;j<NCELLJ;j++){
			for(i=2;i<NCELLI;i++){
				PHI[j][i] = PHI_OLD[j][i]+(u*aw[j][i]*(3*PHI[j][i]+6*PHI[j][i-1]-PHI[j][i-2])/8)-(u*ae[j][i]*(3*PHI[j][i+1]+6*PHI[j][i]-PHI[j][i-1])/8)+(v*as[j][i]*(3*PHI[j][i]+6*PHI[j-1][i]-PHI[j-2][i])/8)-(v*an[j][i]*(3*PHI[j+1][i]+6*PHI[j][i]-PHI[j-1][i])/8);
			}
		}
		for(j=2;j<NCELLJ;j++){
			for(i=2;i<NCELLI;i++){
				sqr=sqr+(PHI[j][i] - PHI_previ[j][i])*(PHI[j][i] - PHI_previ[j][i]);
			}
		}
		rms=sqrt(sqr/(NCELLI*NCELLJ));
		for(j=2;j<NCELLJ;j++){
			for(i=2;i<NCELLI;i++){
				PHI_previ[j][i]=PHI[j][i];
			}
		}
	}
}

void CD_SEIDEL(){
	double cc= 0.00000001;
	double rms=1;
	while(rms>cc){
		double sqr=0;
		for(j=1;j<NCELLJ;j++){
			for(i=1;i<NCELLI;i++){
				PHI[j][i] = PHI_OLD[j][i]+(u*aw[j][i]*(PHI[j][i]+PHI[j][i-1])/2)-(u*ae[j][i]*(PHI[j][i]+PHI[j][i+1])/2)+(v*as[j][i]*(PHI[j][i]+PHI[j-1][i])/2)-(v*an[j][i]*(PHI[j][i]+PHI[j+1][i])/2);
			}
		}
		for(j=1;j<NCELLJ;j++){
			for(i=1;i<NCELLI;i++){
				sqr=sqr+(PHI[j][i] - PHI_previ[j][i])*(PHI[j][i] - PHI_previ[j][i]);
			}
		}
		rms=sqrt(sqr/(NCELLI*NCELLJ));
		for(j=1;j<NCELLJ;j++){
			for(i=1;i<NCELLI;i++){
				PHI_previ[j][i]=PHI[j][i];
			}
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
	ss << "advection_imp_FOU/TEMP_PROFILE_TIME_" << NSTEP*DELT << ".dat";
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
	mkdir("advection_imp_FOU");
	grid_generation();
	initial_condition();
	coefficient();
	for(NSTEP=1;NSTEP<=(Tmax/DELT);NSTEP++){
		boundary_condition();
		FOU_SEIDEL();
		if(NSTEP%200==0){
			write_file();
		}
		update();
	}
	return 0;
}
