// DEVAM SHETH 18BME109
//HEADER FILES
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<dir.h>
#include<sstream>

using namespace std;

//DEFINITION OF CONSTANT
#define NI 102
#define NJ 102
#define pi 3.141592
#define t_max 5.5
#define L 1
#define rho  1000
#define Cp 1
#define H 1
#define alpha 0.1
#define DELT 0.01
#define Sp 0
// USER-DEFINED FUNCTION PROTOTYPE
void grid_generation();
void initial_condition();
void boundary_condition();
void coefficient();
void solve_implicit();
void update();
void write();
//GLOBAL VARIABLES
double x0[NJ][NI],x1[NJ][NI],x2[NJ][NI],x3[NJ][NI];
double y_0[NJ][NI],y_1[NJ][NI],y2[NJ][NI],y3[NJ][NI];

double XCELL_X[NJ][NI], XCELL_Y[NJ][NI];
int NCELLI = NI - 1, NCELLJ = NJ - 1;

double s0[NJ][NI], s1[NJ][NI], s2[NJ][NI], s3[NJ][NI];
double DXF0[NJ][NI], DXF1[NJ][NI], DXF2[NJ][NI], DXF3[NJ][NI];
double VOL[NJ][NI];
double PHI[NJ][NI], PHI_OLD[NJ][NI],PHI_ITERAT[NJ][NI];
double aw[NJ][NI], as[NJ][NI], ae[NJ][NI], an[NJ][NI], ap[NJ][NI];
int NSTEP=1;
double rms,sq,sum;
int i, j,count;
int main()
{
    mkdir("RUN1");
	int step=1;
    grid_generation();
    initial_condition();
    coefficient();


//do
for(NSTEP=1;NSTEP<=(t_max/DELT);NSTEP++)
{
    count=0;
 //NSTEP=NSTEP + 1;
 boundary_condition();
 solve_implicit();
 update();

  if( NSTEP% step == 0)
       {
          write();
          printf("%d \n",count) ;
       }
}
//while(NSTEP<=t_max);
}
void grid_generation()
{
// Constant del_x and del_y
double DELX = L / (NI - 2.0),DELY = H / (NJ - 2.0);

// Co-ordinates of interior nodes
for (j = 1;j < NCELLJ;j++)
{
for (i = 1;i < NCELLI;i++)
 {
    x0[j][i] = DELX * (i - 1);
    x1[j][i] = DELX * i;
    x2[j][i] = DELX * i;
    x3[j][i] = DELX * (i - 1);
    y_0[j][i] = DELX * (j - 1);
    y_1[j][i] = DELX * (j - 1);
    y2[j][i] = DELX * j;
    y3[j][i] = DELX * j;
    XCELL_X[j][i] = 0.25 * (x0[j][i] + x1[j][i] + x2[j][i] + x3[j][i]);
    XCELL_Y[j][i] = 0.25 * (y_0[j][i] + y_1[j][i] + y2[j][i] + y3[j][i]);
 }
}
// Co-ordinates of Boundary nodes : W S E N (anti-clockwise)
for (j = 1;j < NCELLJ;j++)
  {
        // West wall
    XCELL_X[j][0] = 0;
    XCELL_Y[j][0] = XCELL_Y[j][1]; // same as centroids next to west wall
    // East wall
    XCELL_X[j][NCELLI] = L;
    XCELL_Y[j][NCELLI] = XCELL_Y[j][NCELLI-1]; // same as centroids previous to east wall
  }
for (i = 1;i < NCELLI;i++)
{
    // South wall
    XCELL_X[0][i] = XCELL_X[1][i];
    XCELL_Y[0][i] = 0; // same as centroids above to south wall
    // North wall
    XCELL_X[NCELLJ][i] = XCELL_X[NCELLJ-1][i];
    XCELL_Y[NCELLJ][i] = H; // same as centroids below to north wall

    //corners
    XCELL_X[NCELLJ][NCELLI]=L;
    XCELL_Y[NCELLJ][NCELLI]=H;
    XCELL_X[0][NCELLI]=L;
    XCELL_X[NCELLJ][0]=0;
}

// surface area & nodal distance
 for (j = 1;j < NCELLJ;j++)
  {
    for (i = 1;i < NCELLI;i++)
    {
     s0[j][i] = -DELY;
     s1[j][i] = -DELX;
     s2[j][i] = DELY;
     s3[j][i] = DELX;
     DXF0[j][i] = DELX;
     DXF1[j][i] = DELY;
     DXF2[j][i] = DELX;
     DXF3[j][i] = DELY;
     VOL[j][i] = DELX * DELY;
    }
  }
}
void initial_condition()
{
for (j = 0;j <= NCELLJ;j++)
{
for (i = 0;i <= NCELLI;i++)
{
PHI_OLD[j][i] = 100.0;
PHI[j][i] = 100.0;
PHI_ITERAT[j][i]=100.0;
}
}
}
void boundary_condition()
{
for (j = 0;j <=NCELLJ;j++)
{
// West wall
PHI[j][0] = 300.0; // dirichlet boundary condition
PHI_OLD[j][0] = 300.0;
//PHI_ITERAT[j][i]=300.0;
// East wall
PHI[j][NCELLI] = 100.0; // dirichlet boundary condition
PHI_OLD[j][NCELLI] = 100.0;
//PHI_ITERAT[j][NCELLI]=400.0;
}
for (i = 0;i <=NCELLI;i++)
{
// South wall
PHI[0][i] = 200.0; // dirichlet boundary condition
PHI_OLD[0][i] = 200.0;
//PHI_ITERAT[0][i]=200.0;

// North wall
PHI[NCELLJ][i] = 150.0;// dirichlet boundary condition
PHI_OLD[NCELLJ][i] = 150.0;
//PHI_ITERAT[NCELLJ][i]=150.0;
}

}
void coefficient()
{
for (j = 1;j < NCELLJ;j++)
{
for (i = 1;i < NCELLI;i++)
{
aw[j][i] = fabs((alpha*s0[j][i]*DELT)/(VOL[j][i]*DXF0[j][i]));
as[j][i] = fabs((alpha*s1[j][i]*DELT)/(VOL[j][i]*DXF1[j][i]));
ae[j][i] = fabs((alpha*s2[j][i]*DELT)/(VOL[j][i]*DXF2[j][i]));
an[j][i] = fabs((alpha*s3[j][i]*DELT)/(VOL[j][i]*DXF3[j][i]));
ap[j][i] = aw[j][i] + as[j][i] + ae[j][i] + an[j][i];

}
}
}
void solve_implicit()
{

    rms=2.0;

while(rms>0.00000001)
{
    sum=0;
    for(j = 1;j < NCELLJ;j++)
  {
    for(i = 1;i < NCELLI;i++)
    {
      PHI[j][i] = ( (aw[j][i]*PHI[j][i-1]) + (as[j][i]*PHI[j-1][i]) + (ae[j][i]*PHI[j][i+1]) + (an[j][i]*PHI[j+1][i]) + (Sp*DELT) + PHI_OLD[j][i] )/(1+ap[j][i]) ;
    }
  }
    for(j = 1;j < NCELLJ;j++)
    {
        for(i = 1;i < NCELLI;i++)
        {
            sq = (PHI[j][i]- PHI_ITERAT[j][i])*(PHI[j][i]- PHI_ITERAT[j][i]);
            sum += sq;
        }
    }

     rms=(sqrt((sq)/(NCELLI*NCELLJ)));
      //	RMS= (sqrt((sum)/(NCELLI*NCELLJ)))
// for (j = 1;j < NCELLJ;j++)
//  {
//     for (i = 1;i < NCELLI;i++)
//    {
//        if(fabs(PHI[j][i]-PHI_ITERAT[j][i]) > sq)
//                sq=PHI[j][i]-PHI_ITERAT[j][i];
//    }
//  }
      for(j = 1;j < NCELLJ;j++)
       {
        for(i = 1;i < NCELLI;i++)
        {
            PHI_ITERAT[j][i] = PHI[j][i] ;
            //printf("C");
        }
       }
      count+=1;
}

}

void update()
{

      for(j = 1;j < NCELLJ;j++)
      {
        for(i = 1;i < NCELLI;i++)
         {

            PHI_OLD[j][i]=PHI[j][i];

         }
      }

}

void write()
{
	ofstream file,f1;
	stringstream ss;
	ss << "RUN1/TEMP_PROFILE_TIME_" << NSTEP*DELT << ".dat";
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
	f1.open("plot1.csv",ios::out | ios::app);
	f1<<NSTEP*DELT<<","<<count<<"\n";

}



