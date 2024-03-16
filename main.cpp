//Semi explicit Finite volume method
//using Navier stokes solver
// for lid driven cavity
//Both upwind and quick schemes
//---------------------------------------------------------------------------


#include <iostream>
# include <fstream>
# include <cmath>

using namespace std;
const int N=102;         // Grid points
const double Re=100.0;   // Reynolds Number
#include"pressure.h"
#include"uvupwind.h"
#include"uvquick.h"
#include"uvnew.h"

int main()
{
    int i,j, n =0, Nreal, choice;
    Nreal=(N-2)*(N-2);
    float dt=0.001, dx=0.01, dy=0.01, volp, r, rms;  // dt time step dx, dy grid size
    volp=dx*dy;
    double  up[N][N] ={0}, vp[N][N]={0}, unew[N][N]={0},  vnew[N][N] ={0}, u[N][N]={0},  v[N][N]={0};
    double p[N][N]={0}, fe[N][N]={0}, fw[N][N]={0}, fn[N][N]={0}, fs[N][N]={0};
    double U[N],V[N];


    cout<<"Enter 1 for upwind scheme or 2 for quick scheme "<<endl;
    cin>>choice;

     do
    {
     rms=0;

     if (choice == 1)
        uvupwind( u, v, up, vp, fe, fw, fn, fs, dx, dy, volp, dt);
     else if (choice ==2)
        {
          uvquick( u, v, up, vp, fe, fw, fn, fs, dx, dy, volp, dt);
        }

 //pressure Poisson equation using gauss Seidel iterations----------

pressure(p, up, vp, dx, dy, Nreal, dt);

//Flux calculation for each face --------
  for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
        fe[i][j] = (((up[i][j]+up[i+1][j])/2.0)-dt*(p[i+1][j]-p[i][j])/dx)*dy ;
        fw[i][j] = -(((up[i][j]+up[i-1][j])/2.0)-dt*(p[i][j]-p[i-1][j])/dx)*dy ;
        fn[i][j] = (((vp[i][j]+vp[i][j+1])/2.0)-dt*(p[i][j+1]-p[i][j])/dy)*dx ;
        fs[i][j] = -(((vp[i][j]+vp[i][j-1])/2.0)-dt*(p[i][j]-p[i][j-1])/dy)*dx ;
        }
        }
uvnew( unew, vnew, up, vp, p, dx, dy, dt);

    for( i=0;i<N;++i)
    {for(j=0;j<N;++j)
    {   r=pow((unew[i][j]-u[i][j]),2)+pow((vnew[i][j]-v[i][j]),2);
        rms=rms+r;
        }}

rms = sqrt(rms/Nreal);

    for( i=0;i<N;++i)
    for(j=0;j<N;++j)
        {
            u[i][j] = unew[i][j];
            v[i][j]=vnew[i][j];
        }

    n=n+1;
    cout<<"n: "<<n<<" rms: "<<rms<<endl;
    if(n>10000 || rms >10000){
	    cout <<" convergence Error, Maximum iterations" << n <<" reached"<<endl;
	    return 1;
    }
}
while ((rms/dt)>1e-3);        //steady state solution criteria
cout << " Solution Coverged " <<"\n";

//plotting center line velocities

    int cord1 = (N-2)/2;
    int cord2 = cord1+1;
    for(j=0;j<N;++j)
    {
       U[j] = (u[cord1][j]+u[cord2][j])/2.0;
       V[j] = (v[j][cord1]+v[j][cord2])/2.0;
    }
    U[0]=(U[0]+U[1])/2.0;
    U[N-1]=(U[N-2]+U[N-1])/2.0;
    V[0]=(V[0]+V[1])/2.0;
    V[N-1]=(V[N-2]+V[N-1])/2.0;



 //writing data to text file

  ofstream myfile;
  myfile.open ("u_v_data.txt");
  myfile << "Writing U center vertical line velocity to a file.\n";
  for (j=0;j<N;++j)
        myfile << U[j] <<",";
  myfile << "\n Writing V center horizontal line velocity to a file.\n";
  for(j=0;j<N;++j)
        myfile << V[j]<< ",";
    myfile.close();


   myfile.open("u_data.txt");
    myfile << "\n Writing u velocity to file.\n";
    for (j=0;j<N;++j)
    {
        myfile <<"\n";
      for(i=0;i<N;++i)
        myfile <<u[i][j] <<",";
    }
   myfile.close();

       myfile.open("v_data.txt");
    myfile << "\n Writing v velocity to file.\n";
    for (j=0;j<N;++j)
    {
        myfile <<"\n";
      for(i=0;i<N;++i)
        myfile <<v[i][j] <<",";
    }
   myfile.close();

return 0;
    }


























