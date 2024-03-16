#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
const double Re=100;

void uvupwind( double u[N][N], double v[N][N], double up[N][N], double vp[N][N],double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N], float dx, float dy, float volp, float dt)
{
   int i,j;
   double ue,ve,uw,vw,vn,un,us,vs,tfu,tfv;
   double fdue,fdve,fduw,fdvw,fdun,fdvn,fdus,fdvs,fduf,fdvf;

    // upwind scheme
    for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
 if (fe[i][j]>=0 )
    {
    ue=u[i][j];
    ve=v[i][j];}
 else
    {
    ue=u[i+1][j];
    ve=v[i+1][j];}

 if (fw[i][j]>=0 )
    {
    uw=u[i][j];
    vw=v[i][j];}
 else
   {
    uw=u[i-1][j];
    vw=v[i-1][j]; }

 if (fn[i][j]>=0 )
    {
    un=u[i][j];
    vn=v[i][j];}
 else
   {
    un=u[i][j+1];
    vn=v[i][j+1];}

 if (fs[i][j]>=0 )
   {
    us=u[i][j];
    vs=v[i][j];}
 else
   {
    us=u[i][j-1];
    vs=v[i][j-1];}

     tfu=fe[i][j]*ue + fw[i][j]*uw + fn[i][j]*un + fs[i][j]*us;

     tfv=fe[i][j]*ve + fw[i][j]*vw + fn[i][j]*vn + fs[i][j]*vs;

// diffusive flux

fdue = -(dy/(dx*Re))*(u[i+1][j]-u[i][j]);
fduw = -(dy/(dx*Re))*(u[i-1][j]-u[i][j]);
fdun = -(dx/(dy*Re))*(u[i][j+1]-u[i][j]);
fdus = -(dx/(dy*Re))*(u[i][j-1]-u[i][j]);
fduf=fdue+fduw+fdun+fdus;

fdve = -(dy/(dx*Re))*(v[i+1][j]-v[i][j]);
fdvw = -(dy/(dx*Re))*(v[i-1][j]-v[i][j]);
fdvn = -(dx/(dy*Re))*(v[i][j+1]-v[i][j]);
fdvs = -(dx/(dy*Re))*(v[i][j-1]-v[i][j]);
fdvf=fdve+fdvw+fdvn+fdvs;

// u* and v* velocities
up[i][j]=u[i][j]-(dt*(tfu+fduf))/volp;
vp[i][j]=v[i][j]-(dt*(tfv+fdvf))/volp;
}
}
//boundary conditions ---- fictitious cells

    for (j=1;j<N-1;++j)    //left wall
    {
        up[0][j] = -up[1][j];
        vp[0][j] = -vp[1][j];
    }

    for (j=1;j<N-1;++j)      // right wall
    {
        up[N-1][j] = -up[N-2][j];
        vp[N-1][j] = -vp[N-2][j];
    }

    for (j=1;j<N-1;++j)      // bottom wall
    {
        up[j][0] = -up[j][1];
        vp[j][0] = -vp[j][1];
    }

    for (j=1;j<N-1;++j)    // top wall
    {
        up[j][N-1] = 2-up[j][N-2];
        vp[j][N-1] = -vp[j][N-2];
    }
}

