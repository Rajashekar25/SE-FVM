#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
double Re=100;

// corrected velocities

void uvnew( double unew[N][N], double vnew[N][N], double up[N][N], double vp[N][N], double p[N][N], float dx, float dy, float dt)
{
   int i,j;

    for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
            unew[i][j]=up[i][j]-dt*(p[i+1][j]-p[i-1][j])/(2.0*dx);
            vnew[i][j]=vp[i][j]-dt*(p[i][j+1]-p[i][j-1])/(2.0*dy);
        }
    }
         for (j=1;j<N-1;++j)    //left wall
    {
        unew[0][j] = -unew[1][j];
        vnew[0][j] = -vnew[1][j];
    }

    for (j=1;j<N-1;++j)      // right wall
    {
        unew[N-1][j] = -unew[N-2][j];
        vnew[N-1][j] = -vnew[N-2][j];
    }

    for (j=1;j<N-1;++j)      // bottom wall
    {
        unew[j][0] = -unew[j][1];
        vnew[j][0] = -vnew[j][1];
    }

    for (j=1;j<N-1;++j)    // top wall
    {
        unew[j][N-1] = 2.0-unew[j][N-2];
        vnew[j][N-1] = -vnew[j][N-2];
    }
}
