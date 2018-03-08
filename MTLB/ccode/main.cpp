#include <math.h>
#include <stdlib.h>

// g++ -03 ./main.cpp -o file
// check -fbounds-check

int main(int argc, char const *argv[])
{
	// solve -(uxx + uxx) = -(x^2+y^2)*exp(xy) with Dirichlet BC
	// true solution: exp(xy)
	// domain [0,1]x[0,1]
	const int n = 20;
	double h = 1/(n+1);

	static double true[20][20];
	for (int j = 0; j < n; ++j)
	{
		int y = h*j;
		for (int i = 0; i < n; ++i)
		{
			double x = h*i;
			true[j][i] = exp(x*y);
		}
	}
	

	static double grid[20][20]; //grid[x][y]
	for (int i = 0; i < n; ++i)
	{
		double x = h*i;
		double y = x;
		grid[0][i] = exp(x*0); // y=0
		grid[n-1][i] = exp(x*1); // y=1
		
		if (i != 0 && i != n-1)
		{
			grid[i][0] = exp(0*y);
			grid[i][n-1] = exp(1*y);
		}
	}



	return 0;
}