#include "Lentz_Algorithm.hpp"
#define epsilon1 1e-30
#define epsilon2 1e-7

double Lentz_Algorithm(double const x)
{
  double f_prev = epsilon1, C_prev = epsilon1, D_prev = 0, delta = 2+epsilon2, D_curr, C_curr, f_curr;
  double j = 1.0, tmp1, tmp2;
  while (delta-1>=epsilon2 || 1-delta >= epsilon2)
  {
    j++;
    tmp1 = x+2*j-1;
    tmp2 = pow(j-1,2);
    D_curr = 1/(tmp1-tmp2*D_prev);
    C_curr = tmp1-tmp2/C_prev;
    delta = C_curr*D_curr;
    f_curr = f_prev*delta;
    f_prev = f_curr;
    C_prev = C_curr;
    D_prev = D_curr;
  }
  return 1/(x+1+f_curr);
}
