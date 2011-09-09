module calc_square_mod;

double calc_square(in double x)
{
  return(x*x);
}

extern(C) void R_calc_square(ref double x, ref double output)
{
  output = calc_square(x);
} 

