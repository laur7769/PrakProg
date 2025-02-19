using static System.Math;
using static System.Console;
public static class sfuns{
   public static double erf(double x){
    /// single precision error function (Abramowitz and Stegun, from Wikipedia)
    if(x<0) return -erf(-x);
    double[] a={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
    double t=1/(1+0.3275911*x);
    double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
    return 1-sum*Exp(-x*x);
    } 

    public static double sgamma(double x){
    if(x<0)return PI/Sin(PI*x)/sgamma(1-x);
    if(x<9)return sgamma(x+1)/x;
    double lnsgamma=Log(2*PI)/2+(x-0.5)*Log(x)-x
        +(1.0/12)/x-(1.0/360)/(x*x*x)+(1.0/1260)/(x*x*x*x*x);
    return Exp(lnsgamma);
    }

    public static double lngamma(double x){
    if(x<=0) throw new System.ArgumentException("lngamma: x<=0");
    if(x<9) return lngamma(x+1)-Log(x);
    return x*Log(x+1/(12*x-1/x/10))-x+Log(2*PI/x)/2;
    }
}