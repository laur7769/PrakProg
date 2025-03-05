using static System.Console;
using static System.Math;
public class qspline {
	public vector x,y,b,c;
	public qspline(vector xs,vector ys){
		x = xs.copy(); 
        y = ys.copy();
        int n = x.size; 
        b = new vector(n-1);
        c = new vector(n-1);
        c[0] = 0;
        vector dx = new vector(n-1);
        for(int i=0; i<dx.size; i++){
            dx[i] = x[i+1]-x[i];
        }
        vector dy = new vector(n-1);
        for(int i=0; i<dy.size; i++){
            dy[i] = y[i+1]-y[i];
        }
        vector p = new vector(n-1);
        for(int i=0; i<p.size; i++){
            p[i] = dy[i]/dx[i];
        }
        for(int i=0; i<n-2; i++ ){
            c[i+1] = (1/dx[i+1])*(p[i+1]-p[i]-c[i]*dx[i]);
        }//forward substitution
        c[n-2]/=2;
        for(int i=n-3; i>=0; i--){
            c[i]=(p[i+1]-p[i]-c[i+1]*dx[i+1])/dx[i];
        }
        for(int i=0; i<b.size; i++){
            b[i]=p[i]-c[i]*dx[i];
        }
	}
	public double evaluate(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        return y[j]+b[j]*(z-x[j])+c[j]*(Pow(z-x[j],2));
        }
    public static int binsearch(double[] x, double z){ 
	    if( z<x[0] || z>x[x.Length-1] ) throw new System.Exception("binsearch: bad z");
	    int i=0, j=x.Length-1;
	    while(j-i>1){
		    int mid=(i+j)/2;
		    if(z>x[mid]) i=mid; else j=mid;
		}
	    return i;
	}
	public double derivative(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        return b[j]+2*c[j]*(z-x[j]);
    }
	public double integral(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        double integral = 0;
        for(int i=0; i<j; i++){
            integral += y[i]*(x[i+1]-x[i])+b[i]*Pow((x[i+1]-x[i]),2)/2.0+c[i]*Pow((x[i+1]-x[i]),3)/3.0;
        }
        integral += y[j]*(z-x[j])+b[j]*Pow((z-x[j]),2)/2.0+c[j]*Pow((z-x[j]),3)/3.0;
        return integral;
    }
}
class main{
    public static double linterp(double[] x, double[] y, double z){
        int i=binsearch(x,z);
        double dx=x[i+1]-x[i]; if(!(dx>0)) throw new System.Exception("uups...");
        double dy=y[i+1]-y[i];
        return y[i]+dy/dx*(z-x[i]);
        }
    public static int binsearch(double[] x, double z)
	{/* locates the interval for z by bisection */ 
	if( z<x[0] || z>x[x.Length-1] ) throw new System.Exception("binsearch: bad z");
	int i=0, j=x.Length-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}
    public static double linterpInteg(double[] x, double[] y, double z){
        int j=binsearch(x,z);
        double integral = 0;
        for(int i=0; i<j; i++){
            integral += y[i]*(x[i+1]-x[i])+((y[i+1]-y[i])/(x[i+1]-x[i]))*Pow((x[i+1]-x[i]),2)*(1.0/2.0);
        }
        integral += y[j]*(z-x[j])+((y[j+1]-y[j])/(x[j+1]-x[j]))*Pow((z-x[j]),2)*(1.0/2.0);
        return integral;
    }

    static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Linear spline (linear interpolation)");
        WriteLine(" The following table of data is used to test the implemented methods:");
        WriteLine("x   y=Cos(x)");
        var data = new System.IO.StreamWriter("data.txt", append:true);
        double[] x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        double[] y = {Cos(0.0), Cos(1.0), Cos(2.0), Cos(3.0), Cos(4.0), Cos(5.0), Cos(6.0), Cos(7.0), Cos(8.0), Cos(9.0)};
        for(int i=0; i<x.Length; i++){
            WriteLine($"{x[i]} {y[i]}");
            data.WriteLine($"{x[i]} {y[i]}");
        }
        data.WriteLine();
        data.WriteLine();
        for(double i=0; i<1000; i++){
            double z = x[0]+i*(x[x.Length-1]-x[0])/1000;
            data.WriteLine($"{z}    {linterp(x, y, z)} {linterpInteg(x, y, z)}");
        }
        data.WriteLine();
        data.WriteLine();
        WriteLine("The linear interpolation as well as the given data points can be seen plotted in interp.gnuplot.svg");
        WriteLine("The definite integral of the interpolant can be seen plotted in int.gnuplot.svg");
        WriteLine();
        WriteLine("B. Quadratic spline");
        WriteLine("Check quadratic spline:");
        WriteLine("For data {x_i=i, y_i=1}, i=1,...,5 ");
        double[] c_data = {0.0, 0.0, 0.0, 0.0};
        double[] x_data = {1.0, 2.0, 3.0, 4.0, 5.0};
        double[] y_data = {1.0, 1.0, 1.0, 1.0, 1.0};
        double[] b_data = {0.0, 0.0, 0.0, 0.0};
        vector x_1 = new vector(x_data);
        vector y_1 = new vector(y_data);
        vector c_1 = new vector(c_data);
        vector b_1 = new vector(b_data);
        qspline one = new qspline(x_1, y_1);
        WriteLine($"Is manually calculated c's equal to computed c's? {c_1.approx(one.c)}");
        
        WriteLine($"Is manually calculated b's equal to computed b's? {b_1.approx(one.b)}");
        WriteLine("For data {x_i=i, y_i=x_i}, i=1,...,5 ");
        double[] c_data2 = {0.0, 0.0, 0.0, 0.0};
        double[] x_data2 = {1.0, 2.0, 3.0, 4.0, 5.0};
        double[] y_data2 = {1.0, 2.0, 3.0, 4.0, 5.0};
        double[] b_data2 = {1.0, 1.0, 1.0, 1.0};
        vector x_2 = new vector(x_data2);
        vector y_2 = new vector(y_data2);
        vector c_2 = new vector(c_data2);
        vector b_2 = new vector(b_data2);
        qspline two = new qspline(x_2, y_2);
        WriteLine($"Is manually calculated c's equal to computed c's? {c_2.approx(two.c)}");
        
        WriteLine($"Is manually calculated b's equal to computed b's? {b_2.approx(two.b)}");

        WriteLine();
        WriteLine("The table of data where {x_i=i, y_i=Cos(x_i)}, i=1,...,9 is once again used to test the implemented method.");
        WriteLine("The quadratic splines can be sen plotted in interp.gnuplot.svg.");
        WriteLine("The definite integral of the interpolant can be seen in int.gnuplot.svg.");
        qspline cos = new qspline(x, y);
        for(double i=0; i<1000; i++){
            double z = x[0]+i*(x[x.Length-1]-x[0])/1000;
            data.WriteLine($"{z}    {cos.evaluate(z)}   {cos.integral(z)}");
        }
        data.Close();
        return 0;
    }
}
