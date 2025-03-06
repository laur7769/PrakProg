using static System.Console;
using static System.Math;
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
        WriteLine();
        WriteLine("C. cubic spline");
        WriteLine("The table of data where {x_i=i, y_i=Cos(x_i)}, i=1,...,9 is once again used to test the implemented method.");
        WriteLine("The cubic splines of the implemented method can be sen plotted against the built in cubic spline in cubic.gnuplot.svg.");
        cspline cos2 = new cspline(x, y);
        var data2 = new System.IO.StreamWriter("data_cubic.txt", append:true);
        for(int i=0; i<x.Length; i++){
            data2.WriteLine($"{x[i]}    {y[i]}");
        }
        data2.WriteLine();
        data2.WriteLine();
        for(double i=0; i<1000; i++){
            double z = x[0]+i*(x[x.Length-1]-x[0])/1000;
            data2.WriteLine($"{z}    {cos2.evaluate(z)} {cos2.integral(z)}");
        }
        data2.Close();


        return 0;
    }
}
