using static System.Console;
using static System.Math;
using System.Collections.Generic;
public class genlist<T>{
	public T[] data;
	public int size => data.Length; // property
	public T this[int i] => data[i]; // indexer
	public genlist(){ data = new T[0]; }
	public void add(T item){ /* add item to the list */
		T[] newdata = new T[size+1];
		System.Array.Copy(data,newdata,size);
		newdata[size]=item;
		data=newdata;
	}
}
class main{
    
    public static (vector,vector) rkstep23(
        // An implementation of a simple embedded method of order 2 and 3
	System.Func<double,vector,vector> f,/* the f from dy/dx=f(x,y) */
	double x,                    /* the current value of the variable */
	vector y,                    /* the current value y(x) of the sought function */
	double h                     /* the step to be taken */
	){
	    vector k0 = f(x,y);              /* embedded lower order formula (Euler) */
	    vector k1 = f(x+h/2,y+k0*(h/2)); /* lower order formula (midpoint) */
        vector k2 = f(x+(h*3)/4, y+((h*3)/4)*k1);
	    vector yh = y+h*(2.0/9.0)*k0+h*(3.0/9.0)*k1+h*(4.0/9.0)*k2;             /* y(x+h) estimate */
	    vector δy = yh-(y+k1*h);           /* error estimate */
	    return (yh,δy);
        }//stepper 23

public static (vector, vector) rkstep23_ny (System.Func<double, vector, vector> f, double x, vector y, double h) {
        // This is an implementation of the Bogacki-Shampine method
        vector k0 = f(x, y);
        vector k1 = f(x + h/2.0, y + k0*(h/2.0));
        vector k2 = f(x + 3.0/4.0*h, y + 3.0/4.0*h*k1);
        vector k3 = f(x + h, y + 2.0/9.0*h*k0 + 1.0/3.0*h*k1 + 4.0/9.0*h*k2);
        vector k = 2.0/9.0*k0 + 1.0/3.0*k1 + 4.0/9.0*k2;

        vector yh = y + k*h;

        vector dy = (2.0/9.0 - 7.0/24.0)*h*k0 + (1.0/3.0-1.0/4.0)*h*k1 + (4.0/9.0 - 1.0/3.0)*h*k2 -1.0/8.0*h*k3;

        return (yh, dy);
    }


    public static (genlist<double>, genlist<vector>, genlist<double>) driver(
    System.Func<double, vector, vector> F,/* the f from dy/dx=f(x,y) */
    (double, double) interval,    /* (initial-point,final-point) */
    vector yinit,                /* y(initial-point) */
    double h = 0.125,              /* initial step-size */
    double acc = 0.01,             /* absolute accuracy goal */
    double eps = 0.01              /* relative accuracy goal */
    )
    {
        var (a, b) = interval; double x = a; vector y = yinit.copy();
        var xlist = new genlist<double>(); xlist.add(x);
        var ylist = new genlist<vector>(); ylist.add(y);
        var hlist = new genlist<double>(); hlist.add(h);
        do
        {
            if (x >= b) return (xlist, ylist, hlist); /* job done */
            if (x + h > b) h = b - x;               /* last step should end at b */
            var (yh, δy) = rkstep23_ny(F, x, y, h);
            double tol = (acc + eps * yh.norm()) * Sqrt(h / (b - a));
            double err = δy.norm();
            if (err <= tol)
            { // accept step
                x += h; y = yh;
                xlist.add(x);
                ylist.add(y);
            }
            if (err > 0) h *= Min(Pow(tol / err, 0.25) * 0.95, 2); // readjust stepsize
            else h *= 2;
            hlist.add(h);
        } while (true);
    }//driver
    
    public static System.Func<double,vector> make_ode_ivp_qspline
    (System.Func<double,vector,vector> F,
    (double, double) interval,
    vector y,
    double acc=0.01,
    double eps=0.01,
    double hstart=0.01 ){
        System.Func<double, vector> make_qspline = (z) =>
        {
            var (xlist, ylist, hlist) = driver(F, interval, y, acc, eps, hstart);
            vector eval = new vector(ylist[0].size);
            vector x = new vector(xlist.size);
            for (int i = 0; i < xlist.size; i++)
            {
                x[i] = xlist[i];
            }
            for (int t = 0; t < ylist[0].size; t++)
            {
                vector y_ny = new vector(ylist.size);
                for (int i = 0; i < y_ny.size; i++)
                {
                    y_ny[i] = ylist[i][t];
                }
                qspline spline = new qspline(x, y_ny);
                eval[t]=(spline.evaluate(z));
            }
            return eval;
        };
        
	    return make_qspline;
    }
    static int Main(string[] args)
    {
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Embedded rule Runge-Kutta ODE integrator");
        WriteLine("     - Function for debugging: u'=-u, initial value: (0,5)");
        WriteLine("     - The stepper used is an implemented Bogacki-Shampine 23 stepper, and default values for acc, eps, and hstart were used.");
        WriteLine("     - The computed results using the implemented methods are plotted in first_ODE.gnuplot.svg along with the analytical function.");
        WriteLine();
        WriteLine();
        System.Func<double, vector, vector> f = (x, y) => -y;
        genlist<double> xlist = driver(f, (0.0, 1.0), new vector(5.0)).Item1;
        genlist<vector> ylist = driver(f, (0.0, 1.0), new vector(5.0)).Item2;
        WriteLine("x    y");
        for (int i = 0; i < xlist.size; i++)
        {
            WriteLine($"{xlist[i]}  {ylist[i][0]}");
        }
        WriteLine();
        WriteLine("  Oscillator with friction:");
        WriteLine("     -The computed theta(t) and omega(t) can be seen in osc_fric.gnuplot.svg");
        double b = 0.25;
        double c = 5.0;
        System.Func<double, vector, vector> omega_first = (x, theta_omega) => new vector(theta_omega[1], -b * theta_omega[1] - c * Sin(theta_omega[0]));
        genlist<double> tlist = driver(omega_first, (0.0, 10.0), new vector(PI - 0.1, 0.0)).Item1;
        genlist<vector> theta_omega_vec = driver(omega_first, (0.0, 10.0), new vector(PI - 0.1, 0.0)).Item2;
        WriteLine("t    theta(t)    omega(t)");
        WriteLine();
        WriteLine();
        for (int i = 0; i < tlist.size; i++)
        {
            WriteLine($"{tlist[i]}  {theta_omega_vec[i][0]}  {theta_omega_vec[i][1]}");
        }
        WriteLine();
        WriteLine("B. Relativistic precession of planetary orbit");
        WriteLine("     - Once again the implemented Bogacki-Shampine 23 stepper is used for calculations.");
        WriteLine("     - acc was set to 0.0001 and eps was set to 0.0001.");
        WriteLine("     - The orbits are plotted in planet.gnuplot.svg");
        var data = new System.IO.StreamWriter("data_B.txt", append: true);
        double[] eps = { 0.0, 0.0, 0.01 };
        vector[] initials = { new vector(1.0, 0.0), new vector(1.0, -0.5), new vector(1.0, -0.5) };
        for (int i = 0; i < eps.Length; i++)
        {
            System.Func<double, vector, vector> u_double = (x, y_u) => new vector(y_u[1], 1 - y_u[0] + eps[i] * y_u[0] * y_u[0]);
            genlist<double> phi = driver(u_double, (0.0, 32 * PI), initials[i], acc:0.0001, eps:0.0001).Item1;
            genlist<vector> u_ufirst = driver(u_double, (0.0, 32 * PI), initials[i], acc:0.0001, eps:0.0001).Item2;
            for (int t = 0; t < phi.size; t++)
            {
                data.WriteLine($"{phi[t]}   {u_ufirst[t][0]}    {u_ufirst[t][1]}");
            }
            data.WriteLine();
            data.WriteLine();
        }
        data.Close();
        WriteLine();
        WriteLine("C. Test of the order of the method; Alternative interface; Newtonian gravitational three body problem");
        WriteLine("1. The numerical results of the integration of y''=2x using a 23 stepper, as well as the numerical results is plotted in 23.gnuplot.svg ");
        WriteLine("     - From analytical integration: y''=2x,   y'=x^2,  y=1/3x^3");
        WriteLine("     - Based on the above the initial conditions must be y'(0)=0 and y(0)=0");
        WriteLine("     - A 23 stepper has already been implemented, and the Bogacki-Sherpine stepper mentioned earlier is once again used.");
        WriteLine("     - Default acc, eps, and hstart is used.");
        WriteLine("     - Step sizes h are printed below:");

        System.Func<double, vector, vector> C1 = (x, y_yfirst) => new vector(y_yfirst[1], 2 * x);
        var result = driver(C1, (0.0, 10.0), new vector(0, 0));
        genlist<double> xs = result.Item1;
        genlist<vector> f_ffirst = result.Item2;
        var data2 = new System.IO.StreamWriter("data_C.txt", append: true);
        for (int t = 0; t < xs.size; t++)
        {
            data2.WriteLine($"{xs[t]}   {f_ffirst[t][0]}    {f_ffirst[t][1]}");
        }
        data2.WriteLine();
        data2.WriteLine();
        genlist<double> hs = result.Item3;
        WriteLine();
        for (int i = 0; i < hs.size; i++)
        {
            WriteLine($"{hs[i]}");
        }
        WriteLine("     - As seen above, the step size does not double, indicating an error above 0 in each step. However, as can be seen in 23.gnuplot.svg, the numerical calculations show very good agreement with the analytical result. ");
        WriteLine(" 2. A quadratic spline interpolation routine has been implemented in the splines.cs script.");
        WriteLine(" 3. The interface of the driver has been adjusted, so it returns the quadratic spline of the table");
        System.Func<double, vector, vector> three_body = (t, z) => {
            // z=(x_1', y_1', x_2', y_2', x_3', y_2', x_1, y_1, x_2, y_2, x_3, y_3)
            vector z_mark = new vector(12);
            z_mark[0] = (z[8] - z[6]) / Pow(Pow(z[8] - z[6], 2) + Pow(z[9] - z[7], 2), 3.0 / 2.0) + (z[10] - z[6]) / Pow(Pow(z[10] - z[6], 2) + Pow(z[11] - z[7], 2), 3.0 / 2.0);
            z_mark[1] = (z[9] - z[7]) / Pow(Pow(z[8] - z[6], 2) + Pow(z[9] - z[7], 2), 3.0 / 2.0) + (z[11] - z[7]) / Pow(Pow(z[10] - z[6], 2) + Pow(z[11] - z[7], 2), 3.0 / 2.0);
            z_mark[2] = (z[6] - z[8]) / Pow(Pow(z[6] - z[8], 2) + Pow(z[7] - z[9], 2), 3.0 / 2.0) + (z[10] - z[8]) / Pow(Pow(z[10] - z[8], 2) + Pow(z[11] - z[9], 2), 3.0 / 2.0);
            z_mark[3] = (z[7] - z[9]) / Pow(Pow(z[6] - z[8], 2) + Pow(z[7] - z[9], 2), 3.0 / 2.0) + (z[11] - z[9]) / Pow(Pow(z[10] - z[8], 2) + Pow(z[11] - z[9], 2), 3.0 / 2.0);
            z_mark[4] = (z[6] - z[10]) / Pow(Pow(z[6] - z[10], 2) + Pow(z[7] - z[11], 2), 3.0 / 2.0) + (z[8] - z[10]) / Pow(Pow(z[8] - z[10], 2) + Pow(z[9] - z[11], 2), 3.0 / 2.0);
            z_mark[5] = (z[7] - z[11]) / Pow(Pow(z[6] - z[10], 2) + Pow(z[7] - z[11], 2), 3.0 / 2.0) + (z[9] - z[11]) / Pow(Pow(z[8] - z[10], 2) + Pow(z[9] - z[11], 2), 3.0 / 2.0);
            z_mark[6] = z[0];
            z_mark[7] = z[1];
            z_mark[8] = z[2];
            z_mark[9] = z[3];
            z_mark[10] = z[4];
            z_mark[11] = z[5];
            return z_mark;
        };
        vector initial = new vector(0.4662036850, 0.4323657300, -0.93240737, -0.86473146, 0.4662036850, 0.4323657300, -0.97000436, 0.24308753, 0.0, 0.0, 0.97000436, -0.24308753);
        var result_three_body = make_ode_ivp_qspline(three_body, (0.0, 6.32591398 / 2.0), initial, acc: 0.0001, eps: 0.0001, hstart:0.001);
        
        for (double i = 0; i < 1000; i++) {
            double z = 0.0 + i * ((6.32591398 / 2.0) / 1000);
            data2.WriteLine($"{z}    {result_three_body(z)[6]}   {result_three_body(z)[7]}  {result_three_body(z)[8]}   {result_three_body(z)[9]}    {result_three_body(z)[10]}   {result_three_body(z)[11]} ");
        }
        data2.WriteLine();
        data2.WriteLine();
        for (int i = 0; i < 3; i++)
        {
            data2.WriteLine($"{initial[6 + 2*i]}    {initial[7 + 2 * i]}");
        }
        data2.Close();
        WriteLine("4. The periodic solution to the three body problem is plotted in three_body.gnuplot.svg.");
        WriteLine("     - For numerical calculations using the implemeted driver, the acc was set to 0.0001, eps was set to 0.0001 and the starting step soze set to 0.001.");

        return 0;
    }
}
