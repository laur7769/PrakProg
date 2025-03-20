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

    public static (genlist<double>,genlist<vector>) driver(
	System.Func<double,vector,vector> F,/* the f from dy/dx=f(x,y) */
	(double,double) interval,    /* (initial-point,final-point) */
	vector yinit,                /* y(initial-point) */
	double h=0.125,              /* initial step-size */
	double acc=0.01,             /* absolute accuracy goal */
	double eps=0.01              /* relative accuracy goal */
    ){
        var (a,b)=interval; double x=a; vector y=yinit.copy();
        var xlist=new genlist<double>(); xlist.add(x);
        var ylist=new genlist<vector>(); ylist.add(y);
        do{
	        if(x>=b) return (xlist,ylist); /* job done */
	        if(x+h>b) h=b-x;               /* last step should end at b */
	        var (yh,δy) = rkstep23(F,x,y,h);
	        double tol = (acc+eps*yh.norm()) * Sqrt(h/(b-a));
	        double err = δy.norm();
	        if(err<=tol){ // accept step
		        x+=h; y=yh;
		        xlist.add(x);
		        ylist.add(y);
		    }
	        if(err>0) h *= Min( Pow(tol/err,0.25)*0.95 , 2); // readjust stepsize
	        else h*=2;
	    }while(true);
    }//driver
    static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Embedded rule Runge-Kutta ODE integrator");
        WriteLine("Function for debugging: u'=-u, initial value: (0,5)");
        WriteLine("The computed results using the implemented methods are plotted in first_ODE.gnuplot.svg along with the analytical function");
        WriteLine();
        WriteLine();
        System.Func<double,vector,vector> f = (x, y) => -y;
        genlist<double> xlist = driver(f, (0.0, 1.0), new vector(5.0)).Item1;
        genlist<vector> ylist = driver(f, (0.0, 1.0), new vector(5.0)).Item2;
        WriteLine("x    y");
        for(int i=0; i<xlist.size; i++){
            WriteLine($"{xlist[i]}  {ylist[i][0]}");
        }
        WriteLine();
        WriteLine("Oscillator with friction:");
        WriteLine("The computed theta(t) and omega(t) can be seen in osc_fric.gnuplot.svg");
        double b = 0.25;
        double c = 5.0;
        System.Func<double,vector,vector> omega_first = (x, theta_omega) => new vector(theta_omega[1], -b*theta_omega[1]-c*Sin(theta_omega[0]));
        genlist<double> tlist = driver(omega_first, (0.0, 10.0), new vector(PI-0.1, 0.0)).Item1;
        genlist<vector> theta_omega_vec = driver(omega_first, (0.0, 10.0), new vector(PI-0.1, 0.0)).Item2;
        WriteLine("t    theta(t)    omega(t)");
        WriteLine();
        WriteLine();
        for(int i=0; i<tlist.size; i++){
            WriteLine($"{tlist[i]}  {theta_omega_vec[i][0]}  {theta_omega_vec[i][1]}");
        }
        WriteLine();
        WriteLine("B. Relativistic precession of planetary orbit");
        WriteLine("The orbits are plotted in planet.gnuplot.svg");
        var data = new System.IO.StreamWriter("data_B.txt", append:true);
        double[] eps = {0.0, 0.0, 0.01} ;
        vector[] initials = {new vector(1.0, 0.0), new vector(1.0, -0.5), new vector(1.0, -0.5)};
        for(int i=0; i<eps.Length; i++){
            System.Func<double,vector,vector> u_double = (x, y_u) => new vector(y_u[1], 1-y_u[0]+eps[i]*y_u[0]*y_u[0]);
            genlist<double> phi = driver(u_double, (0.0, 32*PI), initials[i]).Item1;
            genlist<vector> u_ufirst = driver(u_double, (0.0, 32*PI), initials[i]).Item2;
            for(int t=0; t<phi.size; t++){
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
        System.Func<double,vector,vector> C1 = (x, y_yfirst) => new vector(y_yfirst[1], 2*x);
        genlist<double> xs = driver(C1, (0.0, 10.0), new vector(0, 0)).Item1;
        genlist<vector> f_ffirst = driver(C1, (0.0, 10.0), new vector(0, 0)).Item2;
        var data2 = new System.IO.StreamWriter("data_C.txt", append:true);
        for(int t=0; t<xs.size; t++){
                data2.WriteLine($"{xs[t]}   {f_ffirst[t][0]}    {f_ffirst[t][1]}");
            }
            data2.WriteLine();
            data2.WriteLine();
        data2.Close();
        /*
        vector stepsize = new vector(xs.size-1);
        stepsize[0]=xs[1];
        for(int i=1; i<stepsize.size; i++){
            stepsize[i] = xs[i]-xs[i-1];
        }
        vector double_step = new vector(xs.size-1);
        double_step[0]=xs[1];
        for(int i=1; i<double_step.size; i++){
            double_step[i]=2*double_step[i-1];
        }
        WriteLine($"    -Does stepsize double for every step? {double_step.approx(stepsize)}");
        stepsize.print();
        double_step.print(); 
        */
        return 0;
    }
}
