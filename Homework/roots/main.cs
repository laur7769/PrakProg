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
//ODE routine
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
//End of ODE routine
    public static matrix jacobian
    (System.Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
	    if(dx == null) dx = x.map(xi => Abs(xi)*Pow(2,-26));
	    if(fx == null) fx = f(x);
	    matrix J=new matrix(x.size);
	    for(int j=0;j < x.size;j++){
		    x[j]+=dx[j];
		    vector df=f(x)-fx;
		    for(int i=0;i < x.size;i++) J[i,j]=df[i]/dx[j];
		    x[j]-=dx[j];
		}
	    return J;
    }
    public static vector newton(
	System.Func<vector,vector>f /* the function to find the root of */
	,vector start        /* the start point */
	,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	,vector δx=null      /* optional δx-vector for calculation of jacobian */
	){
        vector x=start.copy();
        vector fx=f(x); 
        vector z;
        vector fz;
        do{ /* Newton's iterations */
	        if(fx.norm() < acc) break; /* job done */
	        matrix J=jacobian(f,x,fx,δx);
	        var QRJ = QR.decomp(J);
	        vector Dx = QR.solve(QRJ.Item1, QRJ.Item2, -fx); /* Newton's step */
	        double λ=1;
            double λmin = 1.0/128.0;
	        do{ /* linesearch */
		        z=x+λ*Dx;
		        fz=f(z);
		        if( fz.norm() < (1-λ/2)*fx.norm() ) break;
		        if( λ < λmin ) break;
		        λ/=2;
		    }while(true);
	        x=z; fx=fz;
	    }while(true);
        return x;
    }

    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Newton's method with numerical Jacobian and back-tracking line-search");
        WriteLine("     - Test root finding routine on f(x)=(x-2)(x+3):");
        System.Func<vector,vector>f1 = (x) => new vector((x[0]-2)*(x[0]+3));
        WriteLine($"       Computed root with start 1.0 and default accuracy: {newton(f1, new vector(1.0) )[0]}");
        WriteLine($"       Computed root with start -1.0 and default accuracy: {newton(f1, new vector(-1.0) )[0]}");
        WriteLine($"       Actual roots: 2 and -3");
        WriteLine("     - Test root finding routine on f(x,y)=[x^2+y^2-4, x-y]:");
        System.Func<vector,vector>f2 = (vec) => new vector(vec[0]*vec[0]+vec[1]*vec[1]-4, vec[0]-vec[1]);
        vector res1 =newton(f2, new vector(1.0, 1.0) );
        vector res2 =newton(f2, new vector(-1.0, -1.0) );
        WriteLine($"       Computed root with start (1.0, 1.0) and default accuracy: {(res1[0],res1[1])}");
        WriteLine($"       Computed root with start (-1.0, -1.0) and default accuracy: {(res2[0], res2[1])}");
        WriteLine($"       Actual roots: (√2, √2) and (-√2, -√2)");
        WriteLine("     - Find extremum(s) of the Rosenbrock's valley function: f(x,y) = (1-x)^2+100(y-x^2)^2");
        WriteLine("        Analytical gradient: ∂f/∂x=2x-2+400x^3-400xy,   ∂f/∂y=200y-200x^2");
        WriteLine("        It is known that the Rosenbrock's valley function only has one extremum, a global minimum at (1,1)");
        System.Func<vector,vector>Rosen_grad = (vec) => new vector(2*vec[0]-2+400*Pow(vec[0],3)-400*vec[0]*vec[1], 200*vec[1]-200*vec[0]*vec[0]);
        vector ros_res = newton(Rosen_grad, new vector(1.0, 1.0), acc:1e-4);
        WriteLine($"        Computed root with start (7.0, -2.0) and acc=1e-4: {(ros_res[0],ros_res[1])}");
        WriteLine("     - Find minimum(s) ofthe Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
        WriteLine("        Analytical gradient: ∂f/∂x=2x-2+400x^3-400xy,   ∂f/∂y=200y-200x^2");
        System.Func<vector,vector>Himmel_grad = (vec) => new vector(4*Pow(vec[0],3)+4*vec[0]*vec[1]-42*vec[0]+2*vec[1]*vec[1]-14, 2*vec[0]*vec[0]-22+4*vec[0]*vec[1]+4*Pow(vec[1],3)-26*vec[1]);
        WriteLine("        Analytical gradient: 4x^2+4xy-42x+2y^2-14,   ∂f/∂y=2x^2-22+4yx+4y^3-26y");
        WriteLine("        It is known that Himmelblau's funciton has 4 minimums: (-2.805118, 3.131312), (3.0, 2.0), (-3.779310, -3.283186) and (3.584428, -1.848126)") ;
    var points = new List<(double, double)> {(-10, -10), (-10, 10), (10, -10), (10, 10)};
    vector xx = new vector(2);
        foreach ((double p1, double p2) point in points) {
            xx[0] = point.p1; xx[1] = point.p2;
            vector solution = newton(Himmel_grad, xx, acc:1e-10);
            WriteLine($"        Computed root with start {(point.p1, point.p2)} and acc:1e-10: {(solution[0], solution[1])}");
        }
    
    WriteLine("B. Bound states od hydrogen atom with shooting method for boundary value problems");
    double rmin=0.01;
    double rmax=8.0;
    vector y_initial = new vector(rmin-rmin*rmin, 1.0-2.0*rmin);
    System.Func<vector, vector> M_E = (E) => {
        System.Func<double,vector,vector> hydrogen = (r, y) => new vector(y[1], 2.0*(-1.0/r * y[0] - E[0]*y[0]));
        genlist<double> rs = driver(hydrogen, (rmin, rmax), y_initial).Item1;
        genlist<vector> fs = driver(hydrogen, (rmin, rmax), y_initial).Item2;
        return new vector(fs[rs.size-1][0]);
    };
    WriteLine($"E_0 with rmax=8 and rmin=0.01 and start -1.0: {newton(M_E, new vector(-1.0))[0]}"); 
    // finder bølgefunktionen for E_0
    vector E_0 = newton(M_E, new vector(-1.0));
    System.Func<double,vector,vector> hyd = (r, y) => new vector(y[1], 2.0*(-1.0/r * y[0] - E_0[0]*y[0]));
    genlist<double> r_val = driver(hyd, (rmin, rmax), y_initial).Item1;
    genlist<vector> f_val = driver(hyd, (rmin, rmax), y_initial).Item2;  
    var data1 = new System.IO.StreamWriter("data_B.txt", append:true);
    for(int i=0; i<r_val.size; i++){
        data1.WriteLine($"{r_val[i]}   {f_val[i][0]}");
    }
    data1.WriteLine();
    data1.WriteLine();
    WriteLine("The computed f_0 calculated using the computed E_0 above, can be seen plotted against the exact f_0 in A.gnuplot.svg");
    WriteLine("The exact E_0=-1/2 and the exact f_0=r*exp(-r). This is the case as the Schrödinger equation can be rewritten as: ");
    WriteLine("f''=-2f(E+1/r)");
    WriteLine("As f_0''=-2exp(-r)+r*exp(-r) and -2*f_0(E_0+1/r)=-2*r*exp(-r)(-1/2+1/r)=r*exp(-r)-2exp(-r).");
    WriteLine("Therefore, E_0 and f_0 are indeed solutions to the Schrödinger equation");
    WriteLine();
    WriteLine("Test of convergence:");
    WriteLine("Convergence of E_0 with respect to change in rmax with rmin=0.01 can be seen in rmax.gnuplot.svg");
    
    for(int r=1;r<12; r++){
        rmax=r;
        data1.WriteLine($"{rmax} {newton(M_E, new vector(-1.0))[0]}");
    }
    WriteLine("The plot of E_0 as a function of rmin with rmax=8 can be seen in rmin.gnuplot.svg");
    
    data1.WriteLine();
    data1.WriteLine();
    rmax = 8;
    double[] rmins = { 0.0001, 0.0025, 0.0005, 0.0075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5};
    for(int i=0; i<rmins.Length; i++){
        rmin=rmins[i];
        data1.WriteLine($"{rmin} {newton(M_E, new vector(-1.0))[0]}");
    }
    data1.WriteLine();
    data1.WriteLine();
    double ac = 0.01;
    double ep = 0.01;
    System.Func<vector, vector> M_E2 = (E) => {
        System.Func<double,vector,vector> hydrogen = (r, y) => new vector(y[1], 2.0*(-1.0/r * y[0] - E[0]*y[0]));
        genlist<double> rs = driver(hydrogen, (rmin, rmax), y_initial, acc:ac, eps:ep).Item1;
        genlist<vector> fs = driver(hydrogen, (rmin, rmax), y_initial, acc:ac, eps:ep).Item2;
        return new vector(fs[rs.size-1][0]);
    };
    rmin = 0.01;
    rmax = 8;
    double[] accs = {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5};
    double[] epss = {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5};
    for(int i=0; i<accs.Length; i++){
        ac=accs[i];
        data1.WriteLine($"{ac} {newton(M_E2, new vector(-1.0))[0]}");
    }
    ac = 0.01;
    data1.WriteLine();
    data1.WriteLine();
    for(int i=0; i<epss.Length; i++){
        ep=epss[i];
        data1.WriteLine($"{ep} {newton(M_E2, new vector(-1.0))[0]}");
    }
    data1.Close();
    WriteLine("The plot of E_0 as a function of accuracy for the ODE driver can be seen in acc.gnuplot.svg");
    WriteLine("The plot of E_0 as a function of epsilon for the ODE driver can be seen in eps.gnuplot.svg");


    return 0;
    }

}
