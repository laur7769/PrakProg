using static System.Console;
using static System.Math;
using System.Collections.Generic;
public class ODE{
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
}
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
public class integration{
    public static double CC_var(System.Func<double,double> f, double a, double b, double accu=0.001, double epsi=0.001){
        System.Func<double,double> new_f = (theta) => f((a+b)/2+(b-a)/2*Cos(theta))*Sin(theta)*(b-a)/2;
        return integrate(new_f, 0.0, PI, acc:accu, eps:epsi );
    }
    public static double integrate(System.Func<double,double> f, double a, double b,
    double acc=0.001, double eps=0.001, double f2= double.NaN, double f3= double.NaN) // NaN indicates first call
    {
        if(double.IsPositiveInfinity(b) && double.IsNegativeInfinity(a)){
            System.Func<double,double> f_new1 = (t) => f(t/(1-Pow(t,2)))*(1+Pow(t,2))/Pow(1-Pow(t,2),2);
            return CC_var(f_new1,-1.0, 1.0, accu:acc, epsi:eps);
        }
        else if(double.IsPositiveInfinity(b)){
            System.Func<double,double> f_new2 = (t) => f(a+t/(1-t))*1/Pow(1-t,2);
            return CC_var(f_new2, 0.0, 1.0, accu:acc, epsi:eps);
        }
        else if(double.IsNegativeInfinity(a)){
            System.Func<double,double> f_new3 = (t) => f(b+t/(1+t))*1/Pow(1+t,2);
            return CC_var(f_new3, -1.0, 0.0, accu:acc, epsi:eps);
        }
        else {
        double h=b-a;
        if(double.IsNaN(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); } // first call, no points to reuse
        double f1=f(a+h/6), f4=f(a+5*h/6);
        double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
        double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
        double err = Abs(Q-q);
        if (err <= acc+eps*Abs(Q)) return Q;
        else return integrate(f,a,(a+b)/2,acc/Sqrt(2),eps,f1,f2)+ integrate(f,(a+b)/2,b,acc/Sqrt(2),eps,f3,f4);
        }
    }

    public static bool approx(double Q, double exact, double acc=1e-6, double eps=1e-6){
        double tol = acc+eps*Abs(exact);
            if(Abs(exact-Q)<tol)return true;
            return false;
    }

}
public class minimization{

    public static vector gradient(System.Func<vector,double>f, vector x){
    double fx = f(x);
    vector gf = new vector(x.size);
    for(int i=0; i<x.size; i++){
        double dxi = Abs(x[i])*Pow(2,-26);
        x[i] +=dxi;
        gf[i] = (f(x)-fx)/dxi;
        x[i]-=dxi;
    } 
    return gf;
}


public static matrix hessian(System.Func<vector,double>f, vector x){
    matrix H = new matrix(x.size, x.size);
    vector gfx = gradient(f, x);
    for(int j=0; j<x.size; j++){
        double dxj = Abs(x[j])*Pow(2,-26);
        if (dxj == 0) dxj = Pow(2, -26);
        x[j] += dxj;
        vector dgf = gradient(f, x) - gfx;
        for(int i=0; i<x.size; i++){
            H[i, j] = dgf[i]/dxj;
        }
        x[j]-=dxj;   
    } 
    return H;
}

    public static vector newton(
	System.Func<vector,double>f /* the function to find the root of */
	,vector start        /* the start point */
	,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	,vector δx=null      /* optional δx-vector for calculation of jacobian */
	){
        vector x=start.copy();
        double fx=f(x); 
        vector z;
        double fz;
        int counter = 0;
        do{ /* Newton's iterations */
            counter += 1;
            if(counter>=10000) break;
            vector g=gradient(f, x);
	        if(g.norm() < acc) break; /* job done */
            matrix H = hessian(f, x);
	        var QRH = QR.decomp(H);
	        vector Dx = QR.solve(QRH.Item1, QRH.Item2, -g); /* Newton's step */
            if (double.IsNaN(Dx.norm())) {
                //System.Console.WriteLine("dx = NaN. Attempting solution...");
                for (int i = 0; i < H.size1; i++) {
                    H[i, i] += 1e-10;
                }
                QRH = QR.decomp(H);
                Dx = QR.solve(QRH.Item1, QRH.Item2, -g);
                if (double.IsNaN(Dx.norm())) {/*System.Console.WriteLine("dx is still NaN...");*/}
            }

	        double λ=1;
            double λmin = 1.0/1024.0;
	        do{ /* linesearch */
                fz=f(x);
		        z=x+λ*Dx;
		        if( f(z) < fz ) break;
		        if( λ < λmin ) break;
		        λ/=2;
		    }while(true);
	        x=z; fx=fz;
	    }while(true);
        return x;
    }

}
public class ann{
    int n; /* number of hidden neurons */
    System.Func<double,double> f = x => x*Exp(-x*x); /* activation function */

    System.Func<double, double> df = x => Exp(-x*x)*(1-2*x*x);

    System.Func<double, double> ddf = x => 2*x*Exp(-x*x)*(2*x*x-1-2);

    System.Func<double, double> anti = x => -(1.0/2.0)*Exp(-x*x);
    vector p; /* network parameters */

    double alpha;

    double beta;



    public double deriv(double x){
        double sum = 0.0;
        for(int i=0; i<n; i++){
            sum += (1/p[3*i+1])*df((x-p[3*i])/p[3*i+1])*p[3*i+2];
        }
        return sum;
    } 

    public double sec_deriv(double x){
        double sum = 0.0;
        for(int i=0; i<n; i++){
            sum += (p[3*i+2]/Pow(p[3*i+1],2))*ddf((x-p[3*i])/p[3*i+1]);
        }
        return sum;
    }

    public double anti_deriv(double x){
        double sum = 0.0;
        double start = 0.0;
        for(int i=0; i<n; i++){
            sum += p[3*i+1]*p[3*i+2]*anti((x-p[3*i])/p[3*i+1]);
        }
        for(int i=0; i<n; i++){
            start += p[3*i+1]*p[3*i+2]*anti((-1-p[3*i])/p[3*i+1]);
        }
        return sum-start;
    }
    public void print(){
        p.print();
    }
    public ann(int n, bool dif = false){
        if(dif){
            this.n = n;
            p = new vector(3*n);
            var rnd = new System.Random();
            for (int i = 0; i < 3*n; i++) {
                p[i] =100.0*(rnd.NextDouble()-0.5);
            alpha = 1.0;
            beta = 1.0;
            }
        }    
        else{
            this.n = n;
            p = new vector(3*n);
            var rnd = new System.Random();
            for (int i = 0; i < 3*n; i++) {
                p[i] =100.0*(rnd.NextDouble()-0.5);
            }
        }

    }
    public double response(double x){
        double sum = 0;
        for(int i=0; i<n; i++){
            sum += f((x-p[3*i])/p[3*i+1])*p[3*i+2];
        }
        return sum;
    }
    public void train(vector x,vector y){
      System.Func<vector, double>C = (p)=>{
        System.Func<double, double> F = (xs) => {
            double sum = 0;
        for(int i=0; i<n; i++){
            sum += f((xs-p[3*i])/p[3*i+1])*p[3*i+2];
        }
        return sum;
        };
        double sum2 = 0;
        for(int k=0; k<x.size; k++){
            sum2 += Pow(F(x[k])-y[k],2);
        }
        return sum2;
      };
        do{
            var rnd = new System.Random();
                for (int i = 0; i < 3*n; i++) {
                this.p[i] =100.0*(rnd.NextDouble()-0.5);
         }
        }while(C(this.p)>200.0);  
        WriteLine("Before training:");
        WriteLine($"Cost={C(this.p)}");
        vector train_res = minimization.newton(C, this.p);
        this.p = train_res;
        WriteLine("After training:");
        WriteLine($"Cost={C(this.p)}");
   }

    public void train_dif(System.Func<double, double, double, double, double> phi, double a, double b, double c, double yc, double dyc){
        System.Func<vector, double> Cp = (p) => {
            System.Func<double, double> F = (xs) => {
                double sum = 0;
                for(int i=0; i<n; i++){
                    sum += f((xs-p[3*i])/p[3*i+1])*p[3*i+2];
                }
                return sum;
            };
        System.Func<double, double> sec_deriv = (xval) => {
            double sum = 0.0;
            for(int i=0; i<n; i++){
                sum += (p[3*i+2]/Pow(p[3*i+1],2))*ddf((xval-p[3*i])/p[3*i+1]);
            }
            return sum;
        };
        System.Func<double, double> deriv = (xvalue) => {
            double sum = 0.0;
            for(int i=0; i<n; i++){
                sum += (1/p[3*i+1])*df((xvalue-p[3*i])/p[3*i+1])*p[3*i+2];
            }
            return sum;
        };
        System.Func<double, double> integrand = (x) => Pow(phi(sec_deriv(x), deriv(x), F(x), x), 2);
        double term1 = integration.integrate(integrand, a, b);
        return term1 + alpha*Pow(F(c)-yc, 2) + beta*Pow(deriv(c)-dyc, 2);
        };
         do{
            var rnd = new System.Random();
                for (int i = 0; i < 3*n; i++) {
                this.p[i] =100.0*(rnd.NextDouble()-0.5);
         }
        }while(Cp(this.p)>5.0);
        WriteLine("Before training:");
        WriteLine($"Cost={Cp(this.p)}");
        vector train_res = minimization.newton(Cp, this.p, acc:0.0001);
        this.p = train_res;
        if(Cp(this.p)>0.01){train_res = minimization.newton(Cp, this.p);}
        this.p = train_res;
        WriteLine("After training:");
        WriteLine($"Cost={Cp(this.p)}");
        }

    }

class main{
//ODE routine
    
    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. ");
        WriteLine("Test neural network trained on g(x)=Cos(5*x-1)*Exp(-x*x):");
        ann neural_A = new ann(12);
        var data1 = new System.IO.StreamWriter("data_A.txt", append:true);
        vector train_x = new vector(400);
        vector train_y = new vector(400);
        for(int t=1; t<401; t++){
            double xval = -1+t*0.005;
            train_x[t-1] = xval;
            train_y[t-1] = Cos(5*xval-1)*Exp(-xval*xval);
            data1.WriteLine($"{train_x[t-1]}    {train_y[t-1]}");
        }
        data1.WriteLine();
        data1.WriteLine();
        neural_A.train(train_x, train_y);
        //neural_A.print();
        for(int t=1; t<101; t++){
            double xval = -1+t*0.02;
            data1.WriteLine($"{xval}    {neural_A.response(xval)}");
        }
        data1.Close();
        var data2 = new System.IO.StreamWriter("data_B.txt", append:true);
        for(int t=1; t<101; t++){
            double xval = -1+t*0.02;
            data2.WriteLine($"{xval}    {neural_A.deriv(xval)}");
        }
        data2.WriteLine();
        data2.WriteLine();
        for(int t=1; t<101; t++){
            double xval = -1+t*0.02;
            data2.WriteLine($"{xval}    {neural_A.sec_deriv(xval)}");
        }
        data2.WriteLine();
        data2.WriteLine();
        for(int t=1; t<101; t++){
            double xval = -1+t*0.02;
            data2.WriteLine($"{xval}    {neural_A.anti_deriv(xval)}");
        }
        data2.Close();
    WriteLine("The response of the trained neural network with n=12 and trained on 400 data points can be seen plotted against g(x) in A.gnuplot.svg.");
    WriteLine();
    WriteLine("B.  ");
    WriteLine("The neural network with activation function f=x*exp(-x^2) has been modified to include methods giving the first and second derivative as well as the anti-derivative.");
    WriteLine("The first and scond derivative as well as the anti-derivative for a neural network trained on g(x) can be seen plotteed against the analytical functions in B.gnuplot.svg.");
    WriteLine();
    WriteLine("C. ");
    WriteLine("Training the neural network to approximate the solution to d^2y/dx^2+1/2*dy/dx + y = 0");
    WriteLine(" a= 0, b=10, c=0, dy/dx(c)=1, d^2y/dx^2(c)=0");
    WriteLine("The neural network has n=12 and has the same activation function as before");
    ann neural_C = new ann(14, dif:true);
    System.Func<double, double, double, double, double> phi = (double d2y, double dy, double y, double x) => {
            double m = 1;
            double b = 0.5;
            double k = 1;
            return m*d2y + b*dy + k*y;
        };
    neural_C.train_dif(phi, 0, 10, 0, 1, 0);
    var data3 = new System.IO.StreamWriter("data_C.txt", append:true);
    System.Func<double, vector, vector> dt = (val, t) => new vector(t[1], -0.5*t[1]-t[0]);
    genlist<double> xlist = ODE.driver(dt, (0,10), new vector(1, 0)).Item1;
    genlist<vector> res = ODE.driver(dt, (0,10), new vector(1, 0)).Item2;
    for(int t=1; t<xlist.size; t++){
            data3.WriteLine($"{xlist[t]}    {res[t][0]} ");
        }
    data3.WriteLine();
    data3.WriteLine();
    for(int t=1; t<xlist.size; t++){
            data3.WriteLine($"{xlist[t]}    {neural_C.response(xlist[t])} ");
        }
    data3.Close();
    WriteLine("The neural network respone can be seen plotted againd the numerical solution to the differential equation in C.gnuplot.svg.");
    WriteLine("The numerical solution was computed using an ODE routine using a 2,3 stepper.");

    return 0;
    }

}
