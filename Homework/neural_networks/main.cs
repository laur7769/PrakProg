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
    public ann(int n){
        this.n = n;
        p = new vector(3*n);
        var rnd = new System.Random();
         for (int i = 0; i < 3*n; i++) {
            p[i] =100.0*(rnd.NextDouble()-0.5);
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
        WriteLine("Before training:");
        WriteLine($"Cost={C(this.p)}");
        do{
            var rnd = new System.Random();
                for (int i = 0; i < 3*n; i++) {
                this.p[i] =100.0*(rnd.NextDouble()-0.5);
         }
        }while(C(this.p)>200.0);  
        vector train_res = minimization.newton(C, this.p);
        this.p = train_res;
        WriteLine($"Cost={C(this.p)}");
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
        WriteLine("After training:");
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

    return 0;
    }

}
