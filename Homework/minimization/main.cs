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
    public static (vector,int) newton(
	System.Func<vector,double>f /* the function to find the root of */
	,vector start        /* the start point */
	,double acc=1e-3     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	,vector δx=null      /* optional δx-vector for calculation of jacobian */
    , bool central_ax=false
	){
        vector x=start.copy();
        int counter = 0;
        do{ /* Newton's iterations */
            counter += 1;
            if(counter>=2000) break;
            vector g = gradient(f, x);
            matrix H = hessian(f,x);
            if(central_ax){g = central(f,x).Item1; H = central(f,x).Item2;}
	        if(g.norm() < acc) break; /* job done */
	        var QRH = QR.decomp(H);
	        vector Dx = QR.solve(QRH.Item1, QRH.Item2, -g); /* Newton's step */;
	        double λ=1;
            double λmin = 1.0/128.0;
            double fx=f(x); 
	        do{ /* linesearch */
		        if( f(x+λ*Dx) < fx)  break;
		        λ/=2;
		    }while(λ > λmin);
            x+=λ*Dx;
	    }while(true);
        return (x,counter);
    }

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

    public static (vector, matrix) central(System.Func<vector, double> f, vector x)
    {
        vector gf = new vector(x.size);
        matrix H = new matrix(x.size, x.size);
        double fx = f(x);
        for (int i = 0; i < x.size; i++)
        {
            double dxi = Abs(x[i]) * Pow(2, -13);
            x[i] += dxi;
            double fx_plus = f(x);
            x[i] -= 2 * dxi;
            double fx_minus = f(x);
            x[i] += dxi;
            gf[i] = (fx_plus - fx_minus) / (2 * dxi);
            for (int k = 0; k < x.size; k++)
            {
                if (k == i) { H[k, i] = (fx_plus - 2 * fx + fx_minus) / Pow(dxi, 2); }
                else
                {
                    H[i, k] = 0;
                    double dxk = Abs(x[k]) * Pow(2, -13);
                    x[k] += dxk;
                    x[i] += dxi;
                    double f1 = f(x);
                    x[i] -= 2 * dxi;
                    double f2 = f(x);
                    x[k] -= 2 * dxk;
                    double f4 = f(x);
                    x[i] += 2 * dxi;
                    double f3 = f(x);
                    H[i, k] = (f1 - f2 - f3 + f4) / (4 * dxk * dxi);
                    x[i] -= dxi;
                    x[k] += dxk;
                }
            }
        }
        return (gf, H);
    }
    public static bool mini(vector sol, double accu)
        {
            vector min1 = new vector(-2.805118, 3.131312);
            vector min2 = new vector(3.0, 2.0);
            vector min3 = new vector(-3.779310, -3.283186);
            vector min4 = new vector(3.584428, -1.848126);
            if (sol.approx(min1, acc: accu) || sol.approx(min2, acc: accu) || sol.approx(min3, acc: accu) || sol.approx(min4, acc: accu)) return true;
            else return false;
        }
    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Newton's method with numerical gradient, numerical Hessian and back-tracking line-search");
        WriteLine("     - Find minimum of the Rosenbrock's valley function: f(x,y) = (1-x)^2+100(y-x^2)^2");
        WriteLine("        It is known that the Rosenbrock's valley function only has one extremum, a global minimum at (1,1)");
        System.Func<vector,double>Rosen_grad = (vec) => Pow(1-vec[0],2)+100*Pow(vec[1]-vec[0]*vec[0],2);
        vector ros_res = newton(Rosen_grad, new vector(7.0, -2.0), acc:1e-4).Item1;
        WriteLine($"        Computed minimum with start (7.0, -2.0) and acc=1e-4: {(ros_res[0],ros_res[1])}");
        WriteLine($"        Is the comuted root equal to the global minimum within the given accuracy? {ros_res.approx(new vector(1.0, 1.0), acc: 1e-4)}");
        WriteLine($"        Number of newton steps: {newton(Rosen_grad, new vector(7.0, -2.0), acc:1e-4).Item2}");
        vector ros_res2 = newton(Rosen_grad, new vector(7.0, -2.0), acc:1e-4, central_ax:true).Item1;
        WriteLine($"        Computed root with start (7.0, 2.0) and acc=1e-4, with central difference aproximation: {(ros_res2[0],ros_res2[1])}");
        WriteLine($"        Is the computed minimum equal to the known minimum within the given accuracy? {ros_res2.approx(new vector(1.0, 1.0), acc: 1e-4)}");
        WriteLine($"        Number of newton steps: {newton(Rosen_grad, new vector(7.0, -2.0), acc:1e-4, central_ax:true).Item2}");
        WriteLine("     - Find minimum(s) ofthe Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
        System.Func<vector,double>Himmel_grad = (vec) => Pow(vec[0]*vec[0]+vec[1]-11,2)+Pow(vec[0]+vec[1]*vec[1]-7,2);
        WriteLine("        It is known that Himmelblau's funciton has 4 minimums: (-2.805118, 3.131312), (3.0, 2.0), (-3.779310, -3.283186) and (3.584428, -1.848126)") ;

    var points = new List<(double, double)> {(10.0, 10.0), (-10.0, -10.0), (-10.0, 10.0), (10.0, -10.0)};
    vector xx = new vector(2);
        foreach ((double p1, double p2) point in points)
        {
            xx[0] = point.p1; xx[1] = point.p2;
            vector solution = newton(Himmel_grad, xx, acc: 1e-6).Item1;
            vector solution2 = newton(Himmel_grad, xx, acc: 1e-6, central_ax: true).Item1;
            WriteLine($"        Computed root with start {(point.p1, point.p2)} and acc:1e-6: {(solution[0], solution[1])}");
            WriteLine($"        Is computed minimum equal to one of the known minimums within the given accuracy? {mini(solution, 1e-6)}");
            WriteLine($"        Number of newton steps: {newton(Himmel_grad, xx, acc: 1e-6).Item2}");
            WriteLine();
            WriteLine($"        Computed root with start {(point.p1, point.p2)} and acc:1e-6 and central difference approximation: {(solution2[0], solution2[1])}");
            WriteLine($"        Is computed minimum equal to one of the known minimums within given accuracy? {mini(solution2, 1e-6)}");
            WriteLine($"        Number of newton steps: {newton(Himmel_grad, xx, acc: 1e-6, central_ax: true).Item2}");
            WriteLine();
            WriteLine();
        }
    WriteLine();
    WriteLine("B. Higgs bosson discovery");
    var energy = new genlist<double>();
    var signal = new genlist<double>();
    var error  = new genlist<double>();
    var separators = new char[] {' ','\t'};
    var options = System.StringSplitOptions.RemoveEmptyEntries;
    do{
        string line=In.ReadLine();
        if(line==null)break;
        string[] words=line.Split(separators,options);
        energy.add(double.Parse(words[0]));
        signal.add(double.Parse(words[1]));
        error.add(double.Parse(words[2]));
    }while(true);
    System.Func<vector, double>D =(konst) =>{
        double m = konst[0];
        double Gamma = konst[1];
        double A = konst[2];
        //System.Func<double, double>F = (E) => A/(Pow(E-m,2)+Gamma*Gamma/4.0);
        double sum = 0;
        for(int i=0; i<energy.size; i++){
            sum += Pow((A/(Pow(energy[i]-m,2)+Gamma*Gamma/4.0)-signal[i])/(error[i]),2);
        }
        return sum;
        
    };
    WriteLine("     - Computed constants for the Breit-Wigner fucntion with initial guess (126, 2, 10) and and accuracy of 1e-10: ");
    vector Higgs_res = newton(D, new vector(126.0, 2.0, 10.0), acc:1e-10).Item1;
    WriteLine($"         m = {Higgs_res[0]},   Gamma={Higgs_res[1]},  A={Higgs_res[2]}");
    WriteLine($"     - newton steps: {newton(D, new vector(126.0, 2.0, 10.0), acc:1e-10).Item2}");
    var writer = new System.IO.StreamWriter("fit.txt", append:true);
    for(double i=100.0; i<160.0; i+=0.1){
        writer.WriteLine($"{i}  {Higgs_res[2]/(Pow(i-Higgs_res[0],2)+Higgs_res[1]*Higgs_res[1]/4.0)}");
    }
    writer.WriteLine();
    writer.WriteLine();
    WriteLine("     The plot of the experimental data and the fit with the parameters above, can be seen in Higgs.gnuplot.svg");
    WriteLine($"     Even though the width should be around 4.1 MeV, which is magnitudes lower than {Higgs_res[1]} GeV, the computed Gamma is reasonable.");
    WriteLine("     This is the case, as its experimentally very hard to obtain the pysical width, as it is so small.");

    WriteLine();
    WriteLine("C. Central instead of forward finite difference approximation for the derivatives");
    vector Higgs_res2 = newton(D, new vector(126.0, 2.0, 10.0), acc:1e-10, central_ax:true).Item1;
    WriteLine("     - Computed constants for the Breit-Wigner function using the central difference approximation with initial guess (126, 2, 10): ");
    WriteLine($"         m = {Higgs_res2[0]},   Gamma={Higgs_res2[1]},  A={Higgs_res2[2]}");
    WriteLine($"     - newton steps: {newton(D, new vector(126.0, 2.0, 10.0), acc:1e-10, central_ax:true).Item2}");
    for(double i=100.0; i<160.0; i+=0.1){
        writer.WriteLine($"{i}  {Higgs_res2[2]/(Pow(i-Higgs_res2[0],2)+Higgs_res2[1]*Higgs_res2[1]/4.0)}");
    }
    writer.Close();
    WriteLine("     It is seen how rougly the same minimum is found, but the central difference approximation uses significaltly fewer netwon steps.");
    WriteLine("      This is however not the case when the central difference approximation is used to minimize Himmelblau's and Rosenbrock's valley functions.");
    return 0;
    }

}
