using static System.Console;
using static System.Math;
class main{
    static (double,double) plainmc(System.Func<vector,double> f,vector a,vector b,int N){
        int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        double sum=0,sum2=0;
	    var x=new vector(dim);
	    var rnd=new System.Random();
        for(int i=0;i<N;i++){
                for(int k=0;k<dim;k++)x[k]=a[k]+rnd.NextDouble()*(b[k]-a[k]);
                double fx=f(x); sum+=fx; sum2+=fx*fx;
                }
        double mean=sum/N, sigma=Sqrt(sum2/N-mean*mean);
        var result=(mean*V,sigma*V/Sqrt(N));
        return result;
    }
    static double corput(int n, int b=2){
        double q = 0;
        double bk = 1.0/b;
        while(n>0){
            q+= (n % b)*bk;
            n/=b;
            bk /= b;
        }
        return q;
    }
    static vector halton(int n, int d){
        int[] bases = new int[18] {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
        int maxd = bases.Length;
        vector x = new vector(d);
        for(int i=0; i<d; i++){
            x[i]=corput(n, bases[i]);
        }
        return x;
    } 

    static (double,double) haltonmc(System.Func<vector,double> f,vector a,vector b,int N){
        int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        double sum=0, sum2=0;
	    var x=new vector(dim);
        var x2=new vector(dim);
	    var rnd=new System.Random();
        for(int i=0;i<N;i++){
                vector y = halton(rnd.Next(), dim );
                vector y2 = halton(rnd.Next(), dim );
                for(int k=0;k<dim;k++)x[k]=a[k]+y[k]*(b[k]-a[k]);
                for(int k=0;k<dim;k++)x2[k]=a[k]+y2[k]*(b[k]-a[k]);
                double fx=f(x); sum+=fx;
                double fx2=f(x2); sum2+=fx2;
                }
        double mean=sum/N;
        var result=(mean*V,Abs(mean*V-V*sum2/N));
        return result;
    }

    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Plain Monte Carlo integration");
        //int[] Ns = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        System.Func<vector, double> f_A1 = (vec) => Exp(-(Pow(vec[0],2)+Pow(vec[1],2))); //Gaussian function, skal give 0.557746. Grænser 0,1 og 0,1
        System.Func<vector, double> f_A2 = (vec) => 1/(1-Cos(vec[0])*Cos(vec[1])*Cos(vec[2]))*Pow(PI,-3); // 
        //System.Func<vector, double> f_A3 = (vec) => Log(1+vec[0]*vec[1]);// Med singularitet, grænser 0,1 og 0,1. skal give 0.29
        WriteLine(" - Test implemented plain MC on ∫_0^1∫_0^1 exp(-(x^2+y^2)) dx dy=0.557746. The estimated and actual error as a funciton of N is plotted in MC.gnuplot.svg");
        var data1 = new System.IO.StreamWriter("data_A.txt", append:true);
        for(int i=1; i<100000; i+=250){
            data1.WriteLine($"{i}   {plainmc(f_A1, new vector(0, 0), new vector(1, 1), i ).Item2}   {Abs(plainmc(f_A1, new vector(0, 0), new vector(1, 1), i ).Item1-0.557746)}");
        }
        data1.WriteLine();
        data1.WriteLine();
        /*for(int i=0; i<1000; i++){
            data1.WriteLine($"{i}   {plainmc(f_A2, new vector(0, 0), new vector(1, 1), i ).Item2}   {Abs(plainmc(f_A2, new vector(0, 0), new vector(1, 1), i ).Item1-0.0)}");
        }
        data1.WriteLine();
        data1.WriteLine();
        for(int i=0; i<1000; i++){
            data1.WriteLine($"{i}   {plainmc(f_A3, new vector(0, 0), new vector(1, 1), i ).Item2}   {Abs(plainmc(f_A3, new vector(0, 0), new vector(1, 1), i ).Item1-0.29)}");
        }
        data1.Close();
        WriteLine(" - Test implemented plain MC on ∫_0^1∫_0^1 Sin(10x)*Cos(10y) dx dy=0.0. The actual error as a funciton of N is plotted in MC.gnuplot.svg");
        WriteLine(" - Test implemented plain MC on ∫_0^1∫_0^1 ln(1+xy) dx dy = 0.29. The actual error as a funciton of N is plotted in MC.gnuplot.svg"); */
        WriteLine($"    -Is the computed value the same as the table value for ∫0^π  dx/π ∫0^π  dy/π ∫0^π  dz/π [1-cos(x)cos(y)cos(z)]^-1?");
        WriteLine($"       Computed value ved N=500000: {plainmc(f_A2, new vector(0,0,0), new vector(PI, PI, PI), 500000).Item1}");
        WriteLine("        Table value = 1.3932039296856768591842462603255 ");
        WriteLine();
        WriteLine("B. Quasi-random sequences");
        WriteLine(" - Test implemented plain Halton MC on ∫_0^1∫_0^1 exp(-(x^2+y^2)) dx dy=0.557746. The estimated and actual error as a funciton of N is plotted in MC.gnuplot.svg");
        for(int i=1; i<100000; i+=250){
            data1.WriteLine($"{i}   {haltonmc(f_A1, new vector(0, 0), new vector(1, 1), i ).Item2}   {Abs(haltonmc(f_A1, new vector(0, 0), new vector(1, 1), i ).Item1-0.557746)}");
        }
        data1.Close();

        return 0;
    }
}
