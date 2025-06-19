using static System.Console;
using static System.Math;
class main{
    static (double, double) integrate(System.Func<double,double> f, double a, double b,
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
            if (err <= acc + eps * Abs(Q)) return (Q, err);
            else
            {
                var (Q1, err1) = integrate(f, a, (a + b) / 2, acc / Sqrt(2), eps, f1, f2);
                var (Q2, err2) = integrate(f, (a + b) / 2, b, acc / Sqrt(2), eps, f3, f4);
                return (Q1 + Q2, Sqrt(err1 * err1 + err2 * err2));
            }
        }
    }
    static double erf(double z, double accuracy=0.001, double epsilon=0.001){
        if(z<0){
            return -erf(-z);
        }
        else if(0<=z && z<=1){
            System.Func<double,double> F = (x) => Exp(-Pow(x,2));
            return 2.0/Sqrt(PI)*integrate(F, 0, z, acc:accuracy, eps:epsilon ).Item1;
        }
        else{
            System.Func<double,double> F = (t) => Exp(-Pow(z+(1-t)/t, 2))/t/t;
            return 1-2.0/Sqrt(PI)*integrate(F, 0, 1, acc:accuracy, eps:epsilon).Item1;
        }
    }

    public static bool approx(double Q, double exact, double acc=1e-6, double eps=1e-6){
        double tol = acc+eps*Abs(exact);
	    if(Abs(exact-Q)<tol)return true;
	    return false;
    }

    public static (double, double) CC_var(System.Func<double,double> f, double a, double b, double accu=0.001, double epsi=0.001){
        System.Func<double,double> new_f = (theta) => f((a+b)/2+(b-a)/2*Cos(theta))*Sin(theta)*(b-a)/2;
        return integrate(new_f, 0.0, PI, acc:accu, eps:epsi );
    }
    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Recursive open 4-point adaptive integrator");
        WriteLine("Test iimplemented integrator with certain known definite integrals:");
        System.Func<double,double> f1 = (x) => Sqrt(x);
        WriteLine($"    * ∫_0^1 dx √(x) = 2/3. Is integrate(√(x), 0, 1)=2/3 with acc=0.001 and eps=0.001? {approx(integrate(f1, 0, 1).Item1, 2.0/3.0, acc:0.001, eps:0.001)}");
        System.Func<double,double> f2 = (x) => 1.0/Sqrt(x);
        WriteLine($"    * ∫_0^1 dx 1/√(x) = 2. Is integrate(1/√(x), 0, 1)=2 with acc=0.001 and eps=0.001? {approx(integrate(f2, 0, 1).Item1, 2.0, acc:0.001, eps:0.001)}");
        System.Func<double,double> f3 = (x) => 4.0*Sqrt(1-Pow(x,2));
        WriteLine($"    * ∫_0^1 dx 4*√(1-x^2) = π. Is integrate(4*√(1-x^2), 0, 1)=π with acc=0.001 and eps=0.001? {approx(integrate(f3, 0, 1).Item1, PI, acc:0.001, eps:0.001)}");
        System.Func<double,double> f4 = (x) => Log(x)/Sqrt(x);
        WriteLine($"    * ∫_0^1 dx ln(x)/√(x) = -4. Is integrate(ln(x)/√(x), 0, 1)=-4 with acc=0.001 and eps=0.001? {approx(integrate(f4, 0, 1).Item1, -4.0, acc:0.001, eps:0.001)}");
        var data1 = new System.IO.StreamWriter("data_A.txt", append:true);
        double[] erf_tab = {0.0, 0.520499878, 0.842700793, 0.966105146, 0.995322265, 0.999593048, 0.999977910, 0.999999257};
        double[] zs = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
        for(int u=0; u<zs.Length; u++){
            data1.WriteLine($"{zs[u]}   {erf_tab[u]}");
        }
        data1.WriteLine();
        data1.WriteLine();
        for(int i=-50; i<50; i++){
            data1.WriteLine($"{i/10.0} {erf(i/10.0)}");
        }
        data1.WriteLine();
        data1.WriteLine();
        double[] accs = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
        for(int a=0; a<accs.Length; a++){
            data1.WriteLine($"{accs[a]}   {Abs(erf(1, accuracy:accs[a], epsilon:0)-erf_tab[2])}");
        }
        data1.WriteLine();
        data1.WriteLine();
        data1.Close();
        WriteLine("The plot of the implemented error function as well as tabulated values of the error function can be seen in erf.gnuplot.svg ");
        WriteLine("The plot of the absolute difference between the tabulated and compted value of erf(1) as a funciton of accuracy, can be seen in acc.gnuplot.svg");
        WriteLine();
        WriteLine("B. Variable transformation quadratures");
        WriteLine("For all below calculations acc=0.001 and eps=0.001");
        int ncalls = 0;
        System.Func<double,double> f_b1 = z => {ncalls++;return 1.0/Sqrt(z);};
        WriteLine("Integral to evaluate: ∫_0^1 dx 1/√(x) = 2");
        WriteLine($"    - Computed value using implemented Clenshaw-Curtis: {CC_var(f_b1, 0.0, 1.0).Item1}, integrand evaluations: {ncalls}");
        ncalls=0;
        var result1 = integrate(f_b1, 0.0, 1.0);
        WriteLine($"    - Computed value using only 4-point open adaptive quadrature: {result1.Item1}, integrand evaluations: {ncalls} ");
        WriteLine($"    - Computed value using scipy's quad: 1.9999999999999993, integrand evaluations: 231");
        WriteLine("Integral to evaluate: ∫_0^1 dx ln(x)/√(x) = -4");
        ncalls=0;
        System.Func<double,double> f_b2 = z => {ncalls++;return Log(z)/Sqrt(z);};
        WriteLine($"    - Computed value using implemented Clenshaw-Curtis: {CC_var(f_b2, 0.0, 1.0).Item1}, integrand evaluations: {ncalls}");
        ncalls=0;
        var result2 = integrate(f_b2, 0.0, 1.0);
        WriteLine($"    - Computed value using only 4-point open adaptive quadrature: {result2.Item1}, integrand evaluations: {ncalls} ");
        WriteLine($"    - Computed value using scipy's quad: -3.99999999999998273, integrand evaluations: 315");
        ncalls=0;
        WriteLine("Integral to evaluate: ∫_0^∞ dx exp(-x)*x^2 = 2");
        System.Func<double,double> f_b3 = z => {ncalls++;return Exp(-z)*Pow(z,2);};
        var result3 = integrate(f_b3, 0.0, double.PositiveInfinity);
        WriteLine($"    - Computed value using implemented integrator: {result3.Item1}, integrand evaluations: {ncalls}");
        WriteLine($"    - Computed value using scipy's quad: 2.0, integrand evaluations: 165");
        WriteLine();
        WriteLine("C. Error estimate");
        WriteLine(" The implemented integrator has been modified to also return the estimated error.");
        WriteLine(" The same three integrals from exercise B. are evaluated with the same acc and eps.");
        WriteLine("Integral to evaluate: ∫_0^1 dx 1/√(x) = 2");
        WriteLine($"     - Estimated error: {result1.Item2}");
        WriteLine($"     - Actual error: {Abs(result1.Item1 - 2)}");
        WriteLine("Integral to evaluate: ∫_0^1 dx ln(x)/√(x) = -4");
        WriteLine($"     - Estimated error: {result2.Item2}");
        WriteLine($"     - Actual error: {Abs(result2.Item1 + 4)}");
        WriteLine(" Integral to evaluate: ∫_0^∞ dx exp(-x)*x^2 = 2");
        WriteLine($"     - Estimated error: {result3.Item2}");
        WriteLine($"     - Actual error: {Abs(result3.Item1 - 2)}");
        WriteLine(" As seen above, the estimated error is slightly higher than the actual error at this acc and eps level.");

        return 0;

    }
}
