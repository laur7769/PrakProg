using static System.Console;
using static System.Math;
class main{
    static (vector, matrix) lsfit(System.Func<double,double>[] fs, vector x, vector y, vector dy){
            int n = x.size;
            int m = fs.Length;
            matrix A = new matrix(n , m);
            for(int i=0; i<n; i++){
                for(int k=0; k<m; k++){
                    A.set(i, k, fs[k](x[i])/dy[i]);
                }
            }
            vector b = new vector(n);
            for(int i=0; i<n; i++){
                b[i] = y[i]/dy[i];
            }
            matrix Q = QR.decomp(A).Item1;
            matrix R = QR.decomp(A).Item2;
            vector c = QR.solve(Q, R, b);
            matrix cov = QR.inverse(QR.decomp(A.T*A).Item1, QR.decomp(A.T*A).Item2);
            return (c, cov);

        }
    static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        System.Random rand = new System.Random();
        WriteLine("A. Ordinary least-squares fit by QR-decomposition");
        WriteLine("The QR decomposition is checked for a 6x3 matrix with random integer entries between 0 and 20:");
        matrix A = new matrix(6, 3);
        for(int i=0; i<A.size1; i++){
            for(int k=0; k<A.size2; k++){
                A.set(i, k, rand.Next(0, 20));
            }
        }
        matrix Q = QR.decomp(A).Item1;
        matrix R = QR.decomp(A).Item2;
        WriteLine($"Is R upper triangular? {QR.upper_triangular(R)}");
        matrix ID = new matrix(Q.size2, Q.size2); 
        for(int i=0;i<ID.size1;i++){
		ID.set(i,i,1);
		    for(int j=i+1;j<ID.size2;j++){
			ID.set(i,j,0);
            ID.set(j,i,0);
		    }
	    }
        WriteLine($"Is Q^TQ=I? {ID.approx(Q.T*Q)}");
        WriteLine($"Is QR=A? {A.approx(Q*R)}");
        WriteLine();
        WriteLine("Given data:");
        WriteLine("x    log(y)");
        WriteLine();
        WriteLine();
        
        double[] t = {1.0,  2.0,  3.0, 4.0, 6.0, 9.0,   10.0,  13.0,  15.0};
        double [] lnact = {Log(117.0), Log(100.0), Log(88.0), Log(72.0), Log(53.0), Log(29.5), Log(25.2), Log(15.2), Log(11.1)};
        vector x = new vector(t);
        vector y = new vector(lnact);
        double[] delta = {6.0, 5.0, 4.0, 4.0, 4.0, 3.0, 3.0, 2.0, 2.0};
        vector dy = new vector(9);
        for(int i=0; i<dy.size; i++){
            dy[i] = delta[i]/Exp(y[i]);
        }
        var fs = new System.Func<double,double>[] {z => 1.0 , z => -z };
        vector ck = lsfit(fs, x, y, dy).Item1;
        matrix cov = lsfit(fs, x, y, dy).Item2;
        for(int i=0; i<x.size; i++){
            WriteLine($"{x[i]} {y[i]}  {dy[i]}");
        }
        WriteLine();
        WriteLine();
        WriteLine("The function with the calculated least-squares coefficients:");
        WriteLine($"f1(x)={ck[0]}-{ck[1]}*x");
        WriteLine();
        double T_half = Log(2)/ck[1];
        WriteLine("The the data as well as the least squares fit can be seen in fit.gnuplot.svg");
        WriteLine($"Half-life determined from fit= {T_half}");
        WriteLine($"Half-life table value=3.6313(14)");
        WriteLine();
        WriteLine("B. Uncertainties of the fitting coefficients");
        double dT_half = Sqrt(Pow(-Log(2)*Pow(ck[1],-2.0)*Sqrt(cov[1,1]),2));
        WriteLine($"The uncertainty of the half life value from the given data is ∆f=df/λ*∆λ=√(-ln(2)*λ^(-2)*∆λ)^2={dT_half}");
        WriteLine($"Is the value for the half-life based in the given data agree with the table value within the estimated uncertainty? {matrix.approx(T_half, 3.6313, acc:dT_half)}");


        return 0;
    }
}
