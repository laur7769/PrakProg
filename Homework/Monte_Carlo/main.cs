using static System.Console;
using static System.Math;
class main {
    static (double, double) plainmc(System.Func<vector, double> f, vector a, vector b, int N) {
        int dim = a.size; double V = 1; for (int i = 0; i < dim; i++) V *= b[i] - a[i];
        double sum = 0, sum2 = 0;
        var x = new vector(dim);
        var rnd = new System.Random();
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < dim; k++) x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            double fx = f(x); sum += fx; sum2 += fx * fx;
        }
        double mean = sum / N, sigma = Sqrt(sum2 / N - mean * mean);
        var result = (mean * V, sigma * V / Sqrt(N));
        return result;
    }
    static double corput(int n, int b = 2) {
        double q = 0;
        double bk = 1.0 / b;
        while (n > 0) {
            q += (n % b) * bk;
            n /= b;
            bk /= b;
        }
        return q;
    }
    static vector halton(int n, int d) {
        int[] bases = new int[18] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
        int maxd = bases.Length;
        vector x = new vector(d);
        for (int i = 0; i < d; i++) {
            x[i] = corput(n, bases[i]);
        }
        return x;
    }

    static (double, double) haltonmc(System.Func<vector, double> f, vector a, vector b, int N) {
        int dim = a.size; double V = 1; for (int i = 0; i < dim; i++) V *= b[i] - a[i];
        double sum = 0, sum2 = 0;
        var x = new vector(dim);
        var x2 = new vector(dim);
        var rnd = new System.Random();
        for (int i = 0; i < N; i++) {
            vector y = halton(rnd.Next(), dim);
            vector y2 = halton(rnd.Next(), dim);
            for (int k = 0; k < dim; k++) x[k] = a[k] + y[k] * (b[k] - a[k]);
            for (int k = 0; k < dim; k++) x2[k] = a[k] + y2[k] * (b[k] - a[k]);
            double fx = f(x); sum += fx;
            double fx2 = f(x2); sum2 += fx2;
        }
        double mean = sum / N;
        var result = (mean * V, Abs(mean * V - V * sum2 / N));
        return result;
    }

    static (double, double) stratified_sampling(System.Func<vector, double> f, vector a, vector b, int N)
    {
        var rnd = new System.Random();
        int dim = a.size; double V = 1; for (int i = 0; i < dim; i++) V *= b[i] - a[i];
        int nmin = 16 * dim;
        if (N < nmin) return plainmc(f, a, b, N);
        vector n_right = new vector(dim);
        vector n_left = new vector(dim);
        vector x = new vector(dim);
        vector mean_left = new vector(dim);
        vector mean_right = new vector(dim);
        vector sum2_left = new vector(dim);
        vector sum2_right = new vector(dim);
        vector var_left = new vector(dim);
        vector var_right = new vector(dim);
        double mean = 0;
        double sum2 = 0;
        for (int k = 0; k < dim; k++)
        {
            mean_left[k] = 0;
            mean_right[k] = 0;
            n_left[k] = 0;
            n_right[k] = 0;
            sum2_right[k] = 0;
            sum2_left[k] = 0;
        }
        for (int i = 0; i < nmin; i++) {
            for (int k = 0; k < dim; k++) x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            double fx = f(x);
            mean += fx;
            for (int k = 0; k < dim; k++)
            {
                if (x[k] > (a[k] + b[k]) / 2.0)
                {
                    n_right[k]++;
                    mean_right[k] += fx;
                    sum2_right[k] += fx * fx;
                }
                else
                {
                    n_left[k]++;
                    mean_left[k] += fx;
                    sum2_left[k] += fx * fx;
                }
            }
        }
        mean /= nmin;
        for (int k = 0; k < dim; k++)
        {
            mean_left[k] /= n_left[k];
            mean_right[k] /= n_right[k];
            var_left[k] = sum2_left[k] / n_left[k] - mean_left[k] * mean_left[k];
            var_right[k] = sum2_right[k] / n_right[k] - mean_right[k] * mean_right[k];
        }
        int kdiv = 0;
        double maxvar = 0;
        for (int k = 0; k < dim; k++) {
            double var = Abs(mean_right[k] - mean_left[k]);
            if (var > maxvar)
            {
                maxvar = var;
                kdiv = k;
            }
        }
        vector a2 = a.copy();
        vector b2 = b.copy();
        a2[kdiv] = (a[kdiv] + b[kdiv]) / 2.0;
        b2[kdiv] = (a[kdiv] + b[kdiv]) / 2.0;
        int N_remain = N - nmin;
        int N_left = (int)Round(N_remain * (var_left[kdiv] / (var_left[kdiv] + var_right[kdiv])));
        int N_right = N_remain - N_left;
        if (N_left == 0) N_left = 1;
        if (N_right == 0) N_right = 1;
        var integ_left = stratified_sampling(f, a, b2, N_left);
        var integ_right = stratified_sampling(f, a2, b, N_right);
        return (integ_left.Item1 + integ_right.Item1, Sqrt(integ_left.Item2 * integ_left.Item2 + integ_right.Item2 * integ_right.Item2));


    }
    public static int Main(string[] args)
    {
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Plain Monte Carlo integration");
        //int[] Ns = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        System.Func<vector, double> f_A1 = (vec) => Exp(-(Pow(vec[0], 2) + Pow(vec[1], 2))); //Gaussian function, skal give 0.557746. Grænser 0,1 og 0,1
        System.Func<vector, double> f_A2 = (vec) => 1 / (1 - Cos(vec[0]) * Cos(vec[1]) * Cos(vec[2])) * Pow(PI, -3); // 
        //System.Func<vector, double> f_A3 = (vec) => Log(1+vec[0]*vec[1]);// Med singularitet, grænser 0,1 og 0,1. skal give 0.29
        WriteLine("     - Test implemented plain MC on ∫_0^1∫_0^1 exp(-(x^2+y^2)) dx dy=0.557746. The estimated and actual error as a funciton of N is plotted in MCA.gnuplot.svg");
        WriteLine("     - As can be seen in the plot, the estimated error follows the actual error pretty well. As the estimated error as calculated from a*1/√N, with a being a constant, it is confirmed that the actual error scales with 1/√N");

        var data1 = new System.IO.StreamWriter("data_A.txt", append: true);
        for (int i = 1; i < 100000; i += 250)
        {
            data1.WriteLine($"{i}   {plainmc(f_A1, new vector(0, 0), new vector(1, 1), i).Item2}   {Abs(plainmc(f_A1, new vector(0, 0), new vector(1, 1), i).Item1 - 0.557746)}");
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
        WriteLine($"       Computed value with N=500000: {plainmc(f_A2, new vector(0, 0, 0), new vector(PI, PI, PI), 500000).Item1}");
        WriteLine("        Table value = 1.3932039296856768591842462603255 ");
        WriteLine();
        WriteLine("B. Quasi-random sequences");
        WriteLine(" - Test implemented Halton MC on ∫_0^1∫_0^1 exp(-(x^2+y^2)) dx dy=0.557746.");
        var data2 = new System.IO.StreamWriter("data_B.txt", append: true);
        for (int i = 1; i < 100000; i += 250)
        {
            data2.WriteLine($"{i}   {haltonmc(f_A1, new vector(0, 0), new vector(1, 1), i).Item2}   {Abs(haltonmc(f_A1, new vector(0, 0), new vector(1, 1), i).Item1 - 0.557746)}");
        }
        data2.Close();
        WriteLine("     - The estimated and actual error as a funciton of N is plotted in MCB.gnuplot.svg ");
        WriteLine("     - As can be seen, the estimated error varies to a higher degree compared to the actual error.");
        WriteLine("     - The estimated error of plain MC is also plotted agagin in MCB.gnuplot.svg, and is comapred to the actual and estimated error of the halton MC");
        WriteLine("     - As mentioned earlier, the estimated error of the plain MC scales with 1/√N and accurately describes the scaling of the actual plain MC error.");
        WriteLine("     - When comparing this the the actual error of Halton MC, it is seen how the actual Halton MC error decreases faster than the 1/√N plain MC error");
        WriteLine("     - The Halton MC therefor improves the scaling with N compared to the plain MC.");
        WriteLine();
        WriteLine("C. Stratified sampling");
        WriteLine("     - The implemented stratified sampling routine is tested on ∫_0^1∫_0^1 exp(-(x^2+y^2)) dx dy=0.557746 ");
        var resultC = stratified_sampling(f_A1, new vector(0, 0), new vector(1, 1), 5000);
        WriteLine($"     - Computed value with N=5000: {resultC.Item1}");
        WriteLine($"     - Computed estimated error: {resultC.Item2}");
        WriteLine($"     - Actual error {Abs(resultC.Item1 - 0.557746)}");
        WriteLine("     - The actual and estimated error of the stratified sampling is plotted as a function of N in MCC.gnuplot.svg");
        WriteLine("     - When compared to the estimated error for the plain MC, the error scaling with N is significantly improved with the stratified sampling.");
        WriteLine("     - Its also improved with respect to ther Halton MC");
        var data3 = new System.IO.StreamWriter("data_C.txt", append: true);
        for (int i = 1; i < 100000; i += 250)
        {
            data3.WriteLine($"{i}   {stratified_sampling(f_A1, new vector(0, 0), new vector(1, 1), i).Item2}   {Abs(stratified_sampling(f_A1, new vector(0, 0), new vector(1, 1), i).Item1 - 0.557746)}");
        }
        data3.Close();

        return 0;
    }
}

//data_A.txt" index 1 using 1:3 with lines lw 0.5 title "Actual error, Halton MC"
//,"data_A.txt" index 1 using 1:2 with lines lw 0.5 title "Estimated error, Halton MC"\
