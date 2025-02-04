using static System.Math;
static class main{
    static int Main(string[] args){
        double sqrt2 = Sqrt(2.0);
        System.Console.WriteLine($"√2 = {sqrt2}, test: √2^2 = {sqrt2*sqrt2} (should equal 2)");
        double var_2 = Pow(2.0, 1.0/5.0); 
        System.Console.WriteLine($"2^(1/5) = {var_2}, test: (2^(1/5))^5 = {Pow(var_2,5.0)} (should equal 2)");
        double var_3 = Pow(E, PI);
        System.Console.WriteLine($"e^π = {var_3}, test: (e^π)^(1/π) = {Pow(var_3,1/PI)} (should equal {E})");
        double var_4 = Pow(PI, E);
        System.Console.WriteLine($"π^e = {var_4}, test: (π^e)^(1/e) = {Pow(var_4, 1/E)} (should equal {PI})");
                // Creating a list of integers
        System.Console.WriteLine("Results using the fgamma function");
        double[] numbers = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        foreach (double num in numbers){
            double Gamma = sfuns.fgamma(num);
            System.Console.WriteLine($"Gamma({num}) = {Gamma}");
        }
        System.Console.WriteLine("Correct answers: Gamma(1)=1, Gamma(2)=1, Gamma(3)=2, Gamma(4)=6, Gamma(5)=24, Gamma(6)=120, Gamma(7)=720, Gamma(8)=5040, Gamma(9)=40320, Gamma(10)=362880");
        System.Console.WriteLine("Results using the lngamma function");
        foreach (double i in numbers){
            double Gamma = sfuns.lngamma(i);
            System.Console.WriteLine($"lngamma({i}) = {Gamma}");
        }
        System.Console.WriteLine("Correct answers (rounded): lngamma(1)=0, lngamma(2)=0, lngamma(3)=0.693147, lngamma(4)=1.791759, lngamma(5)=3.178054, lngamma(6)=4.787491, lngamma(7)=6.579251, lngamma(8)=8.525161, lngamma(9)=10.604603 lngamma(10)=12.801827");
    return 0;
    }
}
