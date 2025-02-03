//using static System.Console;
using static System.Math;
static class main{
    static string class_variable="class variable";
    static double square(double x){
        double tmp=1;
        double result = x*x*tmp;
        return result;
        }
    static int Main(string[] args){
        System.Console.WriteLine("hello");
        System.Console.WriteLine(class_variable);
        double x=System.Math.Sin(1.0);
        double y=Cos(1.0);
        System.Console.WriteLine($"x={x}, y={y}");
    return 0;
    }
}
