using static System.Math;
using static System.Console;
static class main{
        public static bool approx(double a, double b, double acc=1e-9, double eps=1e-9){
	        if(Abs(b-a) <= acc) return true;
	        if(Abs(b-a) <= Max(Abs(a),Abs(b))*eps) return true;
	        return false;
        }
    static int Main(string[] args){
        WriteLine("1. Maximum/minimum representable integers");
        int i = 1;
        while(i+1>i) {
            i++;
        } 
        WriteLine($"my max int = {i} ");
        WriteLine($"int.MaxValue = {int.MaxValue}");
        WriteLine($"my max int = int.MaxValue? {i == int.MaxValue}"); 

        int t = 1;
        while(t-1<t){
            t--;
        }
        WriteLine($"my min int = {t} ");
        WriteLine($"int.MinValue = {int.MinValue}");
        WriteLine($"my min int = int.MinValue? {t == int.MinValue} \n");

        WriteLine("2. The machine epsilon");
        double x = 1;
        while(1+x != 1){
            x/=2;
        }
        x*=2;
        WriteLine($"my machine epsilon double = {x}");
        WriteLine($"2^-52={System.Math.Pow(2,-52)}");
        WriteLine($"Is my machine epsilon double equal to 2^-52? {x==System.Math.Pow(2,-52)}");
        float y=1F;
        while((float)(1F+y) != 1F){
            y/=2F;
        }
        y*=2F;
        WriteLine($"my machine epsilon float = {y}");
        WriteLine($"2^-23={System.Math.Pow(2,-23)}");
        WriteLine($"Is my machine epsilon flaot equal to 2^-23? True to 6 decimals \n ");  

        WriteLine("3. suppose tiny=epsilon/2");
        double epsilon=Pow(2,-52);
        double tiny=epsilon/2;
        double a = 1+tiny+tiny;
        double b = tiny+tiny+1;
        WriteLine($"tiny={tiny}");
        WriteLine($"a=1+tiny+tiny={a}");
        WriteLine($"b=tiny+tiny+1={b}");
        WriteLine($"a==b ? {a==b}");
        WriteLine($"a>1  ? {a>1}");
        WriteLine($"b>1  ? {b>1}");
        WriteLine("the reason for the above, is the rounding of doubles when the values are saved \n");

        WriteLine("4. Comparing doubles"); 
        double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
        double d2 = 8*0.1; 
        WriteLine($"d1={d1:e15}");
        WriteLine($"d2={d2:e15}");
        WriteLine($"d1==d2 ? => {d1==d2}"); 
        
        WriteLine($"d1==d2 with approx function?  => {approx(d1,d2)}");     
        return 0;
    }
}
