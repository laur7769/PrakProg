using static System.Console;
using static System.Math;
static class main{
    static int Main(){
        WriteLine("1.  ");
        complex minus_one = new complex(-1,0);
        var sqrt_minus_one = cmath.sqrt(minus_one);
        if(sqrt_minus_one.approx(new complex(0,1)) || sqrt_minus_one.approx(new complex(0,-1))){
            WriteLine("√-1 calculated correctly");
        }
        else{
            WriteLine("√-1 calculated wrong");
            WriteLine($"Results = {sqrt_minus_one}");
        }
        complex i = new complex(0,1);
        var sqrt_i = cmath.sqrt(i);
        if(sqrt_i.approx(new complex(1.0/cmath.sqrt(2.0),1.0/cmath.sqrt(2.0)))){
            WriteLine("√i calculated correctly");
        }
        else{
            WriteLine("√i calculated wrong");
            WriteLine($"Results = {sqrt_i}");
        }
        var e_i = cmath.exp(new complex(0,1));
        if(e_i.approx(new complex(cmath.cos(1),cmath.sin(1)))){
            WriteLine("e^i calculated correctly");
        }
        else{
            WriteLine("e^i calculated wrong");
            WriteLine($"Results = {e_i}");
        }
        var e_ipi = cmath.exp(new complex(0,1)*PI);
        if(e_ipi.approx(new complex(-1,0))){
            WriteLine("e^iπ calculated correctly");
        }
        else{
            WriteLine("e^iπ calculated wrong");
            WriteLine($"Results = {e_ipi}");
        }
        var i_i = cmath.pow(new complex(0,1),new complex(0,1));
        if(i_i.approx(cmath.exp(-PI/2.0))){
            WriteLine("i^i calculated correctly");
        }
        else{
            WriteLine("i^i calculated wrong");
            WriteLine($"Results = {i_i}");
        }
        var ln_i = cmath.log(new complex(0,1));
        if(ln_i.approx(new complex(0,PI/2.0))){
            WriteLine("ln(i) calculated correctly");
        }
        else{
            WriteLine("ln(i) calculated wrong");
            WriteLine($"Results = {ln_i}");
        }
        var sin_ipi = cmath.sin(new complex(0,PI));
        if(sin_ipi.approx((new complex(0,cmath.exp(PI)/2.0))-(new complex(0,cmath.exp(-PI)/2.0)))){
            WriteLine("sin(iπ) calculated correctly");
        }
        else{
            WriteLine("sin(iπ) calculated wrong");
            WriteLine($"Results = {sin_ipi}");
        }
        return 0;
    }
}
