using static System.Console;
class main{
    public static int Main(){
        System.Globalization.CultureInfo.DefaultThreadCurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
        for(int i=0; i<=25; i++){
            double x = i * 0.1;
            WriteLine($"{x} {sfuns.erf(x)}");
        }
        WriteLine();
        WriteLine();
        for(int i=-1000; i<=1000; i++){
            if(i!=0){
                double x = i * 0.005;
                WriteLine($"{x} {sfuns.sgamma(x)}");
            }
        }
        WriteLine();
        WriteLine();
        for(int i=0; i<=1000; i++){
            if(i!=0){
                double x = i * 0.005;
                WriteLine($"{x} {sfuns.lngamma(x)}");
            }
        }
    return 0;
    }
}
