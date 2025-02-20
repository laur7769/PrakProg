using static System.Console;
using System.Linq;
class main{
    public static int Main(string[] args){
        WriteLine();

        int nterms = (int)1e8; /* default values */
        foreach(var arg in args) {
            var words = arg.Split(':');
            if(words[0]=="-terms"  ) nterms  =(int)float.Parse(words[1]);
        }
        var sum = new System.Threading.ThreadLocal<double>( ()=>0, trackAllValues:true);
        System.Threading.Tasks.Parallel.For( 1, nterms+1, (int i)=>sum.Value+=1.0/i );
        double totalsum=sum.Values.Sum();
        WriteLine($"Total sum = {totalsum}");
        WriteLine("Returns the correct sum, but is still slower. Racing condition is avoided but this way of managing threads and accessing the sum variable is slower.");

        return 0;
    }    
}
