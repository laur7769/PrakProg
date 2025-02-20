using static System.Console;
class main{
    public static int Main(string[] args){
        WriteLine();
        WriteLine("2. Pitfalls in multiprocessing");

        int nterms = (int)1e8; /* default values */
        foreach(var arg in args) {
            var words = arg.Split(':');
            if(words[0]=="-terms"  ) nterms  =(int)float.Parse(words[1]);
        }
        double sum=0;
        System.Threading.Tasks.Parallel.For( 1, nterms+1, (int i) => sum+=1.0/i );
        WriteLine($"Total sum = {sum}");
        WriteLine("Is slower and returns the wrong results due to 'racing condition'. The variable 'sum' is not local to a thread.");


        return 0;
    }    
}
