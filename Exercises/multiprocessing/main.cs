using static System.Console;
public class data{ public int a,b; public double sum;}
class main{
    public static void harm(object obj){
	var arg = (data)obj;
	arg.sum=0;
	for(int i=arg.a;i<arg.b;i++)arg.sum+=1.0/i;
	}
    public static int Main(string[] args){
        int nthreads = 1, nterms = (int)1e8; /* default values */
    foreach(var arg in args) { 
        var words = arg.Split(':');
        if(words[0]=="-threads") nthreads=int.Parse(words[1]);
        if(words[0]=="-terms"  ) nterms  =(int)float.Parse(words[1]);
    }
    data[] parameters = new data[nthreads];
    for(int i=0;i<nthreads;i++) {
        parameters[i] = new data();
        parameters[i].a = 1 + nterms/nthreads*i;
        parameters[i].b = 1 + nterms/nthreads*(i+1);
    }
    parameters[parameters.Length-1].b=nterms+1; /* the enpoint might need adjustment */
    var threads = new System.Threading.Thread[nthreads];
    for(int i=0;i<nthreads;i++) {
        threads[i] = new System.Threading.Thread(harm); /* create a thread */
        threads[i].Start(parameters[i]); /* run it with params[i] as argument to "harm" */
    }
    foreach(var thread in threads) thread.Join();
    double total=0; foreach(var p in parameters) total+=p.sum;
    WriteLine($"Total sum = {total}");

    return 0;
    }
}
