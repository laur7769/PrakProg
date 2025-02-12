using static System.Console;
using static System.Math;
public class genlist<T>{
	public T[] data;
	public int size => data.Length; // property
	public T this[int i] => data[i]; // indexer
	public genlist(){ data = new T[0]; }
	public void add(T item){ /* add item to the list */
		T[] newdata = new T[size+1];
		System.Array.Copy(data,newdata,size);
		newdata[size]=item;
		data=newdata;
	}
    public void remove(int i){
        if(i<0 || i>size){
            Error.WriteLine("Index is out of range");
        }
        else{
            T[] newdata = new T[size-1];
            for(int t=0; t<size;t++){
                if (t!=i && t<i){
                newdata[t]=data[t]; 
                }
                if (t!=i && t>i){
                    newdata[t-1]=data[t];
                }
            }
            data = newdata;
        }
    }
}

static class main{
    static int Main(string[] args){
        var list = new genlist<double[]>();
        char[] delimiters = {' ','\t'};
        var options = System.StringSplitOptions.RemoveEmptyEntries;
        for(string line = ReadLine(); line!=null; line = ReadLine()){
	        var words = line.Split(delimiters,options);
	        int n = words.Length;
	        var numbers = new double[n];
	        for(int i=0;i<n;i++) numbers[i] = double.Parse(words[i]);
	        list.add(numbers);
       	    }
        for(int i=0;i<list.size;i++){
	        var numbers = list[i];
	        foreach(var number in numbers)Write($"{number : 0.00e+00;-0.00e+00} ");
	        WriteLine();
        }
        WriteLine("The same table but with the first line removed");

        list.remove(0);
        for(int i=0;i<list.size;i++){
	        var numbers = list[i];
	        foreach(var number in numbers)Write($"{number : 0.00e+00;-0.00e+00} ");
	        WriteLine();
        }
        return 0;
    }
}
