using static System.Console;
using static System.Math;
class main{
    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
        string infile=null,outfile=null;
        WriteLine("1. Command-line");
        char[] split_delimiters = {' ','\t','\n'};
        var split_options = System.StringSplitOptions.RemoveEmptyEntries;
        foreach(var arg in args){
	        var words = arg.Split(':');
	        if(words[0]=="-numbers"){
		        var numbers=words[1].Split(',');
                WriteLine("x   Sin(x)               Cos(x)");
		        foreach(var number in numbers){
			        double x = double.Parse(number);
			        WriteLine($"{x}   {Sin(x)}    {Cos(x)}");
			    }
		    }
            if(words[0]=="-input")infile=words[1];
	        if(words[0]=="-output")outfile=words[1]; 
	    }
        if( infile==null || outfile==null) {
	        Error.WriteLine("wrong filename argument");
	        return 1;
	    }
        Error.Write("2. Standard input stream \n \n");
        for( string line = ReadLine(); line != null; line = ReadLine() ){
	        var numbers = line.Split(split_delimiters,split_options);
            Error.WriteLine("x   Sin(x)               Cos(x)");
	        foreach(var number in numbers){
		        double x = double.Parse(number);
		        Error.WriteLine($"{x} {Sin(x)} {Cos(x)}");
            }
        }
        var instream =new System.IO.StreamReader(infile);
        var outstream=new System.IO.StreamWriter(outfile,append:false);
        outstream.WriteLine("3. File Streams");
        outstream.WriteLine("x   Sin(x)               Cos(x)");
        for(string line=instream.ReadLine();line!=null;line=instream.ReadLine()){
	        double x=double.Parse(line);
	        outstream.WriteLine($"{x} {Sin(x)} {Cos(x)}");
        }
        instream.Close();
        outstream.Close();
        return 0;
    }
}
