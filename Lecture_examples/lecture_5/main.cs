class main{
static int Main(string[] args ){
	foreach(string arg in args) System.Console.Out.WriteLine($"stdout: {arg}");
	foreach(string arg in args) System.Console.Error.WriteLine($"stderr: {arg}");
	for(
		string line = System.Console.In.ReadLine(); 
		line!=null; 
		line = System.Console.In.ReadLine()){
			System.Console.Out.WriteLine($"line for stdn: {line}");
		}
		return 0;
	}
}
