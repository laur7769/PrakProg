class main{
    static int Main(){
        double x=0;
        double y=x;
        x=9;
        int n = 3;
        double[] a = new double[n];
        for(int i=0; i<a.Length;++i)a[i]=i+1;
        double[] b=a; //når a ændres ændres b nu også, da de er reference type, og derfor nu peger til på det samme objekt
        b[0]=666;
        System.Console.WriteLine($"a[0]={a[0]}");
        foreach(double ai in a)System.Console.WriteLine(ai);
        vec v = new vec();
        //v.x=0; v.y=0; v.z=0;
        System.Console.WriteLine($"{v.x}, {v.y}, {v.z}");
        return 0;
    }
}
