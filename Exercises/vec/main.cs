using static System.Math;
using static System.Console;
static class main{
    static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
        var rnd=new System.Random(); //En eller anden form for generator
	    var u=new vec(rnd.NextDouble(),rnd.NextDouble(),rnd.NextDouble()); //Laver en vektor med random doubles som komponenter
	    var v=new vec(rnd.NextDouble(),rnd.NextDouble(),rnd.NextDouble()); // Gør det samme som ovenstående
        u.print("u=");
	    v.print("v=");
        WriteLine($"u={u}");
        WriteLine($"v={v}");
        vec t;
        t=new vec(-u.x, -u.y, -u.z);
        (-u).print("-u=");
        t.print("t=");
        if(vec.approx(t,-u)) {
            WriteLine("test -u=[-u.x, -u.y, -u.z] passed \n");
            }
        else {WriteLine("test -u=[-u.x, -u.y, -u.z] not passed \n");
        }

        t=new vec(u.x-v.x,u.y-v.y,u.z-v.z);
	    (u-v).print("u-v =");
	    t.print    ("t   =");
	    if(vec.approx(t,u-v)) {
            WriteLine("test 'operator -' passed\n");
        }
        else {
            WriteLine("test 'operator -'  not passed\n");
        }

        t=new vec(u.x+v.x,u.y+v.y,u.z+v.z);
	    (u+v).print("u+v =");
	    t.print    ("t   =");
	    if(vec.approx(t,u+v)){
            WriteLine("test 'operator +' passed\n");
        }
        else {
            WriteLine("test 'operator +' not passed\n");
        }

        double c=rnd.NextDouble();
	    t=new vec(u.x*c,u.y*c,u.z*c);
	    var tmp=u*c; // bug in mcs
        var tmp2 = c*u;
	    tmp.print("u*c =");
        tmp2.print("c*u =");
	    t.print  ("t   =");
	    if(vec.approx(t,u*c) && vec.approx(t, c*u)) {
            WriteLine("test 'operator*' passed\n");
            }
        else {
            WriteLine("test 'operator *' not passed \n");
        }

        double dot_p = u.x*v.x+u.y*v.y+u.z*v.z;
        WriteLine($"dot_p = u.x*v.x+u.y*v.y+u.z*v.z = {dot_p}");
        WriteLine($"u.dot(v) = {u.dot(v)}");
        WriteLine($"vec.dot(u, v) = {vec.dot(u, v)}");
        if(vec.approx(dot_p,u.dot(v)) && vec.approx(dot_p, vec.dot(u, v))) {
            WriteLine("test 'dot product' passed\n");
            }
        else {
            WriteLine("test 'dot product' not passed \n");
        }

        u.print("u=");
        WriteLine("\n");
        WriteLine("u=" + u.ToString());
        WriteLine("\n If the above two expressions are the same, the ToString funciton test is passed");
        return 0;
    }
}
