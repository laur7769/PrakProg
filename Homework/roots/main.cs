using static System.Console;
using static System.Math;
using System.Collections.Generic;


class main{

    public static matrix jacobian
    (System.Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
	    if(dx == null) dx = x.map(xi => Abs(xi)*Pow(2,-26));
	    if(fx == null) fx = f(x);
	    matrix J=new matrix(x.size);
	    for(int j=0;j < x.size;j++){
		    x[j]+=dx[j];
		    vector df=f(x)-fx;
		    for(int i=0;i < x.size;i++) J[i,j]=df[i]/dx[j];
		    x[j]-=dx[j];
		}
	    return J;
    }
    public static vector newton(
	System.Func<vector,vector>f /* the function to find the root of */
	,vector start        /* the start point */
	,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	,vector δx=null      /* optional δx-vector for calculation of jacobian */
	){
        vector x=start.copy();
        vector fx=f(x); 
        vector z;
        vector fz;
        do{ /* Newton's iterations */
	        if(fx.norm() < acc) break; /* job done */
	        matrix J=jacobian(f,x,fx,δx);
	        var QRJ = QR.decomp(J);
	        vector Dx = QR.solve(QRJ.Item1, QRJ.Item2, -fx); /* Newton's step */
	        double λ=1;
            double λmin = 1.0/128.0;
	        do{ /* linesearch */
		        z=x+λ*Dx;
		        fz=f(z);
		        if( fz.norm() < (1-λ/2)*fx.norm() ) break;
		        if( λ < λmin ) break;
		        λ/=2;
		    }while(true);
	        x=z; fx=fz;
	    }while(true);
        return x;
    }

    public static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        WriteLine("A. Newton's method with numerical Jacobian and back-tracking line-search");
        WriteLine("     - Test root finding routine on f(x)=(x-2)(x+3):");
        System.Func<vector,vector>f1 = (x) => new vector((x[0]-2)*(x[0]+3));
        WriteLine($"       Computed root with start 1.0 and default accuracy: {newton(f1, new vector(1.0) )[0]}");
        WriteLine($"       Computed root with start -1.0 and default accuracy: {newton(f1, new vector(-1.0) )[0]}");
        WriteLine($"       Actual roots: 2 and -3");
        WriteLine("     - Test root finding routine on f(x,y)=[x^2+y^2-4, x-y]:");
        System.Func<vector,vector>f2 = (vec) => new vector(vec[0]*vec[0]+vec[1]*vec[1]-4, vec[0]-vec[1]);
        vector res1 =newton(f2, new vector(1.0, 1.0) );
        vector res2 =newton(f2, new vector(-1.0, -1.0) );
        WriteLine($"       Computed root with start (1.0, 1.0) and default accuracy: {(res1[0],res1[1])}");
        WriteLine($"       Computed root with start (-1.0, -1.0) and default accuracy: {(res2[0], res2[1])}");
        WriteLine($"       Actual roots: (√2, √2) and (-√2, -√2)");
        WriteLine("     - Find extremum(s) of the Rosenbrock's valley function: f(x,y) = (1-x)^2+100(y-x^2)^2");
        WriteLine("        Analytical gradient: ∂f/∂x=2x-2+400x^3-400xy,   ∂f/∂y=200y-200x^2");
        WriteLine("        It is known that the Rosenbrock's valley function only has one extremum, a global minimum at (1,1)");
        System.Func<vector,vector>Rosen_grad = (vec) => new vector(2*vec[0]-2+400*Pow(vec[0],3)-400*vec[0]*vec[1], 200*vec[1]-200*vec[0]*vec[0]);
        vector ros_res = newton(Rosen_grad, new vector(1.0, 1.0), acc:1e-4);
        WriteLine($"        Computed root with start (7.0, -2.0) and acc=1e-4: {(ros_res[0],ros_res[1])}");
        WriteLine("     - Find minimum(s) ofthe Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
        WriteLine("        Analytical gradient: ∂f/∂x=2x-2+400x^3-400xy,   ∂f/∂y=200y-200x^2");
        System.Func<vector,vector>Himmel_grad = (vec) => new vector(4*Pow(vec[0],3)+4*vec[0]*vec[1]-42*vec[0]+2*vec[1]*vec[1]-14, 2*vec[0]*vec[0]-22+4*vec[0]*vec[1]+4*Pow(vec[1],3)-26*vec[1]);
        WriteLine("        Analytical gradient: 4x^2+4xy-42x+2y^2-14,   ∂f/∂y=2x^2-22+4yx+4y^3-26y");
        WriteLine("        It is known that Himmelblau's funciton has 4 minimums: (-2.805118, 3.131312), (3.0, 2.0), (-3.779310, -3.283186) and (3.584428, -1.848126)") ;
    var points = new List<(double, double)> {(-10, -10), (-10, 10), (10, -10), (10, 10)};
    vector xx = new vector(2);
        foreach ((double p1, double p2) point in points) {
            xx[0] = point.p1; xx[1] = point.p2;
            vector solution = newton(Himmel_grad, xx, acc:1e-10);
            WriteLine($"        Computed root with start {(point.p1, point.p2)} and acc:1e-10: {(solution[0], solution[1])}");
        }
    return 0;
    }

}
