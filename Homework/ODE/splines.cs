using static System.Math;
using static System.Console;
public class qspline {
	public vector x,y,b,c;
	public qspline(vector xs,vector ys){
		x = xs.copy(); 
        y = ys.copy();
        int n = x.size; 
        b = new vector(n-1);
        c = new vector(n-1);
        c[0] = 0;
        vector dx = new vector(n-1);
        for(int i=0; i<dx.size; i++){
            dx[i] = x[i+1]-x[i];
        }
        vector dy = new vector(n-1);
        for(int i=0; i<dy.size; i++){
            dy[i] = y[i+1]-y[i];
        }
        vector p = new vector(n-1);
        for(int i=0; i<p.size; i++){
            p[i] = dy[i]/dx[i];
        }
        for(int i=0; i<n-2; i++ ){
            c[i+1] = (1/dx[i+1])*(p[i+1]-p[i]-c[i]*dx[i]);
        }//forward substitution
        c[n-2]/=2;
        for(int i=n-3; i>=0; i--){
            c[i]=(p[i+1]-p[i]-c[i+1]*dx[i+1])/dx[i];
        }
        for(int i=0; i<b.size; i++){
            b[i]=p[i]-c[i]*dx[i];
        }
	}
	public double evaluate(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        return y[j]+b[j]*(z-x[j])+c[j]*(Pow(z-x[j],2));
        }
    public static int binsearch(double[] x, double z){ 
	    if( z<x[0] || z>x[x.Length-1] ) throw new System.Exception("binsearch: bad z");
	    int i=0, j=x.Length-1;
	    while(j-i>1){
		    int mid=(i+j)/2;
		    if(z>x[mid]) i=mid; else j=mid;
		}
	    return i;
	}
	public double derivative(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        return b[j]+2*c[j]*(z-x[j]);
    }
	public double integral(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        double integral = 0;
        for(int i=0; i<j; i++){
            integral += y[i]*(x[i+1]-x[i])+b[i]*Pow((x[i+1]-x[i]),2)/2.0+c[i]*Pow((x[i+1]-x[i]),3)/3.0;
        }
        integral += y[j]*(z-x[j])+b[j]*Pow((z-x[j]),2)/2.0+c[j]*Pow((z-x[j]),3)/3.0;
        return integral;
    }
}

public class cspline {
    public vector x,y,b,c,d;
    public cspline(vector xs,vector ys){
        x = xs.copy(); 
        y = ys.copy();
        int n = x.size;
        b = new vector(n);
        c = new vector(n-1);
        d = new vector(n-1);
        vector h = new vector(n-1);
        //Nu bygges h= x_{i+1}-x_i
        for(int i=0; i<h.size; i++){
            h[i] = x[i+1]-x[i];
        }
        //nu bygges p_i=∆y_i/h_i
        vector p = new vector(n-1);
        for(int i=0; i<p.size;i++){
            p[i] = (y[i+1]-y[i])/h[i];
        }
        //Nu bygges en vektor med diagonalelementer D
        vector D = new vector(n);
        D[0]=2;
        for(int i=0;i<n-2;i++){
            D[i+1]=2*(h[i]/h[i+1])+2;
        }
        D[n-1]=2;
        //Nu bygges en vektor med off-diagonalelementer Q
        vector Q = new vector(n-1);
        Q[0] = 1;
        for(int i=0;i<Q.size-1; i++){
            Q[i+1] = h[i]/h[i+1];
        }
        // Nu bygges vektoren med B
        vector B = new vector(n);
        B[0]=3*p[0];
        for(int i=0; i<n-2; i++){
            B[i+1] = 3*(p[i]+p[i+1]*(h[i]/h[i+1]));
        }
        B[n-1] = 3*p[n-2];
        //Nu bygges D~
        vector D_tilt = new vector(n);
        D_tilt[0] = D[0];
        for(int i=1; i<D_tilt.size; i++){
            D_tilt[i] = D[i]-Q[i-1]/D_tilt[i-1];
        }
        //Nu bygges B~~
        vector B_tilt = new vector(n);
        B_tilt[0] = B[0];
        for(int i=1; i<B_tilt.size; i++){
            B_tilt[i]=B[i]-B_tilt[i-1]/D_tilt[i-1];
        }
        //nu bygges b
        b[n-1]=B_tilt[n-1]/D_tilt[n-1];
        for(int j=n-2; j>=0; j--){
            b[j]=(B_tilt[j]-Q[j]*b[j+1])/D_tilt[j];
        }
        //nu bygges c og 
        for(int i=0; i<c.size;i++){
            c[i]=(-2*b[i]-b[i+1]+3*p[i])/h[i];
        }
        for(int i=0; i<d.size; i++){
            d[i] = (b[i]+b[i+1]-2*p[i])/Pow(h[i],2);
        }

    }
    public static int binsearch(double[] x, double z){ 
	    if( z<x[0] || z>x[x.Length-1] ) throw new System.Exception("binsearch: bad z");
	    int i=0, j=x.Length-1;
	    while(j-i>1){
		    int mid=(i+j)/2;
		    if(z>x[mid]) i=mid; else j=mid;
		}
	    return i;
	}
    public double evaluate(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        return y[j]+b[j]*(z-x[j])+c[j]*Pow((z-x[j]),2)+d[j]*Pow((z-x[j]),3);
        }
    
    public double derivative(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        return b[j]+2*c[j]*(z-x[j])+3*d[j]*Pow((z-x[j]),2);
    }

    public double integral(double z){
        double[] xs = new double[x.size];
        for(int i=0; i<x.size; i++){
            xs[i] = x[i];
        }
        int j=binsearch(xs,z);
        double integral = 0;
        for(int i=0; i<j; i++){
            integral += y[i]*(x[i+1]-x[i])+b[i]*Pow((x[i+1]-x[i]),2)/2.0+c[i]*Pow((x[i+1]-x[i]),3)/3.0 + d[i]*Pow((x[i+1]-x[i]),4)/4.0;
        }
        integral += y[j]*(z-x[j])+b[j]*Pow((z-x[j]),2)/2.0+c[j]*Pow((z-x[j]),3)/3.0 + d[j]*Pow((z-x[j]),4)/4.0;
        return integral;
    }

}
