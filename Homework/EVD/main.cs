using static System.Console;
using static System.Math;
public static class jacobi{
    public static void timesJ(matrix A, int p, int q, double theta){
	    double c=Cos(theta),s=Sin(theta);
	    for(int i=0;i<A.size1;i++){
		double aip=A[i,p],aiq=A[i,q];
		A[i,p]=c*aip-s*aiq;
		A[i,q]=s*aip+c*aiq;
		}
    }
    public static void Jtimes(matrix A, int p, int q, double theta){
	    double c=Cos(theta),s=Sin(theta);
	    for(int j=0;j<A.size1;j++){
		    double apj=A[p,j],aqj=A[q,j];
		    A[p,j]= c*apj+s*aqj;
		    A[q,j]=-s*apj+c*aqj;
		}

    }
    public static (vector,matrix) cyclic(matrix M){
	    matrix A=M.copy();
	    matrix V=matrix.id(M.size1);
	    vector w=new vector(M.size1);
        int n = A.size1;
	    bool changed;
        do{
	        changed=false;
	        for(int p=0;p<n-1;p++)
	        for(int q=p+1;q<n;q++){
		        double apq=A[p,q], app=A[p,p], aqq=A[q,q];
		        double theta=0.5*Atan2(2*apq,aqq-app);
		        double c=Cos(theta),s=Sin(theta);
		        double new_app=c*c*app-2*s*c*apq+s*s*aqq;
		        double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
		        if(new_app!=app || new_aqq!=aqq)/*do rotation */{
			        changed=true;
			        timesJ(A,p,q, theta); // A←A*J 
			        Jtimes(A,p,q,-theta); // A←JT*A 
			        timesJ(V,p,q, theta); // V←V*J
			    }
	        }
        }while(changed);
        for(int i=0; i<w.size; i++){
            w[i] = A[i, i];
        }
	    return (w,V);
	}
}
static class main{
    static int Main(string[] args){
        System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
        int rmax = 0; //default values
        double dr = 0; //default values
        string task = "";
        int N = 0;
        foreach(var arg in args) { 
            var words = arg.Split(':');
            if(words[0]=="-rmax") rmax=int.Parse(words[1]);
            if(words[0]=="-dr"  ) dr  =double.Parse(words[1]);
            if(words[0]=="-task") task=words[1];
            if(words[0]=="-size") N=int.Parse(words[1]);
        }
        if(task ==""){    
            System.Random rnd = new System.Random();
            int n = rnd.Next(3,8);
            matrix A = new matrix(n, n);
            for(int i=0; i<n; i++){
                for(int j=i; j<n; j++){
                    int tal = rnd.Next(1,100);
                    if(i==j){
                        A.set(i, j, tal);
                    }
                    else{
                        A.set(i, j, tal);
                        A.set(j, i, tal);
                    }
                }
            }
            matrix V = jacobi.cyclic(A).Item2;
            vector w = jacobi.cyclic(A).Item1;
            matrix D = new matrix(w);
            WriteLine("A. Jacobi diagonalization with cyclic sweeps");
            WriteLine();
            WriteLine("Check implementation:");
            WriteLine($"V^TV=I? {matrix.id(n).approx(V.T*V)}");
            WriteLine($"VV^T=I? {matrix.id(n).approx(V*V.T)}");
            WriteLine($"VDV^T=A? {A.approx(V*D*V.T)}");
            WriteLine($"V^TAV=D? {D.approx(V.T*A*V)}");
            WriteLine();
            WriteLine("B. Hydrogen Atom, s-wave radial Schrödinger equation on a grid");
            WriteLine();
            WriteLine("Firstly, to prove that f(r → 0) = r-r² → 0, the expression is rewritten.");
            WriteLine("r-r²=r*(1-r)");
            WriteLine("When r becomes smaller in the above expression, the whole expression goes towards r. As r → 0, r*(1-r) → 0.");
            WriteLine("The plots showing convergence with rmax and ∆r can be seen in dr.gnuplot.svg and rmax.gnuplot.svg.");
            WriteLine("Plot of ε[0] as a function of ∆r with rmax=10 can be seen in dr.gnuplot.svg ");
            WriteLine("From this plot ε[0] only comes close to converging at very low ∆r with rmax at this value. ");
            WriteLine("Plot of ε[0] as a function of r_max with ∆r=0.3 can be seen in rmax.gnuplot.svg");
            WriteLine("From this plot ε[0] seems to converge at rmax above 8, however it does not converge to the right ε[0], which should be -0.5 Hartree.");
            WriteLine("The lowest 4 eigenfunctions plotted against the analytical solutions can be seen in fr.gnuplot.svg");
            WriteLine("For calculations rmax was chosen to be 50, as some of the higher quantum number analytical functions only goes towards 0 at this rmax.");
            WriteLine("Furthermore, dr was chosen to be 0.3.");
            WriteLine();
            WriteLine("C. Scaling and optimization");
            WriteLine();
            WriteLine("Here the scaling of the diagonalization time with matrix size is investigated.");
            WriteLine("This is done by running diagonalization of various N size matrices in parallel.");
            WriteLine("The results can be seen in times.gnuplot.svg. By fitting both a*N^3 and b*N^2 to the data, it seems as the time scales with N^2.");
            
            return 0;
            }
        if(task=="fr"){
            int npoints = (int)(rmax/dr)-1;
            vector r = new vector(npoints);
            for(int i=0;i<npoints;i++)r[i]=dr*(i+1);
            matrix H = new matrix(npoints,npoints);
            for(int i=0;i<npoints-1;i++){
                H[i,i]  =-2*(-0.5/dr/dr);
                H[i,i+1]= 1*(-0.5/dr/dr);
                H[i+1,i]= 1*(-0.5/dr/dr);
            }
            H[npoints-1,npoints-1]=-2*(-0.5/dr/dr);
            for(int i=0;i<npoints;i++)H[i,i]+=-1/r[i];
            //H er nu bygget i skal herefter diagonaliseres
            matrix V_H = jacobi.cyclic(H).Item2;
            double Const = 1/Sqrt(dr);
            double sum = 0;
            for(int i=0; i<V_H.size2;i++){
               sum += Pow(V_H[3,i],2); 
            }
            for(int k=0; k<4;k++){
                for(int i=0; i<r.size;i++){
                    WriteLine($"{r[i]}   {Const*V_H[i,k]}");
                }
                WriteLine();
                WriteLine();
            }
            return 0;

        }
        if(task=="dr" || task=="rmax"){
            int npoints = (int)(rmax/dr)-1;
            vector r = new vector(npoints);
            for(int i=0;i<npoints;i++)r[i]=dr*(i+1);
            matrix H = new matrix(npoints,npoints);
            for(int i=0;i<npoints-1;i++){
                H[i,i]  =-2*(-0.5/dr/dr);
                H[i,i+1]= 1*(-0.5/dr/dr);
                H[i+1,i]= 1*(-0.5/dr/dr);
            }
            H[npoints-1,npoints-1]=-2*(-0.5/dr/dr);
            for(int i=0;i<npoints;i++)H[i,i]+=-1/r[i];
            //H er nu bygget i skal herefter diagonaliseres
            vector e = jacobi.cyclic(H).Item1;
            WriteLine($"{e[0]}  {dr}    {rmax}");
            return 0;
        }
        if(task=="time"){
            matrix H = new matrix(N);
            var rnd = new System.Random();
                for (int i = 0; i < H.size1; i++) {
                    for(int t = i; t<H.size1; t++){
                        double q = 100.0*(rnd.NextDouble()-0.5);
                        H[i,t] = q;
                        H[t,i] = q;
                    }
                }
            vector eigen = jacobi.cyclic(H).Item1;
            return 0;
        }
        return 0;
    }
}    
