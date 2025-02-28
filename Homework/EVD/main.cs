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
        foreach(var arg in args) { 
            var words = arg.Split(':');
            if(words[0]=="-rmax") rmax=int.Parse(words[1]);
            if(words[0]=="-dr"  ) dr  =double.Parse(words[1]);
        }
        if(rmax==0 && dr==0){    
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
            return 0;
            }
        if(rmax==10 && dr==0.025){
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
               sum += Pow(V_H[0,i],2); 
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
        else{
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
    }
}    
