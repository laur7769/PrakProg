using static System.Console;
public static class QR{
    public static bool upper_triangular(matrix A){
        int rows = A.size1;
        int cols = A.size2;
        if (rows != cols){
            return false;
        }
        for (int i = 1; i < rows; i++){ 
            for (int j = 0; j < i; j++){ 
                if (A[i, j] != 0)
                    return false; 
            }
        }
        return true;
    }
   public static (matrix,matrix) decomp(matrix A){
      matrix Q=A.copy();
      matrix R=new matrix(A.size2,A.size2);
      for(int i=0; i<A.size2; i++){
        R[i,i] = Q[i].norm();
        Q[i]/=R[i,i];
        for(int j=i+1; j<A.size2;j++){
            R[i,j]=Q[i].dot(Q[j]);
            Q[j]-=Q[i]*R[i,j];
        }
      }
      return (Q,R);
      }
    
    static vector backsub(matrix U, vector c){
            for(int i=c.size-1; i>=0; i--){
                double sum = 0;
                for(int k=i+1; k<c.size; k++){
                    sum +=U[i,k]*c[k];
                }
            c[i]=(c[i] - sum)/U[ i , i ];    
            }
            
            return c;
    }    
    static vector forwsub(matrix U, vector c){
            for(int i=0; i<c.size; i++){
                double sum = 0;
                for(int k=1; k<i-1; k++){
                    sum +=U[i,k]*c[k];
                    c [ i ]=(c[ i]-sum)/U[ i , i ];
                }
            }
            return c;
        }
    public static vector solve(matrix Q, matrix R, vector b){
        if(upper_triangular(R)){
            return backsub(R, Q.T*b);
        }
        else{
            return forwsub(R, Q.T*b);
        }
        
   }
   public static double det(matrix R){
    double determinant = 1;
    for(int i=0; i<R.size1; i++){
        determinant *= R[i,i];
    }
    return determinant;
   }
   public static matrix inverse(matrix Q,matrix R){
    if(Q.size1 == Q.size2){
        matrix A_inv = new matrix(Q.size1, Q.size2);
        for(int n=0; n<Q.size1;n++){
            vector e = new vector(Q.size1);
            for(int i=0; i<e.size; i++){
                if(i==n){
                    e[i]=1;
                }
                else{
                    e[i]=0;
                }
            }
            vector x = solve(Q, R, e);
            for(int i=0; i<A_inv.size1; i++){
                    A_inv.set(i, n, x[i]);
            }

        }
        return A_inv;
    }
    else{
        Error.WriteLine("No inverse as matrix is not square");
        return new matrix(0);
    }
   }
}
class main{

    static int Main(string[] args){
        System.Random rand = new System.Random();
        int n = rand.Next(0,10);
        int m = rand.Next(0,10);
        while(m>n){
            m = rand.Next(0,10);
        }
        matrix A = new matrix(n, m);
        for(int i=0; i<n; i++){
            for(int k=0; k<m; k++){
                A.set(i, k, rand.Next(0, 20));
            }
        }
        WriteLine("A. Solving linear equations using QR-decomposition by modified Gram-Schmidt orthogonalization");
        WriteLine();
        matrix Q = QR.decomp(A).Item1;
        matrix R = QR.decomp(A).Item2;
        WriteLine("Check decomp method:");
        WriteLine($"Is R upper triangular? {QR.upper_triangular(R)}");
        matrix ID = new matrix(Q.size2, Q.size2); 
        for(int i=0;i<ID.size1;i++){
		ID.set(i,i,1);
		    for(int j=i+1;j<ID.size2;j++){
			ID.set(i,j,0);
            ID.set(j,i,0);
		    }
	    }
        WriteLine($"Is Q^TQ=I? {ID.approx(Q.T*Q)}");
        WriteLine($"Is QR=A? {A.approx(Q*R)}");
        WriteLine("Check solve method:");
        matrix A2 = new matrix(3,3);
        for(int i=0; i<3; i++){
            for(int k=0; k<3; k++){
                A2.set(i, k, rand.Next(0, 20));
            }
        }
        var rnd = new System.Random(1); /* or any other seed */
        double x = rnd.NextDouble();
        double y = rnd.NextDouble();
        double z = rnd.NextDouble();
        vector b = new vector(x, y, z);
        matrix Q2 = QR.decomp(A2).Item1;
        matrix R2 = QR.decomp(A2).Item2;
        vector x2 = QR.solve(Q2, R2, b);
        WriteLine($"Is Ax=b?: {b.approx(A2*x2)}");
        WriteLine();
        WriteLine("B. Matrix inverse by Gram-Schmidt QR factorization");
        WriteLine();
        matrix B = QR.inverse(Q2, R2);
        matrix ID2 = new matrix(Q2.size2, Q2.size2); 
        for(int i=0;i<ID2.size1;i++){
		    ID2.set(i,i,1);
		    for(int j=i+1;j<ID2.size2;j++){
			    ID2.set(i,j,0);
                ID2.set(j,i,0);
		    }
	    }
        WriteLine("Check inverse method:");
        WriteLine($"Is AB=I? {ID2.approx(A2*B)}");


        return 0;
    }
}