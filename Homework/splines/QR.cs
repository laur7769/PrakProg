using static System.Console;
using static System.Math;
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

