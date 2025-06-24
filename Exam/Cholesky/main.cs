using static System.Console;
using static System.Math;

public static class Cholesky
{
    public static matrix decomp(matrix A)
    {
        int dim = A.size1;
        matrix L = new matrix(dim);
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double sum = 0;
                for (int k = 0; k < j; k++)
                {
                    sum += L[i, k] * L[j, k];
                }
                if (i == j)
                {
                    L[i, j] = Sqrt(A[i, i] - sum);
                }
                else
                {
                    L[i, j] = (1.0 / L[j, j] * (A[i, j] - sum));
                }
            }
        }

        return L;
    }
    public static bool lower_triangular(matrix A)
    {
        int rows = A.size1;
        int cols = A.size2;
        if (rows != cols)
        {
            return false;
        }
        for (int i = 1; i < rows; i++)
        {
            for (int j = i + 1; j < cols; j++)
            {
                if (A[i, j] != 0)
                    return false;
            }
        }
        return true;
    }
    public static vector forwsub(matrix U, vector c)
    {
        vector y = new vector(c.size);
        for (int i = 0; i < c.size; i++)
        {
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += U[i, k] * y[k];
            }
            y[i] = (c[i] - sum) / U[i, i];
        }
        return y;
    }

    public static vector backsub(matrix U, vector c)
    {
        for (int i = c.size - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int k = i + 1; k < c.size; k++)
            {
                sum += U[i, k] * c[k];
            }
            c[i] = (c[i] - sum) / U[i, i];
        }

        return c;
    }

    public static vector solve(matrix L, vector b)
    {
        vector y = forwsub(L, b);
        vector x = backsub(L.T, y);
        return x;
    }

    public static double det(matrix A)
    {
        matrix L = decomp(A);
        double det = 1.0;
        for (int i = 0; i < L.size1; i++)
        {
            det *= L[i, i];
        }
        return det * det;
    }

    public static matrix inverse(matrix A)
    {
        matrix A_inv = new matrix(A.size1);
        for (int n = 0; n < A.size1; n++)
        {
            vector e = new vector(A.size1);
            for (int i = 0; i < e.size; i++)
            {
                if (i == n)
                {
                    e[i] = 1;
                }
                else
                {
                    e[i] = 0;
                }
            }
            matrix L = decomp(A);
            vector x = solve(L, e);
            for (int i = 0; i < A_inv.size1; i++)
            {
                A_inv.set(i, n, x[i]);
            }
        }
        return A_inv;
    }
}
class main
    {
        public static int Main(string[] args)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-GB");
            int N = 0; //default values
            foreach (var arg in args)
            {
                var words = arg.Split(':');
                if (words[0] == "-size") N = int.Parse(words[1]);
            }
        if (N == 0)
        {
            System.Random rand = new System.Random();
            int n = rand.Next(0, 10);
            matrix A_first = new matrix(n);
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < n; k++)
                {
                    var random = rand.Next(0, 20);
                    A_first.set(i, k, random);
                }
            }
            matrix A = A_first.T * A_first;
            WriteLine("A. Test decomp method");
            WriteLine("     - The Cholesky-Banachiewicz algorithm is implemented for decomposition.");
            WriteLine("     - The deomposittion is tested on random real symmetric positive definite matrices");
            WriteLine("     - To create a random real symmetric positive difinite matrix first at random square real matrix B is generated. Then A = B^TB, and A is symmetric and positive definite.");
            WriteLine("     A:");
            A.print();
            matrix L = Cholesky.decomp(A);
            WriteLine("     L:");
            L.print();
            WriteLine($"     - Check that L is lower triangular using Cholesky.lower_triangular: {Cholesky.lower_triangular(L)}");
            WriteLine($"     - Check that L*L^T = A: {A.approx(L * L.T)}");
            WriteLine();
            WriteLine("B. Implementing linear equation solver, calculation of determinant and calculation of inverse matrix");
            WriteLine("     - Firstly the linear equation solver is implemented. A linear equation Ax=b can be rewritten to LL^Tx=b => Ly=b with L^Tx=y");
            WriteLine("     - As L is lower triangular and L^T therefore is upper triangular, Ly=b can be solved using forward substitution and L^Tx=y can afterwards be solved by back substitution.");
            WriteLine("     - Test of implemented solver on random real symmetric and positive definite matrix A and random matrix b:");
            int m = rand.Next(0, 10);
            matrix B = new matrix(m);
            for (int i = 0; i < m; i++)
            {
                for (int k = 0; k < m; k++)
                {
                    var random = rand.Next(0, 20);
                    B.set(i, k, random);
                }
            }
            matrix B_ny = B.T * B;
            matrix L_B = Cholesky.decomp(B_ny);
            vector b = new vector(m);
            for (int i = 0; i < m; i++)
            {
                b[i] = rand.Next(0, 100);
            }
            WriteLine(" A:");
            B_ny.print();
            WriteLine(" L:");
            L_B.print();
            WriteLine(" b: ");
            b.print();
            vector x = Cholesky.solve(L_B, b);
            WriteLine(" Computed x:");
            x.print();
            WriteLine($"     - Is Ax=b? {b.approx(B_ny * x)} ");
            WriteLine("     - The implemented solve method works as intended.");
            WriteLine();
            WriteLine("     - Secondly a method for calculation of the determinant is implemented.");
            WriteLine("     - det(A)=det(LL^T)=det(L)det(L^T)=det(L)^2. The determinant of L is the product of the diagonal elements, as it is lower triangular.");
            WriteLine("     - Testing the implemented method on random 3x3 matrix A=[[1, 2, 3],[0, 1, 4], [1, 0, 1]] with A^TA = [[2,2,4],[2,5,10],[4, 10, 26]]. det(A^TA) = 36");
            matrix A_det = new matrix(3);
            A_det[0, 0] = 1;
            A_det[0, 1] = 2;
            A_det[0, 2] = 3;
            A_det[1, 0] = 0;
            A_det[1, 1] = 1;
            A_det[1, 2] = 4;
            A_det[2, 0] = 1;
            A_det[2, 1] = 0;
            A_det[2, 2] = 1;
            WriteLine($"     - Is computed determiant for A^TA = 36? {matrix.approx(Cholesky.det(A_det.T * A_det), 36)} ");
            WriteLine();
            WriteLine("     - Lastly a method for calculating the inverse is implemented.");
            WriteLine("     - This is done by using the implemenbted solver to solve n linear equations Ax_i=e_i, where e_i is the i'th unit vector. x_i the make up the columns of the inverse matrix.");
            WriteLine("     - The implemented method is tested on a random square symmetric real positive definite matrix C:");
            int q = rand.Next(0, 10);
            matrix C = new matrix(q);
            for (int i = 0; i < q; i++)
            {
                for (int k = 0; k < q; k++)
                {
                    var random = rand.Next(0, 20);
                    C.set(i, k, random);
                }
            }
            matrix C_ny = C.T * C;
            C_ny.print();
            matrix ID2 = new matrix(C_ny.size1);
            for (int i = 0; i < ID2.size1; i++)
            {
                ID2.set(i, i, 1);
                for (int j = i + 1; j < ID2.size2; j++)
                {
                    ID2.set(i, j, 0);
                    ID2.set(j, i, 0);
                }
            }
            WriteLine("     - Calculated inverse of C is called B and is the following matrix:");
            matrix C_inv = Cholesky.inverse(C_ny);
            C_inv.print();
            WriteLine($"     - Is CB=I? {ID2.approx(C_ny * C_inv)}");
            WriteLine("     - As seen, the routine for calculation of the inverse works as intended.");
            WriteLine();
            WriteLine("C. Operations count");
            WriteLine("     - The measured time of a Cholesky decomposition as a function of matrix size N can be seen plotted in times.gnuplot.svg.");
            WriteLine("     - The function fitted to the time of a QR decomposition as a function of N is also plotted in times.gnuplot.svg for comparison.");
            WriteLine("     - This function was taken from homework 1, linear equations.");
            WriteLine("     - It is seen how the time scaling with N is a lot lower for Cholesky compared to QR decomposition. The scaling is still proportional to N^3 however.");
            WriteLine("     - It should be noted, however, that the QR decompositions in homework 1 was not performed specifically on symmetric  positive definite matrices, which is the case for the Cholesky decomposition.");

        }
        else
        {
            matrix ny = new matrix(N);
            var rnd = new System.Random();
            for (int i = 0; i < ny.size1; i++)
            {
                for (int t = 0; t < ny.size1; t++)
                {
                    ny[i, t] = 100.0 * (rnd.NextDouble() - 0.5);
                }
            }
            matrix Q3 = Cholesky.decomp(ny.T * ny);
        }
            return 0;
        }
    }
