using LegendreFunction;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace Term_Papper_GA
{
    public class ValueContainer
    {
        Dictionary<(double p1, double p2, double p3), (double[][] b_km, double f)> XDictionry;

        public double[] T { get; }
        public double[] A { get; }
        public double Mu { get; }
        public int N_ { get; }

        public (double p1, double p2, double p3)[] X { get; }
        public (double p1, double p2, double p3)[] Y { get; }
        public (double p1, double p2, double p3)[] X_test { get; }

        public double K((double p1, double p2, double p3) x, (double p1, double p2, double p3) y)  ///todo
        {
            return 1/Math.Sqrt(Math.Pow(x.p1-y.p1,2)+ Math.Pow(x.p2 - y.p2, 2)+ Math.Pow(x.p3 - y.p3, 2));
        }
        public double F((double p1, double p2, double p3) x)  ///todo
        {
            //1
            return 10;
            //return Math.Sqrt(Math.Pow(x.p1,2) + Math.Pow(x.p2, 2) + Math.Pow(x.p3, 2));
        }



        public ValueContainer(int n, int n_, int number_x_theta, int number_x_phi)
        {
            X_test = new (double p1, double p2, double p3)[] { (0, 0, 0), (0.15, 0.15, 0.15), (0, 0.5, 0), (-0.1, -0.1, -0.1), (-0.2, 0.2, 0) ,(0,0,1)};

            N_ = n_;

            ////////////////////////////////////////a_s begin
            double[] c = new double[N_ + 2];
            for (int k = 0; k < (N_ + 1) / 2 + 1; k++)
            {
                c[2 * k] = Math.Pow((-1), k) * Factorial(2 * (N_ + 1) - 2 * k) / Factorial((N_ + 1) - 2 * k) / Factorial(N_ + 1 - k) / Factorial(k) / Math.Pow(2, N_ + 1);
                if (k != (N_ + 1) / 2)
                {
                    c[2 * k + 1] = 0;
                }
            }

            var t1 = FindRoots.Polynomial(c.Reverse().ToArray()).Select(el => el.Real).OrderBy(el => el).ToArray();  //octaive sort(roots(c))
            T = new double[n_ + 1];
            for (int i = 0; i < n_ + 1; i++)
            {
                T[i] = MathNet.Numerics.RootFinding.Bisection.FindRootExpand((x) => LegendrePolinom(n_ + 1, 0, x), t1[i] - 0.01, t1[i] + 0.01); //fsolve
            }

            A = new double[n_ + 1];
            for (int i = 0; i < n_ + 1; i++)
            {
                A[i] = 2 * (1 - Math.Pow(T[i], 2)) / Math.Pow( (n_ + 1) * LegendrePolinom(n_, 0, T[i]) , 2);
            }
            //////////////////////////////////////////a_s end

            //////////////////////////////////////////mu_p begin
            Mu = Math.PI / (n_ + 1);
            //////////////////////////////////////////mu_p end

            /////////////////////////////////////////hat(x) begin
            double h_theta = Math.PI / (number_x_theta - 1);  // [0]0..[last]pi  зверху вниз
            double[] x_theta = new double[number_x_theta];
            for (int i = 0; i < number_x_theta; i++)
            {
                x_theta[i] = i * h_theta;
            }

            double h_phi = 2 * Math.PI / (number_x_phi);  // [0]0..[last](ne 2pi) бо 0 = 2pi зліва направо
            double[] x_phi = new double[number_x_phi];
            for (int i = 0; i < number_x_phi; i++)
            {
                x_phi[i] = i * h_phi;
            }

            
            X = new (double p1, double p2, double p3)[number_x_theta * number_x_phi];
            for (int i = 0; i < number_x_theta; i++)
            {
                for (int j = 0; j < number_x_phi; j++)
                {
                    X[i * number_x_phi + j] = Get_X(x_theta[i], x_phi[j]);
                    X[i * number_x_phi + j] =
                        (
                        Math.Round(X[i * number_x_phi + j].p1, 10),
                        Math.Round(X[i * number_x_phi + j].p2, 10),
                        Math.Round(X[i * number_x_phi + j].p3, 10)
                        );  // зверху і внизу по 10 однакових точок (заукруглення потрібні при створенні словника далі)
                }
            }

            /////////////////////////////////////////hat(x) end

            /////////////////////////////////////////hat(y) begin
            int s_length = n_ + 1;   // 1...n'+1   цикл пройде по  0 ..n'
            int p_length = (2 * n_ + 1) + 1; // 0...2n'+1  цикл пройде по 0 .. 2n'+1

            Y = new (double p1, double p2, double p3)[s_length * p_length];
            for (int s = 0; s < s_length; s++)
            {
                for (int p = 0; p < p_length; p++)
                {
                    Y[s * p_length + p] = Get_Y(s, p);         //Y_sp = Y[s * p_length + p];
                }
            }

            /////////////////////////////////////////hat(y) end

            XDictionry = new Dictionary<(double p1, double p2, double p3), (double[][] b, double f)>();

            foreach (var x in X)
            {
                if (!XDictionry.ContainsKey(x))  // на верхній і нижній точках по 10 однакових точок = рахуємо їх як 1 точку
                {
                    double[][] b = new double[n + 1][];
                    for (int k = 0; k < n + 1; k++)
                    {
                        b[k] = new double[2 * k + 1];
                        for (int m = -k; m < k + 1; m++)
                        {
                            b[k][m + k] = B_km(x, s_length, p_length, k, m);
                        }
                    }
                    XDictionry.Add(x, (b, F(x)));
                }
            }


            //print як в octaive
            //Console.WriteLine($"n'= {N_}");

            //Console.Write("c: ");
            //for (int i = 0; i < c.Length; i++)
            //{
            //    Console.Write($"{c[i]} ");
            //}
            //Console.WriteLine();

            //Console.Write("t1: ");
            //for (int i = 0; i < t1.Length; i++)
            //{
            //    Console.Write($"{t1[i]} ");
            //}
            //Console.WriteLine();

            //Console.Write("t: ");
            //for (int i = 0; i < T.Length; i++)
            //{
            //    Console.Write($"{T[i]} ");
            //}
            //Console.WriteLine();

            //Console.Write("a: ");
            //for (int i = 0; i < A.Length; i++)
            //{
            //    Console.Write($"{A[i]} ");
            //}
            //Console.WriteLine();

            //Console.WriteLine($"Y[] ={Y.Length}");
            //for (int i = 0; i < Y.Length; i++)
            //{
            //    Console.WriteLine($"{i+1}: {Y[i].p1}, {Y[i].p2}, {Y[i].p3}");
            //}

        }



        public double GetEval(double[][] psi)  ////add regularyzation
        {
            //double[] evals = new double[XDictionry.Count];
            //int index = 0;
            //foreach (var d in XDictionry)
            //{
            //    evals[index] = 0;
            //    for (int k = 0; k < psi.Length; k++)
            //    {
            //        for (int m = 0; m < psi[k].Length; m++)
            //        {
            //            evals[index] += d.Value.b_km[k][m] * psi[k][m];
            //        }
            //    }
            //    evals[index] = Math.Abs(evals[index] - d.Value.f);

            //    index++;
            //}
            //return evals.Max();

            double alpha = 0.001;
            double[] evals = new double[XDictionry.Count];
            int index = 0;
            foreach (var d in XDictionry)
            {
                double psi_reg = 0;
                evals[index] = 0;
                for (int k = 0; k < psi.Length; k++)
                {
                    for (int m = 0; m < psi[k].Length; m++)
                    {
                        evals[index] += d.Value.b_km[k][m] * psi[k][m];
                        psi_reg += Math.Pow(psi[k][m], 2);
                    }
                }
                evals[index] = Math.Pow(evals[index] - d.Value.f, 2) + alpha * psi_reg;

                index++;
            }
            return evals.Max();
        }
        private double B_km((double p1, double p2, double p3) x, int s_length, int p_length, int k, int m)
        {
            double b_km = 0;
            for (int s = 0; s < s_length; s++)
            {
                for (int p = 0; p < p_length; p++)
                {
                    b_km += Mu * A[s] * K(x, Y[s * p_length + p]) * GetY_km(k, m, s, p);
                }
            }
            return b_km;
        }

        private long Factorial(int a) //12! ok. 13! not ok.
        {
            if (a == 0)
                return 1;
            else
                return a * Factorial(a - 1);
        }
        private double LegendrePolinom(int n, int m, double x)
        {
            return Legendre.ALegendFnmSeri(n, x)[m];
        }
        private (double p1, double p2, double p3) Get_X(double theta, double phi)
        {
            return (Math.Sin(theta) * Math.Cos(phi), Math.Sin(theta) * Math.Sin(phi), Math.Cos(theta));
        }
        private (double p1, double p2, double p3) Get_Y(int s, int p)
        {       
            double th = Theta(s);
            double ph = Phi(p);
            return (Math.Sin(th) * Math.Cos(ph), Math.Sin(th) * Math.Sin(ph), Math.Cos(th));
        }
        private double Theta(int s)
        {
            return Math.Acos(T[s]);
        }
        private double Phi(int p)
        {
            return Math.PI * p / (N_ + 1);
        }
        public double GetY_km(int k, int m, int s, int p)
        {
            double c_km = Math.Pow(-1, ((Math.Abs(m) - m) / 2)) * Math.Sqrt((2 * k + 1) / (4 * Math.PI) * Factorial(k - Math.Abs(m)) / Factorial(k + Math.Abs(m)));

            double cosi1 = 0;
            if (m < 0)
            {
                cosi1 = Math.Cos(Math.Abs(m) * Phi(p));

            }
            else if (m > 0)
            {
                cosi1 = Math.Sin(Math.Abs(m) * Phi(p));
            }
            else
            {
                cosi1 = 1;
            }
            return c_km * LegendrePolinom(k, Math.Abs(m), Math.Cos(Theta(s))) * cosi1;
        }

        public double[] Test(int n_, double[][] psi)
        {
            int s_length = n_ + 1;
            int p_length = 2 * (n_ + 1);

            double[] sum = new double[X_test.Length];
            int i = 0;
            foreach (var x in X_test)
            {
                sum[i] = 0;
                for (int s = 0; s < s_length; s++)
                {
                    for (int p = 0; p < p_length; p++)
                    {
                        sum[i] += Mu * A[s] * K(x, Y[s * p_length + p]) * psi_y(s, p, psi);
                    }
                }
                sum[i] = Math.Abs(sum[i] - F(x));
                i++;
            }

            return sum;
        }
        private double psi_y(int s,int p, double[][] psi)
        {
            double sum = 0;
            for (int k = 0; k < psi.Length; k++)
            {
                for (int m = -k; m < k + 1; m++)
                {
                    sum += psi[k][m + k] * GetY_km(k, m, s, p);
                }
            }
            return sum;
        }
    }
}

