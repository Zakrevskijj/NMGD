using System;

namespace Methods
{
    public class RungeKuttaODU2
    {
        double[] X, YRK4, YdxRK4, Y;
        double eps = 10E-8;
        int N;
        bool full;
        Func<double, double> _func;
        Func<double, double, double> fx;
        Func<double, double, double, double> fxx;


        public RungeKuttaODU2(double x0, double xn, double y0, double ydx0, Func<double, double> func, Func<double, double, double> fx, int N = 100)
        {
            this.N = N;
            double h = (xn - x0) / N;
            this.N += 1;
            X = new double[this.N];
            YRK4 = new double[this.N];
            YdxRK4 = new double[this.N];
            Y = new double[this.N];
            Y[0] = y0;
            YRK4[0] = y0;
            YdxRK4[0] = ydx0;
            this._func = func;
            this.fx = fx;
            full = false;
            for (int i = 0; i < N; i++)
            {
                X[i] = x0 + h * i;
                Y[i] = func(X[i]);
            }
            Calculate();
        }

        public RungeKuttaODU2(double x0, double xn, double y0, double ydx0, Func<double, double> func, Func<double, double, double, double> fxx, int N = 100)
        {
            this.N = N;
            double h = (xn - x0) / N;
            this.N += 1;
            X = new double[this.N];
            YRK4 = new double[this.N];
            YdxRK4 = new double[this.N];
            Y = new double[this.N];
            Y[0] = y0;
            YRK4[0] = y0;
            YdxRK4[0] = ydx0;
            this._func = func;
            this.fxx = fxx;
            full = true;
            for (int i = 0; i < this.N; i++)
            {
                X[i] = x0 + h * i;
                Y[i] = func(X[i]);
            }
            Calculate();
        }

        private void Calculate()
        {
            for (int i = 0; i < N - 1; i++)
            {
                if (full)
                {
                    RK4(X[i], X[i + 1], YRK4[i], YdxRK4[i], eps, fxx, out YRK4[i + 1], out YdxRK4[i + 1]);
                }
                else RK4_2(X[i], X[i + 1], YRK4[i], YdxRK4[i], eps, fx, out YRK4[i + 1], out YdxRK4[i + 1]);
            }
        }

        public static void RK4(double x0, double x1, double y0, double ydx0, double eps, Func<double, double, double, double> fx, out double y1, out double ydx1)
        {
            double ht, xj, yj, ym, yjdx, ymdx, err;
            int m, j;
            double k1, k2, k3, k4;

            m = (int)(x1 - x0) + 1; ym = double.MaxValue; ymdx = double.MaxValue;
            do
            {
                ht = (x1 - x0) / m;
                xj = x0; yj = y0; yjdx = ydx0;
                for (j = 0; j < m; j++)
                {
                    k1 = ht * fx(xj, yj, yjdx);
                    k2 = ht * fx(xj + 0.5 * ht, yj + 0.5 * ht * yjdx + ht / 8 * k1, yjdx + 0.5 * k1);
                    k3 = ht * fx(xj + 0.5 * ht, yj + 0.5 * ht * yjdx + ht / 8 * k1, yjdx + 0.5 * k2);
                    k4 = ht * fx(xj + ht, yj + ht * yjdx + ht / 2 * k3, yjdx + k3);
                    yj = yj + ht * (yjdx + (k1 + k2 + k3) / 6);
                    yjdx = yjdx + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                    xj += ht;
                }
                err = 0.5 * (Math.Abs(ym - yj) + Math.Abs(ymdx - yjdx));
                ym = yj; ymdx = yjdx; m *= 10;
            } while (err > eps);

            y1 = ym;
            ydx1 = ymdx;
        }

        public static void RK4_2(double x0, double x1, double y0, double ydx0, double eps, Func<double, double, double> fx, out double y1, out double ydx1)
        {
            double ht, xj, yj, ym, yjdx, ymdx, err;
            int m, j;
            double k1, k2, k3;

            m = (int)(x1 - x0) + 1; ym = double.MaxValue; ymdx = double.MaxValue;
            do
            {
                ht = (x1 - x0) / m;
                xj = x0; yj = y0; yjdx = ydx0;
                for (j = 0; j < m; j++)
                {
                    k1 = ht * fx(xj, yj);
                    k2 = ht * fx(xj + 0.5 * ht, yj + 0.5 * ht * yjdx + k1 / 8);
                    k3 = ht * fx(xj + ht, yj + ht * yjdx + ht / 2 * k2);
                    yj = yj + ht * (yjdx + (k1 + 2 * k2) / 6);
                    yjdx = yjdx + (k1 + 4 * k2 + k3) / 6;
                    xj += ht;
                }
                err = 0.5 * (Math.Abs(ym - yj) + Math.Abs(ymdx - yjdx));
                ym = yj; ymdx = yjdx; m *= 10;
            } while (err > eps);

            y1 = ym;
            ydx1 = ymdx;
        }

        public void Output()
        {
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    yRK4 = {2:F10}", X[i], Y[i], YRK4[i]);
            }
            Console.WriteLine("Finished");
        }
    }
}
