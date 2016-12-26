using System;

namespace Methods
{
    public class RungeKutta
    {
        double[] X, Ye, YRK3, YRK4, Y;
        double eps = 10E-8;
        int N;
        Func<double, double> func;
        Func<double, double, double> fx;

        public RungeKutta(double x0, double xn, double y0, Func<double, double> func, Func<double, double, double> fx, int N = 100)
        {
            this.N = N;
            double h = (xn - x0) / N;
            this.N += 1;
            X = new double[this.N];
            Ye = new double[this.N];
            YRK3 = new double[this.N];
            YRK4 = new double[this.N];
            Y = new double[this.N];
            Y[0] = y0;
            YRK3[0] = y0;
            YRK4[0] = y0;
            Ye[0] = y0;
            this.func = func;
            this.fx = fx;

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
                Ye[i + 1] = Eylera(X[i], X[i + 1], Ye[i], eps, fx);
                YRK3[i + 1] = RK3(X[i], X[i + 1], YRK3[i], eps, fx);
                YRK4[i + 1] = RK4(X[i], X[i + 1], YRK4[i], eps, fx);

            }
        }

        public static double Eylera(double x0, double x1, double y0, double eps, Func<double, double, double> fx)
        {
            double hx, xj, yj, ym, err; int j, m;
            m = (int)(x1 - x0) + 1; ym = double.MaxValue;
            do
            {
                hx = (x1 - x0) / m; yj = y0; xj = x0;
                for (j = 0; j < m; j++)
                {
                    yj = yj + hx * fx(xj, yj); xj += hx;
                }
                err = Math.Abs(ym - yj); ym = yj; m *= 10;
            }
            while (err > eps);
            return ym;
        }

        public static double RK3(double x0, double x1, double y0, double eps, Func<double, double, double> fx)
        {
            double ht, xj, yj, ym, ymdx, err;
            int m, j;
            double k1, k2, k3;

            m = (int)(x1 - x0) + 1; ym = double.MaxValue; ymdx = ym;
            do
            {
                ht = (x1 - x0) / m;
                xj = x0; yj = y0;
                for (j = 0; j < m; j++)
                {
                    k1 = ht * fx(xj, yj);
                    k2 = ht * fx(xj + 0.5 * ht, yj + 0.5 * k1);
                    k3 = ht * fx(xj + ht, yj + 2 * k2 - k1);
                    yj = yj + (k1 + 4 * k2 + k3) / 6;
                    xj += ht;
                }
                err = (Math.Abs(ym - yj));
                ym = yj; m *= 10;
            } while (err > eps);
            x1 = ym;

            return ym;
        }

        public static double RK4(double x0, double x1, double y0, double eps, Func<double, double, double> fx)
        {
            double ht, xj, yj, ym, ymdx, err;
            int m, j;
            double k1, k2, k3, k4;

            m = (int)(x1 - x0) + 1; ym = double.MaxValue; ymdx = ym;
            do
            {
                ht = (x1 - x0) / m;
                xj = x0; yj = y0;
                for (j = 0; j < m; j++)
                {
                    k1 = ht * fx(xj, yj);
                    k2 = ht * fx(xj + 0.5 * ht, yj + 0.5 * k1);
                    k3 = ht * fx(xj + 0.5 * ht, yj + 0.5 * k2);
                    k4 = ht * fx(xj + ht, yj + k3);
                    yj = yj + 1.0 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
                    xj += ht;
                }
                err = (Math.Abs(ym - yj));
                ym = yj; m *= 10;
            } while (err > eps);
            x1 = ym;

            return ym;
        }

        public void Output()
        {
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    yEyl = {2:F10}    yRK4 = {3:F10}", X[i], Y[i], Ye[i], YRK4[i]);
            }

            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    RK3 = {2:F10}", X[i], Y[i], YRK3[i]);
            }

            Console.WriteLine("Finished");
        }
    }
}
