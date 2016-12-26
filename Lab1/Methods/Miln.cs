using System;

namespace Methods
{
    public class Miln
    {
        double[] X, Ym4, Ym6, Y;
        double eps = 10E-7;
        int N;
        Func<double, double> func;
        Func<double, double, double> fx;

        public Miln(double x0, double xn, double y0, Func<double, double> func, Func<double, double, double> fx)
        {
            N = 100;
            double h = (xn - x0) / N;
            this.N += 1;
            X = new double[this.N];
            Ym4 = new double[this.N];
            Ym6 = new double[this.N];
            Y = new double[this.N];
            this.func = func;
            this.fx = fx;

            Y[0] = y0;
            Ym4[0] = y0;
            Ym6[0] = y0;

            for (int i = 0; i < this.N; i++)
            {
                X[i] = x0 + h * i;
                Y[i] = func(X[i]);
            }


            Calculate();
        }

        private void Calculate()
        {
            Ym4[1] = RungeKutta.RK4(X[0], X[1], Ym4[0], eps, fx);
            Ym4[2] = RungeKutta.RK4(X[1], X[2], Ym4[1], eps, fx);
            Ym4[3] = RungeKutta.RK4(X[2], X[3], Ym4[2], eps, fx);
            Ym6[1] = RungeKutta.RK4(X[0], X[1], Ym6[0], eps, fx);
            Ym6[2] = RungeKutta.RK4(X[1], X[2], Ym6[1], eps, fx);
            Ym6[3] = RungeKutta.RK4(X[2], X[3], Ym6[2], eps, fx);
            Ym6[4] = RungeKutta.RK4(X[3], X[4], Ym6[3], eps, fx);
            for (int i = 0; i < N - 4; i++)
            {
                Ym4[i + 4] = M4(X[i], X[i + 1], X[i + 2], X[i + 3], X[i + 4], Ym4[i], Ym4[i + 1], Ym4[i + 2], Ym4[i + 3], fx);
            }
            for (int i = 0; i < N - 5; i++)
            {
                Ym6[i + 4] = M6(X[i], X[i + 1], X[i + 2], X[i + 3], X[i + 4], X[i + 5], Ym6[i], Ym6[i + 1], Ym6[i + 2], Ym6[i + 3], Ym6[i + 4], fx);
            }
        }

        public static double M4(double x0, double x1, double x2, double x3, double x4, double y0, double y1, double y2, double y3, Func<double, double, double> fx)
        {
            double h = x1 - x0;
            double y = y0 + 4 * h / 3 * (2 * fx(x2, y2)) - fx(x1, y1) + 2 * fx(y0, x0);
            for (int i = 0; i < 4; i++)
            {
                y = y2 + h / 3 * (y + 4 * fx(x3, y3) + fx(x2, y2));
            }
            return y2 + h / 3 * (fx(x4, y) + 4 * fx(x3, y3) + fx(x2, y2));
        }

        public static double M6(double x0, double x1, double x2, double x3, double x4, double x5, double y0, double y1, double y2, double y3, double y4, Func<double, double, double> fx)//Mistake in lecture
        {
            double h = x1 - x0;
            double y = y3 + 3 * h / 10 * (11 * fx(y4, x4) - 14 * fx(x3, y3) + 26 * fx(x2, y2) - 14 * fx(x1, y1) + 11 * fx(x0, y0));
            for (int i = 0; i < 6; i++)
            {
                y = y1 + 4 * h / 90 * (7 * fx(x5, y) + 32 * fx(x4, y4) + 12 * fx(x3, y3) + 32 * fx(x2, y2) + 7 * fx(x1, y1));
            }
            return y1 + 4 * h / 90 * (7 * fx(x5, y) + 32 * fx(x4, y4) + 12 * fx(x3, y3) + 32 * fx(x2, y2) + 7 * fx(x1, y1));
        }

        public void Output()
        {
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    y4 = {2:F10}    y6 = {3:F10}", X[i], Y[i], Ym4[i], Ym6[i]);
            }
            Console.WriteLine("Finished");
        }
    }

}
