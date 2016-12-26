using System;

namespace Methods
{
    public class Adams
    {
        double[] X, Ye, Yi, Y;
        double eps = 10E-7;
        int N;
        Func<double, double> func;
        Func<double, double, double> fx;

        public Adams(double x0, double xn, double y0, Func<double, double> func, Func<double, double, double> fx)
        {
            N = 100;
            double h = (xn - x0) / N;
            this.N += 1;
            X = new double[this.N];
            Ye = new double[this.N];
            Yi = new double[this.N];
            Y = new double[this.N];
            this.func = func;
            this.fx = fx;

            Y[0] = y0;
            Ye[0] = y0;
            Yi[0] = y0;

            for (int i = 0; i < this.N; i++)
            {
                X[i] = x0 + h * i;
                Y[i] = func(X[i]);
            }


            Calculate();
        }

        private void Calculate()
        {
            Ye[1] = RungeKutta.RK4(X[0], X[1], Ye[0], eps, fx);
            Ye[2] = RungeKutta.RK4(X[1], X[2], Ye[1], eps, fx);
            Ye[3] = RungeKutta.RK4(X[2], X[3], Ye[2], eps, fx);
            Yi[1] = RungeKutta.RK4(X[0], X[1], Yi[0], eps, fx);
            Yi[2] = RungeKutta.RK4(X[1], X[2], Yi[1], eps, fx);
            Yi[3] = RungeKutta.RK4(X[2], X[3], Yi[2], eps, fx);
            for (int i = 0; i < N - 4; i++)
            {
                Ye[i + 4] = ecsp(X[i], X[i + 1], X[i + 2], X[i + 3], Ye[i], Ye[i + 1], Ye[i + 2], Ye[i + 3], fx);
                Yi[i + 4] = ecsp(X[i], X[i + 1], X[i + 2], X[i + 3], Yi[i], Yi[i + 1], Yi[i + 2], Yi[i + 3], fx);
                Yi[i + 4] = intr(X[i], X[i + 1], X[i + 2], X[i + 3], Yi[i + 1], Yi[i + 2], Yi[i + 3], Yi[i + 4], fx);
            }
        }

        public static double ecsp(double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3, Func<double, double, double> fx)
        {
            double h = x1 - x0;
            return y3 + h / 24 * (55 * fx(x3, y3) - 59 * fx(x2, y2) + 37 * fx(x1, y1) - 9 * fx(x0, y0));
        }

        public static double intr(double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3, Func<double, double, double> fx)
        {
            double h = x1 - x0;
            return y2 + h / 24 * (9 * fx(x3, y3) + 19 * fx(x2, y2) - 5 * fx(x1, y1) + fx(x0, y0));
        }


        public void Output()
        {
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    yE = {2:F10}    yI = {3:F10}", X[i], Y[i], Ye[i], Yi[i]);
            }
            Console.WriteLine("Finished");
        }
    }

}
