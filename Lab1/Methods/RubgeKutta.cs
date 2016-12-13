using System;

namespace Methods
{
    public class RungeKutta {
        readonly double[] _x;
        readonly double[] _ye;
        readonly double[] _yrk3;
        readonly double[] _yrk4;
        readonly double[] _y;
        readonly double _eps = 10E-7;
        readonly int _n;

        public RungeKutta(double x0, double xn, Func<double, double> func, Func<double, double, double> fx) {
            _n = 100;
            double h = (xn - x0) / _n;
            _x = new double[_n];
            _ye = new double[_n];
            _yrk3 = new double[_n];
            _yrk4 = new double[_n];
            _y = new double[_n];

            for (int i = 0; i < _n; i++) {
                _x[i] = x0 + h * i;
                _y[i] = func(_x[i]);
            }

            Calculate(fx);
        }

        private void Calculate(Func<double, double, double> fx) {
            for (int i = 0; i < _n - 1; i++) {
                _ye[i + 1] = Eylera(_x[i], _x[i + 1], _ye[i], _eps, fx);
                _yrk3[i + 1] = Rk3(_x[i], _x[i + 1], _yrk3[i], _eps, fx);
                _yrk4[i + 1] = Rk4(_x[i], _x[i + 1], _yrk4[i], _eps, fx);
            }
        }

        private static double Eylera(double x0, double x1, double y0, double eps, Func<double, double, double> fx) {
            double err;
            var m = (int)(x1 - x0) + 1;
            var ym = double.MaxValue;

            do {
                var hx = (x1 - x0) / m;
                var yj = y0;
                var xj = x0;
                int j;
                for (j = 0; j < m; j++) {
                    yj = yj + hx * fx(xj, yj);
                    xj += hx;
                }
                err = Math.Abs(ym - yj);
                ym = yj;
                m *= 10;
            }
            while (err > eps);

            return ym;
        }

        private static double Rk3(double x0, double x1, double y0, double eps, Func<double, double, double> fx) {
            double err;

            var m = (int)(x1 - x0) + 1;
            var ym = double.MaxValue;
            do {
                var ht = (x1 - x0) / m;
                var xj = x0;
                var yj = y0;
                int j;
                for (j = 0; j < m; j++) {
                    var k1 = fx(xj, yj);
                    var k2 = fx(xj + 0.5 * ht, yj + 0.5 * k1);
                    var k3 = fx(xj + ht, yj + 2 * k2 - k1);
                    yj = yj + ht * (k1 + 4 * k2 + k3) / 6;
                    xj += ht;
                }
                err = (Math.Abs(ym - yj));
                ym = yj;
                m *= 10;
            } while (err > eps);

            return ym;
        }

        private static double Rk4(double x0, double x1, double y0, double eps, Func<double, double, double> fx) {
            double err;

            var m = (int)(x1 - x0) + 1;
            var ym = double.MaxValue;
            do {
                var ht = (x1 - x0) / m;
                var xj = x0;
                var yj = y0;
                int j;
                for (j = 0; j < m; j++) {
                    var k1 = fx(xj, yj);
                    var k2 = fx(xj + 0.5 * ht, yj + 0.5 * k1);
                    var k3 = fx(xj + 0.5 * ht, yj + 0.5 * k2);
                    var k4 = fx(xj + ht, yj + k3);
                    yj = yj + ht / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
                    xj += ht;
                }
                err = (Math.Abs(ym - yj));
                ym = yj;
                m *= 10;
            } while (err > eps);

            return ym;
        }

        public void Output() {
            for (int i = 0; i < _n; i++) {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    yEyl = {2:F10}    yRK4 = {2:F10}", _x[i], _y[i], _ye[i]);
            }

            for (int i = 0; i < _n; i++) {
                Console.WriteLine("x = {0:F6}    y = {1:F10}    RK3 = {2:F10}", _x[i], _y[i], _yrk3[i]);
            }

            Console.Write("its All!");
        }
    }
}
