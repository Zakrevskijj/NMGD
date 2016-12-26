using System;
using Methods;

namespace Lab4
{
    class Program
    {
        static void Main(string[] args)
        {
            RungeKuttaODU2 rk = new RungeKuttaODU2(0, 2, 0, 0, func, funcX);
            rk.Output();

            Console.ReadKey();
        }

        private static double func(double x)
        {
            return x * x + x * x * x * x;
        }

        private static double funcX(double x, double y)
        {
            return (2 * x * x * x) + (2 * y / x);
        }
    }
}
