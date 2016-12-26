using System;
using Methods;

namespace Lab5
{
    class Program
    {
        static void Main(string[] args)
        {
            RungeKuttaODU2 RK_ODU2 = new RungeKuttaODU2(0, 2, 0.2, 0.2 * Math.Pow(Math.E, 4), func, funcX);
            RK_ODU2.Output();

            Console.ReadKey();
        }

        private static double func(double x)
        {
            return Math.Pow(Math.E, -(x * x)) * (x * x / 2);
        }

        private static double funcX(double x, double y)
        {
            return x * Math.Pow(Math.E, -(x * x)) - 2 * x * y;
        }
    }
}
