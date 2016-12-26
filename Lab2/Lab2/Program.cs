using System;
using Methods;

namespace Lab2
{
    class Program
    {
        static void Main(string[] args)
        {
            Adams A = new Adams(0, 2, 0, func, funcX);
            A.Output();

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
