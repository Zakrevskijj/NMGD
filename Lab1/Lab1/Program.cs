using System;
using Methods;

namespace Lab1 {
    class Program {
        static void Main(string[] args) {//x=0..2 y(0)=0
            RungeKutta RK = new RungeKutta(0, 2, func, funcX);
            RK.Output();
            Console.ReadKey();
        }

        private static double func(double x) {
            return Math.Pow(Math.E, -(x * x)) * (x * x / 2);
        }

        private static double funcX(double x, double y) {
            return x * Math.Pow(Math.E, -(x * x)) - 2 * x * y;
        }
    }
}
