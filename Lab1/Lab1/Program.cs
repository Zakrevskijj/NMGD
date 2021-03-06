﻿using System;
using Methods;

namespace Lab1
{
    class Program
    {
        static void Main(string[] args)
        {//x=0..2 y(0)=0
            RungeKutta RK = new RungeKutta(0.00001, 2, 0, func, funcX);
            RK.Output();

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
