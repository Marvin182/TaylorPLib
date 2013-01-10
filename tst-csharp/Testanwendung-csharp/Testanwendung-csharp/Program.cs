using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using LibMatrix;

namespace Testanwendung_csharp
{
    class Program
    {
        static void Main(string[] args)
        {
            showPolynom();
            Console.ReadKey();
        }

        private static void showPolynom()
        {
            Polynomial P1 = new Polynomial(2, new double[] { 2, 4, 3 });
            Polynomial P2 = new Polynomial(2, new double[] { 2, 4, 3 });

            Polynomial P3 = P1 + P2;
            Polynomial P4 = P1 * P2;
            Polynomial P5 = P4 / P2;

            Console.WriteLine("P1: \n" + P1.ToString());
            Console.WriteLine("P2: \n" + P2.ToString());
            Console.WriteLine("P1 + P2: \n" + P3.ToString());
            Console.WriteLine("P1 * P2: \n" + P4.ToString());
            Console.WriteLine("P1 * P2 / P2: \n" + P5.ToString());
        }
    }
}
