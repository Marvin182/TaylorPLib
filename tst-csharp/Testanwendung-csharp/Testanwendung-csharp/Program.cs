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
            Console.WriteLine("Polynome: \n");
            showPolynom();
            Console.ReadKey();
            Console.Clear();
            Console.WriteLine("Matrizen mit Polynomen: \n");
            showMatrix();
            Console.ReadKey();
            Console.Clear();
            Console.WriteLine("Spezielle Matrizenmultiplikation: \n");
            showSpecialMatrix();
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

        private static void showMatrix()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[5, 5];
            Pa[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pa[0, 2] = new Polynomial(1, new double[] { 5, 6 });
            Pa[0, 3] = new Polynomial(1, new double[] { 1, 2 });
            Pa[0, 4] = new Polynomial(1, new double[] { 5, 6 });

            Pa[1, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[1, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pa[1, 2] = new Polynomial(1, new double[] { 5, 6 });
            Pa[1, 3] = new Polynomial(1, new double[] { 1, 2 });
            Pa[1, 4] = new Polynomial(1, new double[] { 5, 6 });

            Pa[2, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[2, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pa[2, 2] = new Polynomial(1, new double[] { 5, 6 });
            Pa[2, 3] = new Polynomial(1, new double[] { 1, 2 });
            Pa[2, 4] = new Polynomial(1, new double[] { 5, 6 });

            Pa[3, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[3, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pa[3, 2] = new Polynomial(1, new double[] { 5, 6 });
            Pa[3, 3] = new Polynomial(1, new double[] { 1, 2 });
            Pa[3, 4] = new Polynomial(1, new double[] { 5, 6 });

            Pa[4, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[4, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pa[4, 2] = new Polynomial(1, new double[] { 5, 6 });
            Pa[4, 3] = new Polynomial(1, new double[] { 1, 2 });
            Pa[4, 4] = new Polynomial(1, new double[] { 5, 6 });

            Polynomial[,] Pb = new Polynomial[5, 5];
            Pb[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[0, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pb[0, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pb[0, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[1, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pb[1, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pb[1, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[1, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[2, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[2, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pb[2, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pb[2, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[2, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[3, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[3, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pb[3, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pb[3, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[3, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[4, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[4, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pb[4, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pb[4, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[4, 4] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix M1 = new Matrix(5, 5, Pa);
            Matrix M2 = new Matrix(5, 5, Pb);
            Matrix M3 = M1 + M2;
            Matrix M4 = M1 * M2;

            Console.WriteLine("M1: \n" + M1.ToString());
            Console.WriteLine("M2: \n" + M2.ToString());
            Console.WriteLine("M1 + M2: \n" + M3.ToString());
            Console.WriteLine("M1 * M2: \n" + M4.ToString());
        }

        private static void showSpecialMatrix()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[3, 3];
            Pa[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pa[0, 2] = new Polynomial(1, new double[] { 5, 6 });

            Pa[1, 0] = new Polynomial(1, new double[] { 7, 8 });
            Pa[1, 1] = new Polynomial(1, new double[] { 9, 0 });
            Pa[1, 2] = new Polynomial(1, new double[] { 1, 2 });

            Pa[2, 0] = new Polynomial(1, new double[] { 3, 4 });
            Pa[2, 1] = new Polynomial(1, new double[] { 5, 6 });
            Pa[2, 2] = new Polynomial(1, new double[] { 7, 8 });

            Polynomial[,] Pb = new Polynomial[3, 3];
            Pb[0, 0] = new Polynomial(1, new double[] { 3, 4 });
            Pb[0, 1] = new Polynomial(1, new double[] { 5, 6 });
            Pb[0, 2] = new Polynomial(1, new double[] { 7, 8 });

            Pb[1, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pb[1, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pb[1, 2] = new Polynomial(1, new double[] { 5, 6 });

            Pb[2, 0] = new Polynomial(1, new double[] { 7, 8 });
            Pb[2, 1] = new Polynomial(1, new double[] { 9, 0 });
            Pb[2, 2] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] Pc = new Polynomial[3, 3];
            Pc[0, 0] = new Polynomial(1, new double[] { 7, 8 });
            Pc[0, 1] = new Polynomial(1, new double[] { 9, 0 });
            Pc[0, 2] = new Polynomial(1, new double[] { 1, 2 });

            Pc[1, 0] = new Polynomial(1, new double[] { 3, 4 });
            Pc[1, 1] = new Polynomial(1, new double[] { 5, 6 });
            Pc[1, 2] = new Polynomial(1, new double[] { 7, 8 });

            Pc[2, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pc[2, 1] = new Polynomial(1, new double[] { 3, 4 });
            Pc[2, 2] = new Polynomial(1, new double[] { 5, 6 });



            #endregion

            Matrix M1 = new Matrix(3, 3, Pa);
            Matrix M2 = new Matrix(3, 3, Pb);
            Matrix M3 = new Matrix(3, 3, Pc);

            Console.WriteLine("( M1 * M2 * alpha ) + ( M3 * beta ) where: ");
            Console.WriteLine("M1: \n" + M1.ToString());
            Console.WriteLine("M2: \n" + M2.ToString());
            Console.WriteLine("M3: \n" + M3.ToString());
            Console.WriteLine("alpha: 2");
            Console.WriteLine("beta: 2\n");
            DateTime start = DateTime.Now;
            String temp = ( (M1 * M2 * 2) + (M3 * 2) ).ToString();
            Console.WriteLine( temp );
            DateTime end = DateTime.Now;
            Console.WriteLine("Time needed: " + (end - start).ToString());
            Console.WriteLine("This was normal mode...");
            start = DateTime.Now;
            M3.mmCaABbC(2, 2, M1, M2);
            String temp2 = M3.ToString();
            Console.WriteLine(temp2);
            end = DateTime.Now;
            Console.WriteLine("Time needed: " + (end - start).ToString());
            Console.WriteLine("This was method mmCaABbC...");

            if (temp.Equals(temp2))
                Console.WriteLine("Both Are Equal...");
            else
                Console.WriteLine("Both Are NOT Equal. Something went wrong...");
        }
    }
}
