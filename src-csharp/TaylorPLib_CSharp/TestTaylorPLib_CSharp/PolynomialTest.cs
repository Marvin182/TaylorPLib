using LibMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.IO;

namespace TestTaylorPLib_CSharp
{
    
    /// <summary>
    ///Dies ist eine Testklasse für "PolynomialTest" und soll
    ///alle PolynomialTest Komponententests enthalten.
    ///</summary>
    [TestClass()]
    public class PolynomialTest
    {
        private TestContext testContextInstance;

        /// <summary>
        ///Ruft den Testkontext auf, der Informationen
        ///über und Funktionalität für den aktuellen Testlauf bietet, oder legt diesen fest.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Zusätzliche Testattribute
        // 
        //Sie können beim Verfassen Ihrer Tests die folgenden zusätzlichen Attribute verwenden:
        //
        //Mit ClassInitialize führen Sie Code aus, bevor Sie den ersten Test in der Klasse ausführen.
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Mit ClassCleanup führen Sie Code aus, nachdem alle Tests in einer Klasse ausgeführt wurden.
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Mit TestInitialize können Sie vor jedem einzelnen Test Code ausführen.
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Mit TestCleanup können Sie nach jedem einzelnen Test Code ausführen.
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion

        #region Constructor Tests
        /// <summary>
        ///Ein Test für "Polynomial-Konstruktor"
        ///</summary>
        [TestMethod()]
        public void PolynomialConstructorTest()
        {
            Polynomial_Accessor target = new Polynomial_Accessor(0, new double[] {1}, 1);
            Polynomial_Accessor expected = new Polynomial_Accessor();
            Assert.AreEqual(target._order, expected._order);
            CollectionAssert.AreEqual(target._coeffs, expected._coeffs);
            Assert.AreEqual(target._constant, expected._constant);
        }

        /// <summary>
        ///Ein Test für "Polynomial-Konstruktor"
        ///</summary>
        [TestMethod()]
        public void PolynomialConstructorTest1()
        {
            Polynomial_Accessor target = new Polynomial_Accessor();
            Assert.AreEqual(target._order, 0);
            CollectionAssert.AreEqual(target._coeffs, new double[] {1});
            Assert.AreEqual(target._constant, 1);
        }

        /// <summary>
        ///Ein Test für "Polynomial-Konstruktor"
        ///</summary>
        [TestMethod()]
        public void PolynomialConstructorTest2()
        {
            int order = 3;
            Polynomial_Accessor target = new Polynomial_Accessor(order);
            Assert.AreEqual(target._order, order);
            CollectionAssert.AreEqual(target._coeffs, new double[] { 0, 0, 0, 0 });
            Assert.AreEqual(target._constant, -1);
            // da alle zero-initialized wurden, sollte isConst = True liefern
            Assert.IsTrue(target.isConst());
        }

        #endregion

        #region Operator Tests
        /// <summary>
        ///Ein Test für den []-Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorSquareBracketsTest()
        {
            Polynomial target = new Polynomial(3, new double[] {3,2,1,0});
            double expected = target[2];
            double actual = target.getValueAt(2);
            Assert.AreEqual(expected, actual);
            Assert.AreEqual(expected, (double)1);
            Assert.AreEqual(actual, (double)1);
            Assert.AreNotEqual(expected, (double)2);
        }

        /// <summary>
        ///Ein Test für den ==-Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorEqualsTest()
        {
            Polynomial target = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial expected = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial notExpected = new Polynomial(4, new double[] { 4, 3, 2, 1, 0 });
            Assert.AreEqual(expected.ToString(), target.ToString());
            Assert.IsTrue(target == expected);
            Assert.IsFalse(target != expected);

            Assert.AreNotEqual(notExpected.ToString(), target.ToString());
            Assert.IsFalse(target == notExpected);
            Assert.IsTrue(target != notExpected);
        }

        /// <summary>
        ///Ein Test für den < und <= und > und >= -Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorLtGtTest()
        {
            Polynomial target = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial expected = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial greater = new Polynomial(3, new double[] { 4, 2, 1, 0 });
            Polynomial muchGreater = new Polynomial(4, new double[] { 4, 3, 2, 1, 0 });

            Assert.IsTrue(target <= expected);
            Assert.IsFalse(target < expected);

            Assert.IsTrue(target >= expected);
            Assert.IsFalse(target > expected);

            Assert.IsTrue(greater >= expected);
            Assert.IsTrue(greater > expected);

            Assert.IsTrue(expected < greater);
            Assert.IsTrue(expected <= greater);

            Assert.IsTrue(muchGreater > greater);
            Assert.IsTrue(muchGreater >= greater);

            Assert.IsFalse(muchGreater < greater);
            Assert.IsFalse(muchGreater <= greater);

            Assert.IsFalse(greater > muchGreater);
            Assert.IsFalse(greater >= muchGreater);
            Assert.IsFalse(expected >= greater);

            Assert.IsTrue(greater < muchGreater);
            Assert.IsTrue(greater <= muchGreater);
            Assert.IsFalse(greater <= expected);


        }

        /// <summary>
        ///Ein Test für den +-Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorPlusTest()
        {
            Polynomial a1 = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial a2 = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial expected = new Polynomial(3, new double[] { 6, 4, 2, 0 });
            Assert.IsTrue((a1 + a2) == expected);
            a1 += a2;
            Assert.AreEqual(a1, expected);
            Polynomial a3 = new Polynomial(3, new double[] { 6, 0, 0, 0 });
            Polynomial expectedA3 = new Polynomial(3, new double[] { 12, 4, 2, 0 });
            Assert.IsTrue(a1 + a3 == expectedA3);
            a1 += a3;
            Assert.AreEqual(a1, expectedA3);

            Polynomial a4 = new Polynomial(2, new double[] { 3, 2, 1 });
            try
            {
                a1 += a4;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            Polynomial a5 = new Polynomial(2, new double[] { 1, 0, 0 });
            a5 += a4;

        }

        /// <summary>
        ///Ein Test für den unary - -Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorUnaryMinusTest()
        {
            Polynomial a1 = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial a2 = new Polynomial(3, new double[] { -3, -2, -1, 0 });
            Assert.IsTrue(-a1 == a2);
            Assert.AreEqual(-a1, a2);
        }

        /// <summary>
        ///Ein Test für den - -Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorMinusTest()
        {
            Polynomial a1 = new Polynomial(3, new double[] { 6, 4, 2, 0 });
            Polynomial a2 = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Polynomial expected = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Assert.IsTrue((a1 - a2) == expected);
            a1 -= a2;
            Assert.AreEqual(a1, expected);
            Polynomial a3 = new Polynomial(3, new double[] { 6, 0, 0, 0 });
            Polynomial expectedA3 = new Polynomial(3, new double[] { -3, 2, 1, 0 });
            Assert.IsTrue(a1 - a3 == expectedA3);
            Assert.IsTrue(a3 - a1 == expectedA3);
            a1 -= a3;
            Assert.AreEqual(a1, expectedA3);

            Polynomial a4 = new Polynomial(2, new double[] { 6, 0, 0 });
            try
            {
                a3 -= a4;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        ///Ein Test für den *-Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorMultiplyTest()
        {
            Polynomial a1 = new Polynomial(2, new double[] { 1, 2, 3 });
            Polynomial a2 = new Polynomial(2, new double[] { 2, 3, 4 });
            Polynomial expected = new Polynomial(2, new double[] { 2, 7, 16 });
            Assert.IsTrue((a1 * a2) == expected);
            a1 *= a2;
            Assert.AreEqual(a1, expected);

            Polynomial a3 = new Polynomial(2, new double[] { 2, 0, 0 });
            Polynomial expectedA3 = new Polynomial(2, new double[] { 4, 14, 32 });
            Assert.IsTrue(a1 * a3 == expectedA3);
            Assert.IsTrue(a3 * a1 == expectedA3);
            a1 *= a3;
            Assert.AreEqual(a1, expectedA3);

            Polynomial a4 = new Polynomial(2, new double[] { 4, 5, 6 });
            double d = 2;
            Polynomial expectedA4 = new Polynomial(2, new double[] { 8, 10, 12 });
            Assert.IsTrue((a4 * d ) == expectedA4);
            a4 *= d;
            Assert.AreEqual(a4, expectedA4);

            Polynomial a5 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            try
            {
                a5 *= a4; 
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        ///Ein Test für den /-Operator
        ///</summary>
        [TestMethod()]
        public void PolynomialOperatorDivTest()
        {
            Polynomial a1 = new Polynomial(2, new double[] { 2, 4, 6 });
            Polynomial a2 = new Polynomial(2, new double[] { 1, 2, 3 });

            Polynomial a3 = a1 / a2;
            Polynomial a4 = a2 * a3;

            Assert.AreEqual(a1, a4);

            a1 /= a2;

            Polynomial expected = new Polynomial(2, new double[] { 2, 0, 0 });
            Assert.AreEqual(a1, expected);
            Polynomial a5 = new Polynomial(3, new double[] {4, 3, 2, 1 });
            try
            {
                a5 /= a4;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

        }
        #endregion

        #region Function Tests

        /// <summary>
        ///Ein Test für "setSqrt()"
        ///</summary>
        [TestMethod()]
        public void PolynomialGetterTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Assert.AreEqual(P1.ncoeff, 4);
        }

        /// <summary>
        ///Ein Test für "sqr()"
        ///</summary>
        [TestMethod()]
        public void PolynomialSqrTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Polynomial P2 = P1.sqr();
            Polynomial PExpect = new Polynomial(3, new double[] { 1, 4, 10, 20 });
            Assert.AreEqual(P2, PExpect);
            Polynomial PNotExpect = new Polynomial(3, new double[] { 1, 4, 10, 30 });
            Assert.AreNotEqual(P2, PNotExpect);
        }

        /// <summary>
        ///Ein Test für "setSqr()"
        ///</summary>
        [TestMethod()]
        public void PolynomialSetSqrTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Polynomial P2 = P1.sqr();
            P1.setSqr();
            Assert.AreEqual(P2, P1);
        }

        /// <summary>
        ///Ein Test für "sqrt()"
        ///</summary>
        [TestMethod()]
        public void PolynomialSqrtTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 4, 10, 20 });
            Polynomial PExpect = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Polynomial P2 = P1.sqrt();
            Assert.AreEqual(P2, PExpect);
            Polynomial PNotExpect = new Polynomial(3, new double[] { 1, 2, 3, 5 });
            Assert.AreNotEqual(P2, PNotExpect);
        }

        /// <summary>
        ///Ein Test für "setSqrt()"
        ///</summary>
        [TestMethod()]
        public void PolynomialSetSqrtTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 4, 10, 20 });
            Polynomial P2 = P1.sqrt();
            P1.setSqrt();
            Assert.AreEqual(P2, P1);
        }

        /// <summary>
        ///Ein Test für "eval()"
        ///</summary>
        [TestMethod()]
        public void PolynomialEvalTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Assert.AreEqual(10, P1.eval(2,1));
            Assert.AreNotEqual(10, P1.eval(2, 0.5));
        }

        /// <summary>
        ///Ein Test für "feval()"
        ///</summary>
        [TestMethod()]
        public void PolynomialFEvalTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 4, 3, 2, 1});
            Assert.AreEqual(4, P1.feval());
            Assert.AreNotEqual(5, P1.feval());
        }

        /// <summary>
        ///Ein Test für "shift()"
        ///</summary>
        [TestMethod()]
        public void PolynomialShiftTest()
        {

            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Polynomial P2 = new Polynomial(3, new double[] { 2, 6, 12, 0 });
	        P1.shift();
	        Assert.AreEqual(P1, P2);
            Polynomial P3 = new Polynomial(3, new double[] { 6, 24, 0, 0 });
	        P1.shift();
            Assert.AreEqual(P1, P3);
        }

        /// <summary>
        ///Ein Test für "isConst"
        ///</summary>
        [TestMethod()]
        public void PolynomialIsConstTest()
        {
            Polynomial target = new Polynomial();
            Assert.IsTrue(target.isConst());
            Polynomial_Accessor target2 = new Polynomial_Accessor(2, new double[] { 0, 1, 2 });
            Assert.IsFalse(target2.isConst());
            target2.unsetConst();
            Assert.AreEqual(target2._constant, 0);
        }

        /// <summary>
        ///Ein Test für "isConst(eps)"
        ///</summary>
        [TestMethod()]
        public void PolynomialIsConstEpsTest()
        {
            Polynomial P1 = new Polynomial(2, new double[] { 4, 2, 3 });
            Polynomial P2 = new Polynomial(2, new double[] { 5, 3, 1 });

            Assert.IsTrue(P1.isConst(3));
            Assert.IsFalse(P2.isConst(2));
        }

        /// <summary>
        ///Ein Test für "isId()"
        ///</summary>
        [TestMethod()]
        public void PolynomialIsIdTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 0, 0, 0 });
            Polynomial P2 = new Polynomial(3, new double[] { 1, 2, 1, 2 });

            Assert.IsTrue(P1.isId());
            Assert.IsFalse(P2.isId());
        }

        /// <summary>
        ///Ein Test für "isId(eps)"
        ///</summary>
        [TestMethod()]
        public void PolynomialIsIdEpsTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });

            Assert.IsTrue(P1.isId(4));
            Assert.IsFalse(P1.isId(3));
        }

        /// <summary>
        ///Ein Test für "isZero()"
        ///</summary>
        [TestMethod()]
        public void PolynomialIsZeroTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 0, 0, 0, 0 });
            Polynomial P2 = new Polynomial(3, new double[] { 0, 0, 4, 0 });
            Polynomial P3 = new Polynomial(3, new double[] { 6, 0, 0, 0 });

            Assert.IsTrue(P1.isZero());
            Assert.IsFalse(P2.isZero());
            Assert.IsFalse(P3.isZero());
        }

        /// <summary>
        ///Ein Test für "isZero(eps)"
        ///</summary>
        [TestMethod()]
        public void PolynomialIsZeroEpsTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 2, 3, 4 });
            Assert.IsTrue(P1.isZero(4));
            Assert.IsFalse(P1.isZero(3));
        }
        
        /// <summary>
        ///Ein Test für "set2Zero() and set2zero(int)"
        ///</summary>
        [TestMethod()]
        public void PolynomialSet2zeroTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 1, 2, 1 });
            Assert.IsFalse(P1.isZero());
            Assert.IsFalse(P1.isZero(1));
            Assert.IsTrue(P1.isZero(2));
            P1.set2Zero();
            Assert.IsTrue(P1.isZero());
            P1 = new Polynomial(3, new double[] { 1, 1, 2, 1 });
            P1.set2Zero(2);
            Assert.IsTrue(P1.isZero(1));
        }

        /// <summary>
        ///Ein Test für "set2const()"
        ///</summary>
        [TestMethod()]
        public void PolynomialSet2ConstTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 1, 2, 1 });
            Assert.IsFalse(P1.isConst());
            P1.set2const(6);
            Assert.IsTrue(P1.isConst());
            Polynomial P2 = new Polynomial(3, new double[] { 6, 0, 0, 0 });
            Assert.AreEqual(P1, P2);
        }

        /// <summary>
        ///Ein Test für "setCoeffs(double[])"
        ///</summary>
        [TestMethod()]
        public void PolynomialSetCoeffsTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 1, 2, 1 });
            double[] d;
            try
            {
                d = new double[] { 1, 2, 3, 4, 5 };
                P1.setCoeffs(d);
                Assert.Fail(); // If it gets to this line, no exception was thrown
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            d = new double[] { 4, 2, 4, 2 };
            P1.setCoeffs(d);
            Polynomial P2 = new Polynomial(3, new double[] { 4, 2, 4, 2 });
            Assert.AreEqual(P1, P2);
            P1 = new Polynomial(3, new double[] { 1, 1, 2, 1 });
            d = new double[] { 4, 2, 4 };
            P2 = new Polynomial(3, new double[] { 4, 2, 4, 0 });
            P1.setCoeffs(d);
            Assert.AreEqual(P1, P2);            
        }

        /// <summary>
        ///Ein Test für "GetValue()"
        ///</summary>
        [TestMethod()]
        public void PolynomialGetValueTest()
        {
            Polynomial target = new Polynomial(3, new double[] { 3, 2, 1, 0 });
            Assert.IsTrue(target[2] == target.getValueAt(2));
            Assert.IsTrue(target.getValueAt(0) == (double)3);
            Assert.IsFalse(target.getValueAt(0) == target[2]);
            try
            {
                target.getValueAt(-1);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                target.getValueAt(6);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                double temp = target[-1];
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                target[-1] = 6;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                double temp = target[6];
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                target[6] = 6; 
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        ///Ein Test für "ToString"
        ///</summary>
        [TestMethod()]
        public void PolynomialToStringTest()
        {
            Polynomial target = new Polynomial(3);
            string expected = "0x^3\t + 0x^2\t + 0x^1\t + 0";
            string notExpected = "0x^3 + 0x^2 + 0x^1 + 0";
            string actual = "";
            actual = target.ToString();
            Assert.AreEqual(expected, actual);
            Assert.AreNotEqual(notExpected, actual);
            
            Polynomial target2 = new Polynomial();
            expected = "1";
            notExpected = "\t1";
            actual = target2.ToString();
            Assert.AreEqual(expected, actual);
            Assert.AreNotEqual(notExpected, actual);

        }

        /// <summary>
        ///Ein Test für "Equals and GetHashCode"
        ///</summary>
        [TestMethod()]
        public void PolynomialEqualsHashTest()
        {
            double[] d = new double[] { 1, 1, 2, 1 };
            Polynomial_Accessor P1 = new Polynomial_Accessor(3, d);
            Assert.IsFalse(P1.Equals(null));
            Assert.AreEqual(P1._coeffs.GetHashCode(), P1.GetHashCode());

        }

        /// <summary>
        ///Ein Test für "Print"
        ///</summary>
        [TestMethod()]
        public void PolynomialPrintTest()
        {
            Polynomial P1 = new Polynomial(3, new double[] { 1, 1, 2, 1 });

            string file = Path.GetTempFileName();
            FileStream fs = new FileStream(file, FileMode.Create);
            TextWriter tmp = Console.Out;
            StreamWriter sw = new StreamWriter(fs);
            Console.SetOut(sw);
            P1.print();
            Console.SetOut(tmp);
            sw.Close();
            string actual = File.ReadAllText(file);
            string expected = "1\t1\t2\t1\t" + System.Environment.NewLine;

            Assert.AreEqual(actual, expected);

            File.Delete(file);

            P1.print(file);
            actual = File.ReadAllText(file);
            Assert.AreEqual(actual, expected);
        }

        #endregion
    }
}
