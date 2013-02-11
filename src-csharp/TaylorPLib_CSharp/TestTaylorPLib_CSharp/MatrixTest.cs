using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using LibMatrix;

namespace TestTaylorPLib_CSharp
{
    /// <summary>
    ///Dies ist eine Testklasse für "MatrixTest" und soll
    ///alle MatrixTest Komponententests enthalten.
    ///</summary>
    [TestClass]
    public class MatrixTest
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
        ///Ein Test für "Matrix-Konstruktor"
        ///</summary>
        [TestMethod()]
        public void MatrixConstructorTest()
        {
            Matrix_Accessor target = new Matrix_Accessor(1, 1, 0);
            Matrix_Accessor expected = new Matrix_Accessor();
            Assert.AreEqual(target._rows, expected._rows);
            Assert.AreEqual(target._cols, expected._cols);
            Assert.AreEqual(target._dimT, expected._dimT);
            CollectionAssert.AreEqual(target._data, expected._data);
        }

        /// <summary>
        ///Ein Test für "Matrix-Konstruktor"
        ///</summary>
        [TestMethod()]
        public void MatrixConstructorTest1()
        {
            Matrix_Accessor target = new Matrix_Accessor(2, 2);
            Matrix_Accessor expected = new Matrix_Accessor(2, 2, 0);
            Assert.AreEqual(target._rows, expected._rows);
            Assert.AreEqual(target._cols, expected._cols);
            Assert.AreEqual(target._dimT, expected._dimT);
            CollectionAssert.AreEqual(target._data, expected._data);
        }

        /// <summary>
        ///Ein Test für "Matrix-Konstruktor"
        ///</summary>
        [TestMethod()]
        public void MatrixConstructorTest2()
        {
            Matrix forCopy = new Matrix(2, 2);
            Matrix_Accessor target = new Matrix_Accessor(2, 2);
            Matrix_Accessor expected = new Matrix_Accessor(forCopy);
            Assert.AreEqual(target._rows, expected._rows);
            Assert.AreEqual(target._cols, expected._cols);
            Assert.AreEqual(target._dimT, expected._dimT);
            CollectionAssert.AreEqual(target._data, expected._data);

            Polynomial[,] p = new Polynomial[2, 2];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    p[i, j] = new Polynomial();

            Matrix_Accessor target2 = new Matrix_Accessor(2, 2, p);
            Assert.AreEqual(target2._rows, expected._rows);
            Assert.AreEqual(target2._cols, expected._cols);
            Assert.AreEqual(target2._dimT, expected._dimT);
        }
        
        #endregion

        #region Getter Tests

        /// <summary>
        ///Ein Test für "Getter"
        ///</summary>
        [TestMethod()]
        public void MatrixGetterTest()
        {
            Polynomial[,] p = new Polynomial[2, 3];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 3; j++)
                    p[i, j] = new Polynomial(i);

            Matrix_Accessor initM = new Matrix_Accessor(2, 3, p);
            Assert.AreEqual(initM.nrows(), 2);
            Assert.AreEqual(initM.ncols(), 3);
            CollectionAssert.AreEqual(initM._data, p);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(p[i, j], initM.get(i, j));

            Assert.AreEqual(initM.dimT(), 0);
            try
            {
                initM.get(-1, 1);
                Assert.Fail(); // If it gets to this line, no exception was thrown
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                initM.get( 1,-1);
                Assert.Fail(); // If it gets to this line, no exception was thrown
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); } 
            try
            {
                initM.get( 5, 1);
                Assert.Fail(); // If it gets to this line, no exception was thrown
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                initM.get( 1, 51);
                Assert.Fail(); // If it gets to this line, no exception was thrown
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        #endregion

        #region Operator Tests

        /// <summary>
        ///Ein Test für den []-Operator
        ///</summary>
        [TestMethod()]
        public void MatrixOperatorSquareBracketsTest()
        {
            Polynomial[,] P1 = new Polynomial[2,2];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    P1[i,j] = new Polynomial(i);

            Matrix m1 = new Matrix(2, 2, P1);
            Polynomial actual = m1[1, 1];
            Polynomial expected = new Polynomial(1);
            Assert.AreEqual(actual, expected);

            m1 = new Matrix(2, 2);
            m1[1, 1] = expected;
            Assert.AreEqual(m1.get(1, 1), expected);
            Assert.AreEqual(m1[1,1], expected);

            #region Exception Catching

            try
            {
                Polynomial p = m1[-1, 1];
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                Polynomial p = m1[ 1,-1];
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                Polynomial p = m1[ 6, 1];
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                Polynomial p = m1[ 1, 6];
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                m1[-1, 1] = new Polynomial();
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                m1[ 1,-1] = new Polynomial();
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                m1[ 6, 1] = new Polynomial();
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                m1[ 1, 6] = new Polynomial();
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            #endregion
        }

        /// <summary>
        ///Ein Test für den == und != -Operator
        ///</summary>
        [TestMethod()]
        public void MatrixOperatorEqualsTest()
        {
            Matrix target = new Matrix(2,2);
            Matrix expected = new Matrix(2, 2);
            Matrix notExpected = new Matrix(3,3);

            Assert.AreEqual(expected.ToString(), target.ToString());
            Assert.IsTrue(target == expected);
            Assert.IsFalse(target != expected);

            Assert.AreNotEqual(notExpected.ToString(), target.ToString());
            Assert.IsFalse(target == notExpected);
            Assert.IsTrue(target != notExpected);

            expected[1, 1] = new Polynomial(4);
            Assert.IsFalse(target == expected);
        }

        /// <summary>
        ///Ein Test für den + -Operator
        ///</summary>
        [TestMethod()]
        public void MatrixOperatorPlusTest()
        {
            Polynomial[,] P1 = new Polynomial[2, 2];
            Polynomial[,] P2 = new Polynomial[2, 2];
            Polynomial[,] P3 = new Polynomial[2, 2];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    P1[i, j] = new Polynomial(3);
                    P2[i, j] = new Polynomial(3);
                    P3[i, j] = P1[i, j] + P2[i, j];
                }

            Matrix a = new Matrix(2, 2, P1);
            Matrix b = new Matrix(2, 2, P2);

            Matrix c = a + b;
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    Assert.AreEqual(c[i,j], P3[i,j]);

            a += b;

            string aS = a.ToString();
            string cS = c.ToString();
            Assert.AreEqual(aS, cS);
            Assert.IsTrue(a == c);

            
            try
            {
                Matrix d = new Matrix(2, 4);
                d += a;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                Matrix d = new Matrix(4, 2);
                d += a;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        ///Ein Test für den - -Operator
        ///</summary>
        [TestMethod()]
        public void MatrixOperatorMinusTest()
        {
            Polynomial[,] P1 = new Polynomial[2, 2];
            Polynomial[,] P2 = new Polynomial[2, 2];
            Polynomial[,] P3 = new Polynomial[2, 2];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    P1[i, j] = new Polynomial(3);
                    P2[i, j] = new Polynomial(3);
                    P3[i, j] = P1[i, j] - P2[i, j];
                }

            Matrix a = new Matrix(2, 2, P1);
            Matrix b = new Matrix(2, 2, P2);

            Matrix c = a - b;
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    Assert.AreEqual(c[i, j], P3[i, j]);

            a -= b;

            string aS = a.ToString();
            string cS = c.ToString();
            Assert.AreEqual(aS, cS);
            Assert.IsTrue(a == c);


            try
            {
                Matrix d = new Matrix(2, 4);
                d -= a;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                Matrix d = new Matrix(4, 2);
                d -= a;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        ///Ein Test für den unary - -Operator
        ///</summary>
        [TestMethod()]
        public void MatrixOperatorUnaryMinusTest()
        {
            Polynomial[,] P1 = new Polynomial[1,1];
            Polynomial[,] P2 = new Polynomial[1,1];

            P1[0, 0] = new Polynomial(2, new double[] { 1, 2, 3 });
            P2[0, 0] = new Polynomial(2, new double[] { -1, -2, -3 });

            Matrix M1 = new Matrix(1, 1, P1);
            Matrix M2 = new Matrix(1, 1, P2);

            Matrix M3 = -M1;
            string m2S = M2.ToString();
            string m3S = M3.ToString();
            Assert.AreEqual(m2S, m3S);
            Assert.IsTrue(M2 == M3);
        }

        /// <summary>
        ///Ein Test für den unary - -Operator
        ///</summary>
        [TestMethod()]
        public void MatrixOperatorMultiplyTest()
        {

            Polynomial[,] P1 = new Polynomial[1, 1];
            Polynomial[,] P2 = new Polynomial[1, 1];

            P1[0, 0] = new Polynomial(2, new double[] { 1, 2, -3 });
            P2[0, 0] = new Polynomial(2, new double[] { 4, 8, -12 });

            Matrix target = new Matrix(1, 1, P1);
            Matrix expected = new Matrix(1, 1, P2);

            Assert.AreEqual((target * 4).ToString(), expected.ToString());

            target *= 4;

            Assert.AreEqual(target.ToString(), expected.ToString());

            Polynomial[,] P3 = new Polynomial[1, 2];
            Polynomial[,] P4 = new Polynomial[2, 1];
            Polynomial[,] P5 = new Polynomial[1, 1];

            P3[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P3[0, 1] = new Polynomial(1, new double[] { 1, 2 });
            P4[0, 0] = new Polynomial(1, new double[] { 2, 3 });
            P4[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P5[0, 0] = new Polynomial(1, new double[] { 4, 14 });

            Matrix M1 = new Matrix(1, 2, P3);
            Matrix M2 = new Matrix(2, 1, P4);

            Matrix expectedByMultiply = new Matrix(1, 1, P5);

            Assert.AreEqual((M1 * M2).ToString(), expectedByMultiply.ToString());
            M1 *= M2;
            Assert.AreEqual(M1.ToString(), expectedByMultiply.ToString());

            try
            {
                M1 = new Matrix(1, 2, P3);
                M1 *= expectedByMultiply;
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        #endregion

        #region Function Tests

        ///<summary>
        ///Ein Test für mmCaABbC(double alpha, double beta, Matrix A, Matrix B) 
        ///</summary>
        [TestMethod()]
        public void MatrixFunction_mmCaABbC()
        {
            Polynomial[,] p = new Polynomial[2, 2];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    p[i, j] = new Polynomial(2, new double[] { 1, 2, 3 });
                }

            Matrix A = new Matrix(2, 2, p);
            Matrix B = new Matrix(2, 2, p);
            Matrix C = new Matrix(2, 2, 2);
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * B * alpha) + (C * beta);
            C.mmCaABbC(alpha, beta, A, B);
            Assert.AreEqual(C.ToString(), expect.ToString());

            try
            {
                Matrix M1 = new Matrix(1, 2, 2);
                M1.mmCaABbC(alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix M1 = new Matrix(2, 1, 2);
                C.mmCaABbC(alpha, beta, M1, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        ///<summary>
        ///Ein Test für bmmCaABbC(int r, int c, double alpha, double beta, Matrix A, Matrix B)
        ///</summary>
        [TestMethod()]
        public void MatrixFunction_bmmCaABbC()
        {
            Polynomial[,] P1 = new Polynomial[5, 5];
            P1[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            P1[0, 1] = new Polynomial(1, new double[] { 2, 2 });
            P1[0, 2] = new Polynomial(1, new double[] { 2, 2 });
            P1[0, 3] = new Polynomial(1, new double[] { 0, 0 });
            P1[0, 4] = new Polynomial(1, new double[] { 0, 0 });

            P1[1, 0] = new Polynomial(1, new double[] { 2, 2 });
            P1[1, 1] = new Polynomial(1, new double[] { 2, 2 });
            P1[1, 2] = new Polynomial(1, new double[] { 2, 2 });
            P1[1, 3] = new Polynomial(1, new double[] { 0, 0 });
            P1[1, 4] = new Polynomial(1, new double[] { 0, 0 });
            
            P1[2, 0] = new Polynomial(1, new double[] { 2, 2 });
            P1[2, 1] = new Polynomial(1, new double[] { 2, 2 });
            P1[2, 2] = new Polynomial(1, new double[] { 2, 2 });
            P1[2, 3] = new Polynomial(1, new double[] { 0, 0 });
            P1[2, 4] = new Polynomial(1, new double[] { 0, 0 });

            P1[3, 0] = new Polynomial(1, new double[] { 0, 0 });
            P1[3, 1] = new Polynomial(1, new double[] { 0, 0 });
            P1[3, 2] = new Polynomial(1, new double[] { 0, 0 });
            P1[3, 3] = new Polynomial(1, new double[] { 0, 1 });
            P1[3, 4] = new Polynomial(1, new double[] { 0, 0 });

            P1[4, 0] = new Polynomial(1, new double[] { 0, 0 });
            P1[4, 1] = new Polynomial(1, new double[] { 0, 0 });
            P1[4, 2] = new Polynomial(1, new double[] { 0, 0 });
            P1[4, 3] = new Polynomial(1, new double[] { 0, 0 });
            P1[4, 4] = new Polynomial(1, new double[] { 0, 1 });

            Matrix A = new Matrix(5, 5, P1);
            Matrix B = new Matrix(5, 5, P1);
            Matrix C = new Matrix(5, 5, 1);
            Matrix D = new Matrix(5, 5, 1);
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * B * alpha) + (C * beta);
            D.mmCaABbC(alpha, beta, A, B);
            C.bmmCaABbC(3, 3, alpha, beta, A, B);
            Assert.AreEqual(C.ToString(), expect.ToString());

            try
            {
                C.bmmCaABbC(7, 3, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                C.bmmCaABbC(3, 7, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
            try
            {
                Matrix E = new Matrix(4, 5, 1);
                E.bmmCaABbC(3, 3, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix E = new Matrix(4, 5, 1);
                A.bmmCaABbC(3, 3, alpha, beta, A, E);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        ///<summary>
        ///Ein Test für mmCasABbC(int r, double alpha, double beta, Matrix A, Matrix B)
        ///</summary>
        [TestMethod()]
        public void MatrixFunction_mmCasABbC()
        {
            Polynomial[,] Pb = new Polynomial[5, 5];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[0, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pb[0, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pb[0, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pb[0, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pb[1, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[1, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pb[1, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pb[1, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pb[1, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pb[2, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[2, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pb[2, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pb[2, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pb[2, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pb[3, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[3, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pb[3, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pb[3, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pb[3, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pb[4, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[4, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pb[4, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pb[4, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pb[4, 4] = new Polynomial(1, new double[] { 0, 0 });

            Polynomial[,] Pa = new Polynomial[5, 5];
            Pa[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pa[0, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pa[0, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pa[0, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pa[0, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pa[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pa[1, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pa[1, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pa[1, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pa[1, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pa[2, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pa[2, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pa[2, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pa[2, 3] = new Polynomial(1, new double[] { 0, 0 });
            Pa[2, 4] = new Polynomial(1, new double[] { 0, 0 });

            Pa[3, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pa[3, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pa[3, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pa[3, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[3, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pa[4, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pa[4, 1] = new Polynomial(1, new double[] { 2, 2 });
            Pa[4, 2] = new Polynomial(1, new double[] { 2, 2 });
            Pa[4, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[4, 4] = new Polynomial(1, new double[] { 2, 1 });

            Matrix A = new Matrix(5, 5, Pa);
            Matrix B = new Matrix(5, 5, Pb);
            double alpha = 2;
            double beta = 2;
            Matrix C = new Matrix(5, 5, 1);
            Matrix D = new Matrix(5, 5, 1);

            Matrix expect = (A * B * alpha) + (C * beta);
            D.mmCaABbC(alpha, beta, A, B);
            C.mmCasABbC(2, alpha, beta, A, B);
            Assert.AreEqual(C.ToString(), expect.ToString());
            Assert.AreEqual(D.ToString(), expect.ToString());

            try
            {
                C.mmCasABbC(7, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix E = new Matrix(5, 4, 1);
                C.mmCasABbC(2, alpha, beta, E, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix E = new Matrix(4, 4, 1);
                E.mmCasABbC(2, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

        }

        ///<summary>
        ///Ein Test für mmCaAsBbC(int r, double alpha, double beta, Matrix A, Matrix B)
        ///</summary>
        [TestMethod()]
        public void MatrixFunction_mmCaAsBbC()
        {
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

            Matrix A = new Matrix(5, 5, Pa);
            Matrix B = new Matrix(5, 5, Pb);
            double alpha = 2;
            double beta = 2;
            Matrix C = new Matrix(5, 5, 1);
            Matrix D = new Matrix(5, 5, 1);

            Matrix expect = (A * B * alpha) + (C * beta);
            D.mmCaABbC(alpha, beta, A, B);
            C.mmCaAsBbC(2, alpha, beta, A, B);
            Assert.AreEqual(C.ToString(), expect.ToString());
            Assert.AreEqual(D.ToString(), expect.ToString());

            try
            {
                C.mmCaAsBbC(7, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix E = new Matrix(5, 4, 1);
                C.mmCaAsBbC(2, alpha, beta, E, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix E = new Matrix(4, 4, 1);
                E.mmCaAsBbC(2, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

        }

        #endregion

    }
}
