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
            #region definitions

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

            #endregion

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
            #region definitions

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

            #endregion

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
            #region definitions

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

            #endregion 

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

        ///<summary>
        ///Ein Test für mmCaAUTBPbC(double alpha, double beta, Matrix A, Matrix B, int[] piv)
        ///</summary>
        [TestMethod()]
        public void MatrixFunction_mmCaAUTBPbC()
        {
            #region definitions
            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            Pb[0, 1] = new Polynomial(1, new double[] { 2, 2 });

            Pb[1, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] PbUT = new Polynomial[5, 5];
            PbUT[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            PbUT[0, 1] = new Polynomial(1, new double[] { 1, 1 });

            PbUT[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            PbUT[1, 1] = new Polynomial(1, new double[] { 2, 2 });
            #endregion

            Matrix A = new Matrix(2, 2, Pa);
            Matrix B = new Matrix(2, 2, Pb);
            Matrix BUT = new Matrix(2,2,PbUT);
            Matrix C = new Matrix(2, 2, 1);
            Matrix D = new Matrix(2, 2, 1);
            int[] piv = { 1, 0 };
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * BUT * alpha) + (C * beta);
            D.mmCaABbC(alpha, beta, A, BUT);
            C.mmCaAUTBPbC(alpha, beta, A, B, piv);
            Assert.AreEqual(C.ToString(), expect.ToString());
            Assert.AreEqual(D.ToString(), expect.ToString());

            try
            {
                Matrix E = new Matrix(2, 1, 1);
                C.mmCaAUTBPbC(alpha, beta, E, B, piv);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix E = new Matrix(1, 1, 1);
                E.mmCaAUTBPbC(alpha, beta, A, B, piv);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaAATbC(double alpha, double beta, Matrix A) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaAATbC()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            P_A[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] P_C = new Polynomial[2, 2];
            P_C[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[0, 1] = new Polynomial(1, new double[] { 1, 1 });

            P_C[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            #endregion

            Matrix A = new Matrix(2, 2, P_A);
            Matrix AT = A.asTranspose();
            Matrix C = new Matrix(2, 2, P_C);
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * AT * alpha) + (C * beta);
            C.mmCaAATbC(alpha, beta, A);

            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                Matrix D = new Matrix(1, 1, 1);
                D.mmCaAATbC(alpha, beta, A);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaATAbC(double alpha, double beta, Matrix A) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaATAbC()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            P_A[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] P_C = new Polynomial[2, 2];
            P_C[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[0, 1] = new Polynomial(1, new double[] { 1, 1 });

            P_C[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            #endregion

            Matrix A = new Matrix(2, 2, P_A);
            Matrix AT = A.asTranspose();
            Matrix C = new Matrix(2, 2, P_C);
            double alpha = 2;
            double beta = 2;

            Matrix expect = ( AT * A * alpha) + (C * beta);
            C.mmCaATAbC(alpha, beta, A);

            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                Matrix D = new Matrix(1, 1, 1);
                D.mmCaATAbC(alpha, beta, A);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaATBbC(double alpha, double beta, Matrix A, Matrix B) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaATBbC()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            P_A[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });

            Polynomial[,] P_C = new Polynomial[2, 2];
            P_C[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[0, 1] = new Polynomial(1, new double[] { 1, 1 });

            P_C[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            #endregion

            Matrix A = new Matrix(2, 2, P_A);
            Matrix B = new Matrix(2, 2, P_B);
            Matrix C = new Matrix(2, 2, P_C);
            Matrix AT = A.asTranspose();
            double alpha = 2;
            double beta = 2;

            Matrix expect = (AT * B * alpha) + (C * beta);
            C.mmCaATBbC(alpha, beta, A, B);

            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                Matrix D = new Matrix(1, 1, 1);
                D.mmCaATBbC(alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                A = new Matrix(1, 2, 1);
                C.mmCaATBbC(alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaATBPbC(double alpha, double beta, Matrix A, Matrix B, int[] piv) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaATBPbC()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 3];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_A[0, 2] = new Polynomial(1, new double[] { 4, 3 });

            P_A[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 2 });
            P_A[1, 2] = new Polynomial(1, new double[] { 4, 3 });

            Polynomial[,] P_C = new Polynomial[3, 3];
            P_C[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[0, 1] = new Polynomial(1, new double[] { 1, 1 });
            P_C[0, 2] = new Polynomial(1, new double[] { 2, 2 });

            P_C[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[1, 1] = new Polynomial(1, new double[] { 1, 1 });
            P_C[1, 2] = new Polynomial(1, new double[] { 2, 2 });

            P_C[2, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_C[2, 1] = new Polynomial(1, new double[] { 1, 1 });
            P_C[2, 2] = new Polynomial(1, new double[] { 2, 2 });

            Polynomial[,] P_B = new Polynomial[2, 3];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });
            P_B[0, 2] = new Polynomial(1, new double[] { 4, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });
            P_B[1, 2] = new Polynomial(1, new double[] { 4, 3 });

            Polynomial[,] P_B2 = new Polynomial[2, 3];
            P_B2[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B2[0, 1] = new Polynomial(1, new double[] { 4, 3 });
            P_B2[0, 2] = new Polynomial(1, new double[] { 2, 3 });

            P_B2[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B2[1, 1] = new Polynomial(1, new double[] { 4, 3 });
            P_B2[1, 2] = new Polynomial(1, new double[] { 3, 2 });

            #endregion

            Matrix A = new Matrix(2, 3, P_A);
            Matrix B = new Matrix(2, 3, P_B);
            Matrix B2 = new Matrix(2, 3, P_B2);
            Matrix C = new Matrix(3, 3, P_C);
            Matrix AT = A.asTranspose();

            double alpha = 2;
            double beta = 2;
            int[] piv = { 0, 2, 1 };

            Matrix expect = (AT * B * alpha) + (C * beta);
            C.mmCaATBPbC(alpha, beta, A, B2, piv);

            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                Matrix D = new Matrix(3, 3, 1);
                C.mmCaATBPbC(alpha, beta, D, B, piv);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                C = new Matrix(2, 2, 1);
                C.mmCaATBPbC(alpha, beta, A, B, piv);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaABTbC(double alpha, double beta, Matrix A, Matrix B) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaABTbC()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            P_A[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });

            #endregion

            Matrix A = new Matrix(2, 2, P_A);
            Matrix B = new Matrix(2, 2, P_B);
            Matrix C = new Matrix(2, 2, 1);
            Matrix BT = B.asTranspose();
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * BT * alpha) + (C * beta);
            Assert.AreNotEqual(expect.ToString(), C.ToString());
            C.mmCaABTbC(alpha, beta, A, B);
            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                Matrix D = new Matrix(3, 3, 1);
                C.mmCaABTbC(alpha, beta, A, D);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix D = new Matrix(3, 3, 1);
                D.mmCaABTbC(alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaABTbC(int r, bool up, double alpha, double beta, Matrix A, Matrix B)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaABTbC_2()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            P_A[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 0, 0 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 0, 0 });

            Polynomial[,] P_B2 = new Polynomial[2, 2];
            P_B2[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_B2[0, 1] = new Polynomial(1, new double[] { 3, 2 });

            P_B2[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_B2[1, 1] = new Polynomial(1, new double[] { 2, 3 });

            #endregion

            Matrix A = new Matrix(2, 2, P_A);
            Matrix B = new Matrix(2, 2, P_B);
            Matrix B2 = new Matrix(2, 2, P_B2);
            Matrix C = new Matrix(2, 2, 1);

            Matrix BT = B.asTranspose();
            Matrix BT2 = B2.asTranspose();
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * BT * alpha) + (C * beta);
            Assert.AreNotEqual(expect.ToString(), C.ToString());
            C.mmCaABTbC(1, true, alpha, beta, A, B);
            Assert.AreEqual(expect.ToString(), C.ToString());

            C = new Matrix(2, 2, 1);
            expect = (A * BT2 * alpha) + (C * beta);
            Assert.AreNotEqual(expect.ToString(), C.ToString());
            C.mmCaABTbC(1, false, alpha, beta, A, B2);
            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                A = new Matrix(2, 1, 1);
                B = new Matrix(2, 2, 1);
                C.mmCaABTbC(1, true, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                A = new Matrix(2, 2, 1);
                B = new Matrix(2, 2, 1);
                C.mmCaABTbC(3, true, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                A = new Matrix(3, 2, 1);
                B = new Matrix(2, 2, 1);
                C.mmCaABTbC(1, true, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für bmmCaABTbC(int r, int c, double alpha, double beta, Matrix A, Matrix B) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_bmmCaABTbC()
        {

            #region definitions

            Polynomial[,] P_A = new Polynomial[5, 5];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_A[0, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_A[0, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_A[0, 3] = new Polynomial(1, new double[] { 0, 0 });
            P_A[0, 4] = new Polynomial(1, new double[] { 0, 0 });

            P_A[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_A[1, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_A[1, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_A[1, 3] = new Polynomial(1, new double[] { 0, 0 });
            P_A[1, 4] = new Polynomial(1, new double[] { 0, 0 });

            P_A[2, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_A[2, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_A[2, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_A[2, 3] = new Polynomial(1, new double[] { 0, 0 });
            P_A[2, 4] = new Polynomial(1, new double[] { 0, 0 });

            P_A[3, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_A[3, 1] = new Polynomial(1, new double[] { 0, 0 });
            P_A[3, 2] = new Polynomial(1, new double[] { 0, 0 });
            P_A[3, 3] = new Polynomial(1, new double[] { 1, 0 });
            P_A[3, 4] = new Polynomial(1, new double[] { 0, 0 });

            P_A[4, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_A[4, 1] = new Polynomial(1, new double[] { 0, 0 });
            P_A[4, 2] = new Polynomial(1, new double[] { 0, 0 });
            P_A[4, 3] = new Polynomial(1, new double[] { 0, 0 });
            P_A[4, 4] = new Polynomial(1, new double[] { 1, 0 });

            Polynomial[,] P_B = new Polynomial[5, 5];
            P_B[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_B[0, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_B[0, 3] = new Polynomial(1, new double[] { 4, 1 });
            P_B[0, 4] = new Polynomial(1, new double[] { 5, 1 });

            P_B[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_B[1, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_B[1, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_B[1, 3] = new Polynomial(1, new double[] { 4, 1 });
            P_B[1, 4] = new Polynomial(1, new double[] { 5, 1 });
            
            P_B[2, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_B[2, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_B[2, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_B[2, 3] = new Polynomial(1, new double[] { 4, 1 });
            P_B[2, 4] = new Polynomial(1, new double[] { 5, 1 });
            
            P_B[3, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_B[3, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_B[3, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_B[3, 3] = new Polynomial(1, new double[] { 4, 1 });
            P_B[3, 4] = new Polynomial(1, new double[] { 5, 1 });

            P_B[4, 0] = new Polynomial(1, new double[] { 1, 1 });
            P_B[4, 1] = new Polynomial(1, new double[] { 2, 1 });
            P_B[4, 2] = new Polynomial(1, new double[] { 3, 1 });
            P_B[4, 3] = new Polynomial(1, new double[] { 4, 1 });
            P_B[4, 4] = new Polynomial(1, new double[] { 5, 1 });

            #endregion

            Matrix A = new Matrix(5, 5, P_A);
            Matrix B = new Matrix(5, 5, P_B);
            Matrix BT = B.asTranspose();
            Matrix C = new Matrix(5, 5, 1);
            double alpha = 2;
            double beta = 2;

            Matrix expect = (A * BT * alpha) + (C * beta);
            Assert.AreNotEqual(expect.ToString(), C.ToString());
            C.bmmCaABTbC(3, 3, alpha, beta, A, B);
            Assert.AreEqual(expect.ToString(), C.ToString());

            try
            {
                C.bmmCaABTbC(7, 3, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                C.bmmCaABTbC(3, 7, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                B = new Matrix(5, 4, 1);
                C.bmmCaABTbC(3, 3, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                A = new Matrix(6, 6, 1);
                B = new Matrix(6, 6, 1);
                C.bmmCaABTbC(3, 3, alpha, beta, A, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaIBbC(double alpha, double beta, Matrix B) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaIBbC()
        {
            #region definitions 

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 0 });
            P_A[0, 1] = new Polynomial(1, new double[] { 0, 0 });

            P_A[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 0 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });

            #endregion

            Matrix I = new Matrix(2, 2, P_A);
            Matrix B = new Matrix(2, 2, P_B);
            Matrix C = new Matrix(2, 2, 1);
            double alpha = 2;
            double beta = 2;

            Matrix expected = (I * B * alpha) + (C * beta);
            Assert.AreNotEqual(expected.ToString(), C.ToString());
            C.mmCaIBbC(alpha, beta, B);
            Assert.AreEqual(expected.ToString(), C.ToString());

            try
            {
                B = new Matrix(3, 3, 1);
                C.mmCaIBbC(alpha, beta, B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaIBbC(double alpha, double beta, int[] piv, bool rows, Matrix B)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaIBbC_2()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_A[0, 1] = new Polynomial(1, new double[] { 1, 0 });

            P_A[1, 0] = new Polynomial(1, new double[] { 1, 0 });
            P_A[1, 1] = new Polynomial(1, new double[] { 0, 0 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });

            #endregion

            Matrix I = new Matrix(2, 2, P_A);
            Matrix B = new Matrix(2, 2, P_B);
            Matrix C = new Matrix(2, 2, 1);
            double alpha = 2;
            double beta = 2;

            Matrix expected = (I * B * alpha) + (C * beta);
            Assert.AreNotEqual(expected.ToString(), C.ToString());
            C.mmCaIBbC(alpha, beta, new int[] { 1, 0 }, true, B);
            Assert.AreEqual(expected.ToString(), C.ToString());

            C = new Matrix(2, 2, 1);
            Assert.AreNotEqual(expected.ToString(), C.ToString());
            C.mmCaIBbC(alpha, beta, new int[] { 1, 0 }, false, B);
            Assert.AreEqual(expected.ToString(), C.ToString());
        }

        /// <summary>
        /// Ein Test für mmCaAIbC(double alpha, double beta, Matrix A) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaAIbC()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 1, 0 });
            P_A[0, 1] = new Polynomial(1, new double[] { 0, 0 });

            P_A[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_A[1, 1] = new Polynomial(1, new double[] { 1, 0 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });

            #endregion

            Matrix I = new Matrix(2, 2, P_A);
            Matrix A = new Matrix(2, 2, P_B);
            Matrix C = new Matrix(2, 2, 1);
            double alpha = 2;
            double beta = 2;

            Matrix expected = (A * I * alpha) + (C * beta);
            Assert.AreNotEqual(expected.ToString(), C.ToString());
            C.mmCaAIbC(alpha, beta, A);
            Assert.AreEqual(expected.ToString(), C.ToString());

            try
            {
                A = new Matrix(3, 3, 1);
                C.mmCaAIbC(alpha, beta, A);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für mmCaAIbC()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_mmCaAIbC_2()
        {
            #region definitions

            Polynomial[,] P_A = new Polynomial[2, 2];
            P_A[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            P_A[0, 1] = new Polynomial(1, new double[] { 1, 0 });

            P_A[1, 0] = new Polynomial(1, new double[] { 1, 0 });
            P_A[1, 1] = new Polynomial(1, new double[] { 0, 0 });

            Polynomial[,] P_B = new Polynomial[2, 2];
            P_B[0, 0] = new Polynomial(1, new double[] { 3, 2 });
            P_B[0, 1] = new Polynomial(1, new double[] { 2, 3 });

            P_B[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            P_B[1, 1] = new Polynomial(1, new double[] { 3, 2 });

            #endregion

            Matrix I = new Matrix(2, 2, P_A);
            Matrix A = new Matrix(2, 2, P_B);
            Matrix C = new Matrix(2, 2, 1);
            double alpha = 2;
            double beta = 2;

            Matrix expected = (A * I * alpha) + (C * beta);
            Assert.AreNotEqual(expected.ToString(), C.ToString());
            C.mmCaAIbC(alpha, beta, A, new int[] { 1, 0 }, true);
            Assert.AreEqual(expected.ToString(), C.ToString());

            C = new Matrix(2, 2, 1);
            Assert.AreNotEqual(expected.ToString(), C.ToString());
            C.mmCaAIbC(alpha, beta, A, new int[] { 1, 0 }, false);
            Assert.AreEqual(expected.ToString(), C.ToString());
        }

        /// <summary>
        /// Ein Test für utsolve(Matrix A) 
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_utsolve()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[3, 3];
            Pa[0, 0] = new Polynomial(1, new double[] { 3, 0 });
            Pa[0, 1] = new Polynomial(1, new double[] { 1, 0 });
            Pa[0, 2] = new Polynomial(1, new double[] { 0, 0 });

            Pa[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pa[1, 1] = new Polynomial(1, new double[] { 4, 0 });
            Pa[1, 2] = new Polynomial(1, new double[] { 6, 0 });

            Pa[2, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pa[2, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pa[2, 2] = new Polynomial(1, new double[] { 2, 0 });

            Polynomial[,] Pb = new Polynomial[3, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 1] = new Polynomial(1, new double[] { 1, 2 });

            Pb[1, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pb[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Pb[2, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[2, 1] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix A = new Matrix(3, 3, Pa);
            Matrix X = new Matrix(3, 2, Pb);
            Matrix B = A * X;
            A.utsolve(B);
            Assert.AreEqual(B.ToString(), X.ToString());

            try
            {
                X.utsolve(A);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                Matrix C = new Matrix(2, 2, 1);
                C.utsolve(B);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            A = new Matrix(3, 3, Pa);
            X = new Matrix(3, 2, Pb);
            B = A * X;
            Matrix actual = new Matrix(3, 2, 1);

            // X = new Matrix(3, 2);
            A.utsolve(B, actual, new int[] { 0, 1, 2 });

            Assert.AreEqual(actual.ToString(), X.ToString());

            A = new Matrix(3, 3, Pa);

            Polynomial[,] x = new Polynomial[3,1];
            x[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            x[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            x[2, 0] = new Polynomial(1, new double[] { 2, 1 });

            Polynomial[] x2 = new Polynomial[3];
            x2[0] = new Polynomial(1, new double[] { 8, 4 });
            x2[1] = new Polynomial(1, new double[] { 20, 10 });
            x2[2] = new Polynomial(1, new double[] { 4, 2 });

            Matrix XTest = new Matrix(3, 1, x);

            Matrix BTest = A * XTest;

            A.utsolve(x2);
            Assert.AreEqual(x2[0].ToString(), x[0, 0].ToString());
            Assert.AreEqual(x2[1].ToString(), x[1, 0].ToString());
            Assert.AreEqual(x2[2].ToString(), x[2, 0].ToString());

            Matrix U = new Matrix(3, 3, Pa);
            Polynomial[,] XSolve = new Polynomial[2, 3];
            XSolve[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            XSolve[0, 1] = new Polynomial(1, new double[] { 2, 1 });
            XSolve[0, 2] = new Polynomial(1, new double[] { 2, 1 });

            XSolve[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            XSolve[1, 1] = new Polynomial(1, new double[] { 2, 1 });
            XSolve[1, 2] = new Polynomial(1, new double[] { 2, 1 });

            X = new Matrix(2, 3, XSolve);

            B = X * U;

            U.utxsolve(B);

            Assert.AreEqual(B.ToString(), X.ToString());

        }
        
        /// <summary>
        /// Ein Test für cpermutem(int[] piv, bool trans = false)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_cpermutem()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 1] = new Polynomial(1, new double[] { 1, 2 });

            Pb[1, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pb[1, 1] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix expect = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);

            actual.cpermutem(new int[] { 1, 0 });
            Assert.AreEqual(actual.ToString(), expect.ToString());

            actual = new Matrix(2, 2, Pa);
            actual.cpermutem(new int[] { 1, 0 }, true);
            Assert.AreEqual(actual.ToString(), expect.ToString());

            try
            {
                actual.cpermutem(new int[] { 7, 0 });
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                actual.cpermutem(new int[] { -7, 0 });
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für rpermutem(int[] piv)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_rpermutem()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 1] = new Polynomial(1, new double[] { 1, 2 });

            Pb[1, 0] = new Polynomial(1, new double[] { 1, 2 });
            Pb[1, 1] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix expect = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);

            actual.rpermutem(new int[] { 1, 0 });
            Assert.AreEqual(actual.ToString(), expect.ToString());

            try
            {
                actual.rpermutem(new int[] { 7, 0 });
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                actual.rpermutem(new int[] { -7, 0 });
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für transpose()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_transpose()
        {
            #region definitions 

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 1, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[0, 1] = new Polynomial(1, new double[] { 2, 2 });

            Pb[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            Pb[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            Polynomial[,] Pc = new Polynomial[3, 2];
            Pc[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            Pc[0, 1] = new Polynomial(1, new double[] { 2, 2 });

            Pc[1, 0] = new Polynomial(1, new double[] { 3, 3 });
            Pc[1, 1] = new Polynomial(1, new double[] { 4, 4 });

            Pc[2, 0] = new Polynomial(1, new double[] { 5, 5 });
            Pc[2, 1] = new Polynomial(1, new double[] { 6, 6 });

            Polynomial[,] Pd = new Polynomial[2, 3];
            Pd[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            Pd[0, 1] = new Polynomial(1, new double[] { 3, 3 });
            Pd[0, 2] = new Polynomial(1, new double[] { 5, 5 });

            Pd[1, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pd[1, 1] = new Polynomial(1, new double[] { 4, 4 });
            Pd[1, 2] = new Polynomial(1, new double[] { 6, 6 });

            #endregion

            Matrix expect = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);

            actual.transpose();
            Assert.AreEqual(expect.ToString(), actual.ToString());

            expect = new Matrix(3, 2, Pc);
            actual = new Matrix(2, 3, Pd);
            actual.transpose();

            Assert.AreEqual(expect.ToString(), actual.ToString());

        }

        /// <summary>
        /// Ein Test für asTranspose()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_asTranspose()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pa[0, 1] = new Polynomial(1, new double[] { 1, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 2 });
            Pb[0, 1] = new Polynomial(1, new double[] { 2, 2 });

            Pb[1, 0] = new Polynomial(1, new double[] { 1, 1 });
            Pb[1, 1] = new Polynomial(1, new double[] { 1, 1 });

            #endregion

            Matrix expect = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);

            Matrix c = actual.asTranspose();
            Assert.AreEqual(expect.ToString(), c.ToString());
            Assert.AreNotEqual(expect.ToString(), actual.ToString());

        }
        
        /// <summary>
        /// Ein Test für shift()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_shift()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 2, 3 });
            Pa[0, 1] = new Polynomial(1, new double[] { 1, 2 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 3 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 2 });

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 3, 0 });
            Pb[0, 1] = new Polynomial(1, new double[] { 2, 0 });

            Pb[1, 0] = new Polynomial(1, new double[] { 3, 0 });
            Pb[1, 1] = new Polynomial(1, new double[] { 2, 0 });

            #endregion
            
            Matrix expect = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);

            actual.shift();
            Assert.AreEqual(expect.ToString(), actual.ToString());
        }

        /// <summary>
        /// Ein Test für isId()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_isId()
        {
            #region definitions

            Polynomial[,] Pa = new Polynomial[2, 2];
            Pa[0, 0] = new Polynomial(1, new double[] { 1, 0 });
            Pa[0, 1] = new Polynomial(1, new double[] { 0, 0 });

            Pa[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pa[1, 1] = new Polynomial(1, new double[] { 1, 0 });

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 1, 0 });
            Pb[0, 1] = new Polynomial(1, new double[] { 0, 1 });

            Pb[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[1, 1] = new Polynomial(1, new double[] { 1, 0 });

            Polynomial[,] Pc = new Polynomial[2, 2];
            Pc[0, 0] = new Polynomial(1, new double[] { 1, 1 });
            Pc[0, 1] = new Polynomial(1, new double[] { 0, 0 });

            Pc[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pc[1, 1] = new Polynomial(1, new double[] { 1, 0 });

            #endregion

            Matrix idMatrix = new Matrix(2, 2, Pa);
            Matrix nonIdMatrix = new Matrix(2, 2, Pb);
            Matrix nonIdMatrix2 = new Matrix(2, 2, Pc);
            Matrix nonIdMatrix3 = new Matrix(2, 3, 3);
            Assert.IsTrue(idMatrix.isId());
            Assert.IsFalse(nonIdMatrix.isId());
            Assert.IsFalse(nonIdMatrix2.isId());
            Assert.IsFalse(nonIdMatrix3.isId());
        }

        /// <summary>
        /// Ein Test für isZero()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_isZero()
        {
            Matrix isZeroMatrix = new Matrix(2, 2, 2);
            Assert.IsTrue(isZeroMatrix.isZero());
            isZeroMatrix[0, 0] = new Polynomial(2, new double[] { 0, 2, 4 });
            Assert.IsFalse(isZeroMatrix.isZero());
        }

        /// <summary>
        /// Ein Test für set2Id()
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_set2Id()
        {
            Matrix isZeroMatrix = new Matrix(2, 2, 2);
            Assert.IsFalse(isZeroMatrix.isId());
            isZeroMatrix.set2Id();
            Assert.IsTrue(isZeroMatrix.isId());

            try
            {
                isZeroMatrix = new Matrix(3, 2, 2);
                isZeroMatrix.set2Id();
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für set2IdFromIndices(int firstRow, int lastRow, int firstCol, int lastCol)
        /// und set2Id(int top, int bottom, int left, int right)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_set2IdFromIndices()
        {
            #region definitions

            Polynomial[,] Pb = new Polynomial[5, 5];
            Pb[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[0, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[1, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pb[1, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pb[1, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[1, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[2, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[2, 1] = new Polynomial(1, new double[] { 1, 0 });
            Pb[2, 2] = new Polynomial(1, new double[] { 0, 0 });
            Pb[2, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[2, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[3, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[3, 1] = new Polynomial(1, new double[] { 0, 0 });
            Pb[3, 2] = new Polynomial(1, new double[] { 1, 0 });
            Pb[3, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[3, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pb[4, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pb[4, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pb[4, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pb[4, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pb[4, 4] = new Polynomial(1, new double[] { 2, 1 });

            Polynomial[,] Pa = new Polynomial[5, 5];
            Pa[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[0, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pa[0, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pa[0, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[0, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pa[2, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[2, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pa[2, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pa[2, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[2, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pa[3, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[3, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pa[3, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pa[3, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[3, 4] = new Polynomial(1, new double[] { 2, 1 });

            Pa[4, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[4, 1] = new Polynomial(1, new double[] { 2, 1 });
            Pa[4, 2] = new Polynomial(1, new double[] { 2, 1 });
            Pa[4, 3] = new Polynomial(1, new double[] { 2, 1 });
            Pa[4, 4] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix expect = new Matrix(5, 5, Pb);
            Matrix actual = new Matrix(5, 5, Pa);

            Assert.AreNotEqual(expect.ToString(), actual.ToString());
            actual.set2IdFromIndices(2, 3, 1, 2);
            Assert.AreEqual(expect.ToString(), actual.ToString());

            actual = new Matrix(5, 5, Pa);
            Assert.AreNotEqual(expect.ToString(), actual.ToString());
            actual.set2Id(2, 1, 1, 2);
            Assert.AreEqual(expect.ToString(), actual.ToString());

            actual = new Matrix(5, 5, Pa);
            try
            {
                actual.set2IdFromIndices(-2, 3, 1, 2);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                actual.set2IdFromIndices( 2, 3, -1, 2);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                actual.set2IdFromIndices(2, -3, 1, 3);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        /// <summary>
        /// Ein Test für set2Zero()
        /// und set2Zero(int top, int bottom, int left, int right)
        /// und set2ZeroFromIndices(int firstRow, int lastRow, int firstCol, int lastCol)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_set2Zero()
        {
            #region definitions 

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[0, 1] = new Polynomial(1, new double[] { 0, 0 });

            Pb[1, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pb[1, 1] = new Polynomial(1, new double[] { 0, 0 });

            Polynomial[,] Pa = new Polynomial[5, 5];
            Pa[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 1] = new Polynomial(1, new double[] { 2, 1 });

            Polynomial[,] Pc = new Polynomial[5, 5];
            Pc[0, 0] = new Polynomial(1, new double[] { 0, 0 });
            Pc[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pc[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pc[1, 1] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix expected = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);
            Assert.AreNotEqual(expected.ToString(), actual.ToString());
            actual.set2Zero();
            Assert.AreEqual(expected.ToString(), actual.ToString());

            expected = new Matrix(2, 2, Pc);
            actual = new Matrix(2, 2, Pa);
            Assert.AreNotEqual(expected.ToString(), actual.ToString());
            actual.set2Zero(0, 1, 0, 1);
            Assert.AreEqual(expected.ToString(), actual.ToString());

            expected = new Matrix(2, 2, Pc);
            actual = new Matrix(2, 2, Pa);
            Assert.AreNotEqual(expected.ToString(), actual.ToString());
            actual.set2ZeroFromIndices(0, 0, 0, 0);
            Assert.AreEqual(expected.ToString(), actual.ToString());

        }

        /// <summary>
        /// Ein Test für set2Val(double v)
        /// und set2Val(int top, int bottom, int left, int right, double v)
        /// und set2ValFromIndices(int firstRow, int lastRow, int firstCol, int lastCol, double v)
        /// </summary>
        [TestMethod()]
        public void MatrixFunction_set2val()
        {
            #region definitions

            Polynomial[,] Pb = new Polynomial[2, 2];
            Pb[0, 0] = new Polynomial(1, new double[] { 5, 0 });
            Pb[0, 1] = new Polynomial(1, new double[] { 5, 0 });

            Pb[1, 0] = new Polynomial(1, new double[] { 5, 0 });
            Pb[1, 1] = new Polynomial(1, new double[] { 5, 0 });

            Polynomial[,] Pa = new Polynomial[5, 5];
            Pa[0, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pa[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pa[1, 1] = new Polynomial(1, new double[] { 2, 1 });

            Polynomial[,] Pc = new Polynomial[5, 5];
            Pc[0, 0] = new Polynomial(1, new double[] { 4, 0 });
            Pc[0, 1] = new Polynomial(1, new double[] { 2, 1 });

            Pc[1, 0] = new Polynomial(1, new double[] { 2, 1 });
            Pc[1, 1] = new Polynomial(1, new double[] { 2, 1 });

            #endregion

            Matrix expected = new Matrix(2, 2, Pb);
            Matrix actual = new Matrix(2, 2, Pa);
            Assert.AreNotEqual(expected.ToString(), actual.ToString());
            actual.set2Val(5);
            Assert.AreEqual(expected.ToString(), actual.ToString());

            
            expected = new Matrix(2, 2, Pc);
            actual = new Matrix(2, 2, Pa);
            Assert.AreNotEqual(expected.ToString(), actual.ToString());
            actual.set2Val(0, 1, 0, 1, 4);
            Assert.AreEqual(expected.ToString(), actual.ToString());

            expected = new Matrix(2, 2, Pc);
            actual = new Matrix(2, 2, Pa);
            Assert.AreNotEqual(expected.ToString(), actual.ToString());
            actual.set2ValFromIndices(0, 0, 0, 0, 4);
            Assert.AreEqual(expected.ToString(), actual.ToString());

            try
            {
                actual.set2ValFromIndices( -1, 0, 0, 0, 4);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }

            try
            {
                actual.set2ValFromIndices( 0, 0, -1, 0, 4);
                Assert.Fail();
            }
            catch (MathException) { }
            catch (Exception) { Assert.Fail(); }
        }

        #endregion

    }
}
