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
    }
}
