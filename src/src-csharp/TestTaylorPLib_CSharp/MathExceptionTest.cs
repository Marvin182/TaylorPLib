using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using LibMatrix;

namespace TestTaylorPLib_CSharp
{
    /// <summary>
    ///Dies ist eine Testklasse für "MathException" und soll
    ///alle MathExceptionTest Komponententests enthalten.
    ///</summary>
    [TestClass]
    public class MathExceptionTest
    {
        [TestMethod]
        public void MathExceptionConstructorTest()
        {
            MathException me1 = new MathException();
            MathException me2 = new MathException();
            Assert.AreEqual(me1.Message, me2.Message);

            MathException me3 = new MathException("Exception of Type MathException...");
            Assert.AreEqual(me3.Message, "Exception of Type MathException...");
            Assert.AreEqual(me3.what(), "Exception of Type MathException...");
        }
    }
}
