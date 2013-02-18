using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LibMatrix
{
    /// <summary>
    /// Exception Class for own Exceptions in Matrix and Polynomial
    /// </summary>
    public class MathException : System.Exception
    {
        /// <summary>
        /// Base Constructor
        /// </summary>
        public MathException() : base() { }

        /// <summary>
        /// Base Constructor with initializing the Message
        /// </summary>
        /// <param name="message">Message the Exception should throw</param>
        public MathException(string message) : base(message) { }

        /// <summary>
        /// Returns the Exception description 
        /// </summary>
        /// <returns>String of the exception description</returns>
        public String what()
        {
            return Message;
        }
    }
}
