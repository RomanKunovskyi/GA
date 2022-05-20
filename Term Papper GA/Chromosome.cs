using System.Collections;

namespace Term_Papper_GA
{
    public class Chromosome : ICloneable
    {
        public BitArray Array { get; private set; }
        public double Number { get; private set; }
        public Chromosome(BitArray array, int a, int b)
        {
            Array = array;
            Number = ToDouble(a, b);
        }
        public Chromosome(int chromosomeLength, int a, int b)
        {
            var rand = new Random();
            Array = new BitArray(chromosomeLength);
            for (int i = 0; i < chromosomeLength; i++)
            {
                Array[i] = rand.Next(2) == 1;
            }
            Number = ToDouble(a, b);
        }

        private Chromosome(BitArray array, double number)
        {
            Array = (BitArray)array.Clone();
            Number = number;
        }
        private double ToDouble(int a, int b)
        {
            double sum = 0;
            for (int i = 0; i < Array.Length; i++)
            {
                sum += Convert.ToInt32(Array[Array.Length - i - 1]) * Math.Pow(2, i);
            }
            return a + sum * ((b - a) / (Math.Pow(2, Array.Length) - 1));
        }
        public override string ToString()
        {
            string str = "";
            foreach (bool bit in Array)
            {
                if (bit)
                {
                    str += "1";
                }
                else
                {
                    str += "0";
                }
            }
            if (Number < 0)
            {
                return str + $" ={Number:f3}";
            }
            else
            {
                return str + $" = {Number:f3}";
            }
        }

        public object Clone()
        {
            return new Chromosome(Array, Number);
        }
    }
}
